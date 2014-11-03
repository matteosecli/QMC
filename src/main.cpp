/*          FYS3150@UiO - QMC Project
 *          Matteo Seclì, Autumn 2014
 *          <secli.matteo@gmail.com>
 *
 * This program performs a Variational Monte Carlo
 * calculation for Quantum Dots. It works both in 2
 * and 3 dimensions, but you have to use a number
 * of electrons equal to the degeneracy of the energy
 * level you want to consider.
 * The program is structured in an object-oriented way,
 * in order to perform further generalization after
 * this project. As an example, you can easily remove
 * electron-electron interaction (just comment one line)
 * or change the program to perform calculations over
 * an atom instead of an harmonic oscillator (again,
 * just one line to change).
 *
 * The main obecjt structure is thanks to the work of
 * Jørgen Høgberget, from whom I've borrowed many lines
 * of code. He's done a great job with his project.
 * You can check his work here:
 *  <https://github.com/jorgehog/QMC2>
 * and here:
 *  <http://urn.nb.no/URN:NBN:no-38645>
 */

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include <armadillo>
#include "Potential.h"
#include "System.h"
#include "Structs.h"
#include "AlphaHarmonicOscillator.h"
#include "Jastrow.h"
#include <time.h>
#include <omp.h>

using namespace  std;
using namespace arma;
using namespace QMC2;

// output file as global variable
ofstream ofile;

// the step length and its squared inverse for the second derivative
#define h 0.001
#define h2 1000000


/*--- Declaraton of functions ---*/

// Function to read in data from screen, note call by reference
void initialise(GeneralParams&, VariationalParams&, VMCparams&) ;

// The Mc sampling for the variational Monte Carlo
void  mc_sampling(GeneralParams& , VariationalParams&, VMCparams&, mat &, mat &, System&, Orbitals*, Jastrow*);

// The variational wave function
double  wave_function(mat&, VariationalParams&, GeneralParams&, Orbitals*, Jastrow*);

// The local energy
double  local_energy(mat&, double, GeneralParams&, VariationalParams&, System&, Orbitals*, Jastrow*);

// prints to screen the results of the calculations
void  output(GeneralParams&, VariationalParams&, VMCparams&, mat &, mat &, bool&);


/*--- Begin of main program ---*/
int main(int argc, char* argv[])
{
    // Start timing
    //boost::timer t;

    // Initialize variables
    char *outfilename;
    mat cumulative_e, cumulative_e2;

    // Read in output file, abort if there are too few command-line arguments
    if( argc <= 1 ){
        cout << "Bad Usage: " << argv[0] << " read also output file on same line" << endl;
        exit(1);
    }
    else{
        outfilename=argv[1];
    }
    ofile.open(outfilename);

    // Initialize structs
    struct GeneralParams gP;
    struct VariationalParams vP;
    struct VMCparams vmcParams;

    // Read in data and filling structs
    vP.beta = 0;    // Default value for No-Jastrow, can be overridden
    vP.beta_old = vP.beta;
    initialise(gP, vP, vmcParams);
    cumulative_e = zeros<mat>(vmcParams.max_variations+1,vmcParams.max_variations+1);
    cumulative_e2 = zeros<mat>(vmcParams.max_variations+1,vmcParams.max_variations+1);

    // Build the required system
    System MySystem;
    MySystem = System(gP.number_particles, gP.dimension);
    MySystem.add_potential(new Harmonic_osc(gP));    // Harmonic Oscillator potential
    MySystem.add_potential(new eRepulsion(gP));      // Electron-Electron Repulsion
    /* To be generalized into (e.g.) MySystem = new Fermions(gP).
     * After generalization, substitute with
     *      MySystem->add_potential(new Harmonic_osc(gP))
     */

    // Begin parallelization
    omp_set_num_threads(gP.num_threads);
    #pragma omp parallel shared(cumulative_e, cumulative_e2)
    {
        // Make a private copy of all the stuff
        mat cumulative_e_local = zeros<mat>(vmcParams.max_variations+1,vmcParams.max_variations+1);
        mat cumulative_e2_local = zeros<mat>(vmcParams.max_variations+1,vmcParams.max_variations+1);
        struct GeneralParams gP_local = gP;
        struct VariationalParams vP_local = vP;
        struct VMCparams vmcParams_local = vmcParams;

        // Build the orbitals
        Orbitals* MyWaveFunction_local = new AlphaHarmonicOscillator(gP_local, vP_local);

        // Set up the Jastrow factor
        Jastrow* MyJastrowFactor_local;
        if (vmcParams_local.jF_active) MyJastrowFactor_local = new Pade_Jastrow(gP,vP);
        else MyJastrowFactor_local = new No_Jastrow();

        // Make a local copy of the system for convenience
        System MySystem_local = MySystem;

        // Do the mc sampling
        mc_sampling(gP_local, vP_local, vmcParams_local,
                    cumulative_e_local, cumulative_e2_local,
                    MySystem_local, MyWaveFunction_local, MyJastrowFactor_local);

        // Be sure to have all the contributions
        #pragma omp barrier
        #pragma omp critical
        {
            // Add the contributions
            cumulative_e += cumulative_e_local;
            cumulative_e2 += cumulative_e2_local;
        }
    }

    // Normalize to the number of threads
    cumulative_e = cumulative_e/((double) gP.num_threads);
    cumulative_e2 = cumulative_e2/((double) gP.num_threads);

    // Print out results
    output(gP, vP, vmcParams, cumulative_e, cumulative_e2, vmcParams.jF_active);

    // Do some cleaning
    cumulative_e.reset();
    cumulative_e2.reset();
    ofile.close();  // close output file

    // Print out elapsed time
    //double time = t.elapsed();
    //cout <<"Elapsed time:\t\t" << time <<"s." <<endl;

    return 0;
}


// Monte Carlo sampling with the Metropolis algorithm
void mc_sampling(GeneralParams& gP, VariationalParams& vP, VMCparams& vmcParams,
                 mat& cumulative_e, mat& cumulative_e2, System& InputSystem, Orbitals* wF, Jastrow* jF)
{
    int cycles, variate_alpha, variate_beta, accept, i, j;
    long idum = -1;
    double wfnew, wfold, energy, energy2, delta_e;
    mat r_old, r_new;

    vP.alpha_old = vP.alpha;
    vP.beta_old = vP.beta;

    // allocate matrices which contain the position of the particles
    r_old = zeros<mat>(gP.number_particles, gP.dimension);
    r_new = zeros<mat>(gP.number_particles, gP.dimension);

    // loop over variational parameters
    for (variate_beta = 1; variate_beta <= vmcParams.max_variations; variate_beta++){
//        #pragma omp parallel for
        for (variate_alpha = 1; variate_alpha <= vmcParams.max_variations; variate_alpha++){

            // initialisations of energies
            energy = energy2 = 0; accept =0; delta_e=0;

            //  initial trial position, note calling with alpha
            //  and in three dimensions
            for (i = 0; i < gP.number_particles; i++) {
                for ( j=0; j < gP.dimension; j++) {
                    r_old(i,j) = vmcParams.step_length*(ran1(&idum)-0.5);
                }
            }

            wfold = wave_function(r_old, vP, gP, wF, jF);

            // loop over monte carlo cycles
            for (int cycles = 1; cycles <= vmcParams.number_cycles+vmcParams.thermalization; cycles++){

                // new position
                for (int i = 0; i < gP.number_particles; i++) {
                    for ( int j=0; j < gP.dimension; j++) {
                        r_new(i,j) = r_old(i,j)+vmcParams.step_length*(ran1(&idum)-0.5);
                    }
                }
                wfnew = wave_function(r_new, vP, gP, wF, jF);

                // Metropolis test
                if(ran1(&idum) <= wfnew*wfnew/wfold/wfold ) {
                    for (int i = 0; i < gP.number_particles; i++) {
                        for (int  j=0; j < gP.dimension; j++) {
                            r_old(i,j)=r_new(i,j);
                        }
                    }
                    wfold = wfnew;
                    accept = accept+1;
                }

                // compute local energy
                if ( cycles > vmcParams.thermalization ) {
                    delta_e = local_energy(r_old, wfold, gP, vP, InputSystem, wF, jF);
                    // update energies
                    energy += delta_e;
                    energy2 += pow(delta_e,2);
                }

            }   // end of loop over MC trials

            cout << "alpha = " << vP.alpha << "\t"
                 << "beta = " << vP.beta << "\t"
                 << " accepted steps = " << accept << endl;

            // update the energy average and its squared
            cumulative_e(variate_beta, variate_alpha) = energy/vmcParams.number_cycles;
            cumulative_e2(variate_beta, variate_alpha) = energy2/vmcParams.number_cycles;

            // update alpha
            vP.alpha += vP.alpha_step;

        }    // end of loop over variational steps alpha

        // Reset alpha
        vP.alpha = vP.alpha_old;

        // Update beta if Jastrow is active
        if (jF->active) {
            vP.beta += vP.beta_step;
        } else {
            break;
        }

    } // end of loop over variational steps beta

    // Reset beta
    vP.beta = vP.beta_old;

    r_old.reset(); // free memory
    r_new.reset(); // free memory
}   // end mc_sampling function


// Function to compute the squared wave function, simplest form
double  wave_function(mat& r, VariationalParams & vP, GeneralParams & gP, Orbitals *wF, Jastrow *jF)
{

    // Update the pointers
    wF->set_parameter(vP.alpha, 0);
    jF->set_parameter(vP.beta, 0);

    // Just call the Slater Determinant and the Jastow Factor calculators
    double wavefunction = wF->SlaterD(r) * jF->get_val(r);

    return wavefunction;
}

// Function to calculate the local energy with num derivative
double  local_energy(mat& r, double wfold, GeneralParams & gP, VariationalParams & vP, System & InputSystem, Orbitals *wF, Jastrow *jF)
{
  int i, j , k;
  double e_local, wfminus, wfplus, e_kinetic, e_potential, r_12,
    r_single_particle;
  mat r_plus, r_minus;

  // allocate matrices which contain the position of the particles
  // the function matrix is defined in the progam library
  r_plus = mat(gP.number_particles, gP.dimension);
  r_minus = mat(gP.number_particles, gP.dimension);
  r_plus = r_minus = r;
  // compute the kinetic energy
  e_kinetic = 0;
  for (i = 0; i < gP.number_particles; i++) {
    for (j = 0; j < gP.dimension; j++) {
      r_plus(i,j) = r(i,j)+h;
      r_minus(i,j) = r(i,j)-h;
      wfminus = wave_function(r_minus, vP, gP, wF, jF);
      wfplus  = wave_function(r_plus, vP, gP, wF, jF);
      e_kinetic -= (wfminus+wfplus-2*wfold);
      r_plus(i,j) = r(i,j);
      r_minus(i,j) = r(i,j);
    }
  }
  // include electron mass and hbar squared and divide by wave function
  e_kinetic = 0.5*h2*e_kinetic/wfold;

  // Add the potential energy from the system object
  e_potential = InputSystem.get_potential_energy(r);

  r_plus.reset(); // free memory
  r_minus.reset();
  e_local = e_potential+e_kinetic;
  return e_local;
}

void initialise(GeneralParams & gP, VariationalParams & vP, VMCparams & vmcParams)
{
    int jF_active_resp;

    cout << "Number of particles = ";
    cin >> gP.number_particles;

    cout << "Dimensionality = ";
    cin >> gP.dimension;
    gP.dimension = 2;   // Override user dimensionality, waiting for 3D implementation

    cout << "Number of threads (max " << omp_get_max_threads() << ") = ";
    cin >> gP.num_threads;

    cout << "Do you want to use a Jastrow factor? (1-y / 0-n) = ";
    cin >> jF_active_resp;
    if (jF_active_resp != 0) vmcParams.jF_active = true;

    cout << "Spatial variational parameter start point = ";
    cin >> vP.alpha;    //0.91;
    vP.alpha_old = vP.alpha;

    cout << "Spatial variational parameter step length = ";
    cin >> vP.alpha_step;

    if (jF_active_resp != 0) {

        cout << "Jastrow variational parameter start point = ";
        cin >> vP.beta;     //0.32;
        vP.beta_old = vP.beta;

        cout << "Jastrow variational parameter step length = ";
        cin >> vP.beta_step;
    }

    cout << "Maximum variations per parameter = ";
    cin >> vmcParams.max_variations;

    cout << "# Monte Carlo steps = ";
    cin >> vmcParams.number_cycles;

    cout << "# Step length = ";
    cin >> vmcParams.step_length;
}  // end of function initialise


// UPDATE OUTPUT FUNCTION! (RESET VAR PARAMS, THEN CALL THE FUNCTION)
void output(GeneralParams& gP, VariationalParams& vP, VMCparams& vmcParams, mat& cumulative_e, mat& cumulative_e2, bool& jF_active)
{
  int i, j;
  double variance, error;
  for( i=1; i <= vmcParams.max_variations; i++){
      for( j=1; j <= vmcParams.max_variations; j++){
        variance = (cumulative_e2(i,j)-pow(cumulative_e(i,j),2))/gP.num_threads;
        error=sqrt(variance/vmcParams.number_cycles);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << vP.alpha;
        ofile << setw(15) << setprecision(8) << vP.beta;
        ofile << setw(15) << setprecision(8) << cumulative_e(i,j);
        ofile << setw(15) << setprecision(8) << variance;
        ofile << setw(15) << setprecision(8) << error << endl;
        vP.alpha += vP.alpha_step;
      }
      vP.alpha = vP.alpha_old;
      if (!jF_active) break;
      vP.beta += vP.beta_step;
  }
//  fclose (output_file);
}  // end of function output
