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
void initialise(GeneralParams&, VMCparams&) ;

// The Mc sampling for the variational Monte Carlo
void  mc_sampling(GeneralParams& , VariationalParams&, VMCparams&, mat &, mat &, System&, Orbitals*, Jastrow*);

// The variational wave function
double  wave_function(mat&, VariationalParams&, GeneralParams&, Orbitals*, Jastrow*);

// The local energy
double  local_energy(mat&, double, GeneralParams&, VariationalParams&, System&, Orbitals*, Jastrow*);

// prints to screen the results of the calculations
void  output(VMCparams&, mat &, mat &, bool&);


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
    initialise(gP, vmcParams);
    vP.alpha = 0.91;
    vP.beta = 0.32;
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

    // Build the orbitals
    Orbitals* MyWaveFunction = new AlphaHarmonicOscillator(gP, vP);

    // Set up the Jastrow factor
    Jastrow* MyJastrowFactor = new Pade_Jastrow(gP,vP);
//    Jastrow* MyJastrowFactor = new No_Jastrow();

    //  Do the mc sampling
    mc_sampling(gP, vP, vmcParams, cumulative_e, cumulative_e2, MySystem, MyWaveFunction, MyJastrowFactor);

    // Print out results
    bool jF_active = MyJastrowFactor->active;
    output(vmcParams, cumulative_e, cumulative_e2, jF_active);
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
    double wfnew, wfold, alpha, alpha_old, energy, energy2, delta_e;
    mat r_old, r_new;

    alpha_old = vP.alpha;

    // allocate matrices which contain the position of the particles
    r_old = zeros<mat>(gP.number_particles, gP.dimension);
    r_new = zeros<mat>(gP.number_particles, gP.dimension);

    // loop over variational parameters
    for (variate_beta = 1; variate_beta <= vmcParams.max_variations; variate_beta++){
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
    //        #pragma omp parallel for reduction(+: accept, energy, energy2) private(wfold, wfnew, delta_e) //schedule(dynamic,2)
            for (cycles = 1; cycles <= vmcParams.number_cycles+vmcParams.thermalization; cycles++){

                // new position
                for (i = 0; i < gP.number_particles; i++) {
                    for ( j=0; j < gP.dimension; j++) {
                        r_new(i,j) = r_old(i,j)+vmcParams.step_length*(ran1(&idum)-0.5);
                    }
                }
                wfnew = wave_function(r_new, vP, gP, wF, jF);

                // Metropolis test
                if(ran1(&idum) <= wfnew*wfnew/wfold/wfold ) {
                    for (i = 0; i < gP.number_particles; i++) {
                        for ( j=0; j < gP.dimension; j++) {
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
            vP.alpha += 0.01;

        }    // end of loop over variational steps alpha

        // Reset alpha
        vP.alpha = alpha_old;

        // Update beta if Jastrow is active
        if (jF->active) {
            vP.beta += 0.01;
        } else {
            break;
        }

    } // end of loop over variational steps beta

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

void initialise(GeneralParams & gP, VMCparams & vmcParams)
{
    cout << "number of particles = ";
    cin >> gP.number_particles;
//    cout << "charge of nucleus = ";
//    cin >> gP.nuclear_charge;     // Default: 1
    cout << "dimensionality = ";
    cin >> gP.dimension;
    cout << "maximum variational parameters = ";
    cin >> vmcParams.max_variations;
//    cout << "# Thermalization  steps= ";
//    cin >> vmcParams.thermalization;     // Default: 100
    cout << "# MC steps= ";
    cin >> vmcParams.number_cycles;
    cout << "# step length= ";
    cin >> vmcParams.step_length;
}  // end of function initialise


// UPDATE OUTPUT FUNCTION! (RESET VAR PARAMS, THEN CALL THE FUNCTION)
void output(VMCparams & vmcParams, mat& cumulative_e, mat& cumulative_e2, bool& jF_active)
{

  int i, j;
  double alpha, alpha_old, beta, variance, error;
  //alpha = 0.5*charge;
  alpha = 0.91;
  alpha_old = alpha;
  beta = 0.32;
  for( i=1; i <= vmcParams.max_variations; i++){
      for( j=1; j <= vmcParams.max_variations; j++){
        variance = cumulative_e2(i,j)-pow(cumulative_e(i,j),2);
        error=sqrt(variance/vmcParams.number_cycles);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setw(15) << setprecision(8) << alpha;
        ofile << setw(15) << setprecision(8) << beta;
        ofile << setw(15) << setprecision(8) << cumulative_e(i,j);
        ofile << setw(15) << setprecision(8) << variance;
        ofile << setw(15) << setprecision(8) << error << endl;
        alpha += 0.01;
      }
      alpha = alpha_old;
      if (!jF_active) break;
      beta += 0.01;
  }
//  fclose (output_file);
}  // end of function output
