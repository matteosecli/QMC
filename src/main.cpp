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
void  mc_sampling(GeneralParams& , VMCparams&, vec&, vec&, System&);

// The variational wave function
double  wave_function(mat&, double, GeneralParams&);

// The local energy
double  local_energy(mat&, double, double, GeneralParams&, System&);

// prints to screen the results of the calculations
void  output(VMCparams&, vec&, vec&);


/*--- Begin of main program ---*/
int main(int argc, char* argv[])
{
    // Start timing
    //boost::timer t;

    // Initialize variables
    char *outfilename;
    vec cumulative_e, cumulative_e2;

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
    vP.alpha = 0.5;
    vP.beta = 0;
    cumulative_e = zeros<vec>(vmcParams.max_variations+1);
    cumulative_e2 = zeros<vec>(vmcParams.max_variations+1);

    // Build the required system
    System MySystem;
    MySystem = System(gP.number_particles, gP.dimension);
    MySystem.add_potential(new Harmonic_osc(gP));    // Harmonic Oscillator potential
    MySystem.add_potential(new eRepulsion(gP));      // Electron-Electron Repulsion
    /* To be generalized into (e.g.) MySystem = new Fermions(gP).
     * After generalization, substitute with
     *      MySystem->add_potential(new Harmonic_osc(gP))
     */

    //  Do the mc sampling
    mc_sampling(gP, vmcParams, cumulative_e, cumulative_e2, MySystem);

    // Print out results
    output(vmcParams, cumulative_e, cumulative_e2);
    cumulative_e.reset();
    cumulative_e2.reset();
    ofile.close();  // close output file

    // Print out elapsed time
    //double time = t.elapsed();
    //cout <<"Elapsed time:\t\t" << time <<"s." <<endl;

    return 0;
}


// Monte Carlo sampling with the Metropolis algorithm
void mc_sampling(GeneralParams& gP,VMCparams& vmcParams,
                 vec& cumulative_e, vec& cumulative_e2, System& InputSystem)
{
    int cycles, variate, accept, i, j;
    long idum;
    double wfnew, wfold, alpha, energy, energy2, delta_e;
    mat r_old, r_new;
    alpha = 0.5*gP.nuclear_charge;
    idum=-1;

    // allocate matrices which contain the position of the particles
    r_old = zeros<mat>(gP.number_particles, gP.dimension);
    r_new = zeros<mat>(gP.number_particles, gP.dimension);

    // loop over variational parameters
    for (variate=1; variate <= vmcParams.max_variations; variate++){

        // initialisations of variational parameters and energies
        alpha += 0.01;
        energy = energy2 = 0; accept =0; delta_e=0;

        //  initial trial position, note calling with alpha
        //  and in three dimensions
        for (i = 0; i < gP.number_particles; i++) {
            for ( j=0; j < gP.dimension; j++) {
                r_old(i,j) = vmcParams.step_length*(ran1(&idum)-0.5);
            }
        }
        wfold = wave_function(r_old, alpha, gP);

        // loop over monte carlo cycles
//        #pragma omp parallel for reduction(+: accept, energy, energy2) private(wfold, wfnew, delta_e) //schedule(dynamic,2)
        for (cycles = 1; cycles <= vmcParams.number_cycles+vmcParams.thermalization; cycles++){

            // new position
            for (i = 0; i < gP.number_particles; i++) {
                for ( j=0; j < gP.dimension; j++) {
                    r_new(i,j) = r_old(i,j)+vmcParams.step_length*(ran1(&idum)-0.5);
                }
            }
            wfnew = wave_function(r_new, alpha, gP);

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
                delta_e = local_energy(r_old, alpha, wfold, gP, InputSystem);
                // update energies
                energy += delta_e;
                energy2 += pow(delta_e,2);
            }
        }   // end of loop over MC trials

        cout << "variational parameter= " << alpha
        << " accepted steps= " << accept << endl;

        // update the energy average and its squared
        cumulative_e[variate] = energy/vmcParams.number_cycles;
        cumulative_e2[variate] = energy2/vmcParams.number_cycles;
    }    // end of loop over variational  steps

    r_old.reset(); // free memory
    r_new.reset(); // free memory
}   // end mc_sampling function


// Function to compute the squared wave function, simplest form
double  wave_function(mat& r, double alpha, GeneralParams & gP)
{
  int i, j;
  double wf, argument, r_single_particle;

  argument = wf = 0;
//#pragma omp parallel for private(j) shared(r,argument) reduction(+:r_single_particle) firstprivate(dimension)
//  #pragma omp parallel for firstprivate(j,dimension) lastprivate(i,number_particles) reduction(+: r_single_particle, argument) //schedule(dynamic,2)
  for (i = 0; i < gP.number_particles; i++) {
    r_single_particle = 0;
    for (j = 0; j < gP.dimension; j++) {
      r_single_particle  += pow(r(i,j),2);
    }
    argument += r_single_particle/2.0;
  }
  wf = exp(-argument*alpha) ;
  return wf;
}

// Function to calculate the local energy with num derivative
double  local_energy(mat& r, double alpha, double wfold, GeneralParams & gP, System & InputSystem)
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
      wfminus = wave_function(r_minus, alpha, gP);
      wfplus  = wave_function(r_plus, alpha, gP);
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



void output(VMCparams & vmcParams, vec & cumulative_e, vec & cumulative_e2)
{
  int i;
  double alpha, variance, error;
  //alpha = 0.5*charge;
  alpha = 0.5;
  for( i=1; i <= vmcParams.max_variations; i++){
    alpha += 0.01;
    variance = cumulative_e2(i)-pow(cumulative_e(i),2);
    error=sqrt(variance/vmcParams.number_cycles);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << alpha;
    ofile << setw(15) << setprecision(8) << cumulative_e(i);
    ofile << setw(15) << setprecision(8) << variance;
    ofile << setw(15) << setprecision(8) << error << endl;
//    fprintf(output_file, "%12.5E %12.5E %12.5E %12.5E \n", alpha,cumulative_e[i],variance, error );
  }
//  fclose (output_file);
}  // end of function output
