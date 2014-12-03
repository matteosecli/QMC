#ifndef STRUCTS_H
#define STRUCTS_H

#include <time.h>

typedef int seed_type;

namespace QMC2
{
    //! Struct used to initialize VMC parameters.       // WARNING: VMC and DMC are mixed here!
    struct VMCparams {
        int max_variations = 10; //!< The number of variations for each variational parameter.
        int number_cycles = 1000000; //!< The number of Monte Carlo cycles.
        int thermalization = 100; //!< The number of thermalization steps (so obscure!).
        double step_length = 2; //!< The step length used in Monte Carlo walker.
        double dt = 0.001; //!< Time step
        bool jF_active = false; //!< Disable Jastrow by default.
        bool imp_active = false; //!< Disable importance sampling by default.
    };

    //! Struct used to initialize the varational parameters.
    struct VariationalParams {
        double alpha = 0.987; //!< The spatial variational parameter.
        double beta = 0.398; //!< The Jastrow variational parameter.
        double alpha_step = 0.001;
        double beta_step = 0.001;
        double alpha_old;
        double beta_old;
    };

    //! Struct used to initialize general parameters.
    struct GeneralParams {
        int number_particles = 2; //!< The number of particles.
        int dimension = 2; //!< The dimension.
        double frequency = 1.0;
        int nuclear_charge = 1;
        seed_type random_seed = -time(NULL); //!< The random number generator's seed.
        double systemConstant = 1;
        int num_threads = 1; //!< The maximum number of threads
    };
}


#endif // STRUCTS_H
