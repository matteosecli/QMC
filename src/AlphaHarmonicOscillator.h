#include <armadillo>
#include <string>

using namespace arma;

namespace QMC2
{

    class BasisFunctions;
    struct GeneralParams;
    struct VariationalParams;

    /*! \brief Superclass for the single particle orbital classes.
     * Handles everything specific regarding choice of single particle basis.
     */
    class Orbitals {
        protected:
            int n_p;
            int n2;
            int dim;

            int nCap; //!< The highest number of states accessible

            int max_implemented; //!< The maximum number basis size supported for any system ##RYDD OPP

            arma::imat qnums; //!< Quantum number matrix needed by Hartree-Fock and the variational derivatives.

            BasisFunctions** basis_functions; //!< A vector maping a quantum number index to a single particle wave function.

            std::string name;

            virtual double get_sp_energy(int qnum) const {
                (void) qnum;

                //Do nothing
                return 0;
            }

            void set_n_p(int n_p) {
                this->n_p = n_p;
            }

            void set_dim(int dim) {
                this->dim = dim;
            }


        public:

            Orbitals(GeneralParams &gP, VariationalParams &vP);
            Orbitals(int n_p, int dim);
            Orbitals();

            //! A method for retrieving variational parameters.
            /*!
             * @param n Index of the sought variational parameter.
             */
            virtual double get_parameter(int n) {
                (void) n;

                std::cout << "ATTEMPT TO GET PARAMETER FROM ASSUMINGLY PARAMETER FREE ORBITAL " << name << std::endl;

                return 1.0;
            }

            //! A method for setting variational parameters.
            /*!
             * @param n Index of the sought variational parameter.
             * @param parameter The new value of the variational parameter.
             */
            virtual void set_parameter(double parameter, int n) {
                (void) parameter;
                (void) n;

                std::cout << "ATTEMPT TO SET PARAMETER IN ASSUMINGLY PARAMETER FREE ORBITAL " << name << std::endl;
            }

            //! Calculates single particle wave function terms which are independent of the quantum numbers

            /*!
             * If a term in the single particle functions are independent of the quantum number,
             * this function can be overridden to calculate them beforehand (for each particle),
             * and rather extract the value instead of recalculating.
             * @param i Particle number.
             */
            virtual void set_qnum_indie_terms(mat& r, int i) {
                (void) r;
                (void) i;
            }

            //! Calculates the single particle wave function for a given walker's particle.
            /*!
             * @param q_num The quantum number index.
             */
            virtual double phi(const mat& r, int particle, int q_num);

            std::string getName() const {
                return name;
            }

    };


//    struct GeneralParams;
//    struct VariationalParams;

    /*! \brief Harmonic Oscillator single particle wave function class.
     * Uses the HarmonicOscillator BasisFunction subclasses auto-generated by SymPy
     * through the orbitalsGenerator tool.
     */
    class AlphaHarmonicOscillator : public Orbitals {

        public:
            AlphaHarmonicOscillator(GeneralParams &, VariationalParams &);

            /*!
             * Calculates the exponential term shared by all oscillator function once pr. particle
             * to save CPU-time.
             * \see Orbitals::set_qnum_indie_terms()
             */
            void set_qnum_indie_terms(mat& r, int i);

            void set_w(double w)
            {
                this->w = w;
                set_parameter(*alpha, 0);
            }

        protected:

            double w; //!< The oscillator frequency.

            double *alpha; //!< Pointer to the variational parameter alpha. Shared address with all the BasisFunction subclasses.
            double *k; //!< Pointer to sqrt(alpha*w). Shared address with all the BasisFunction subclasses.
            double *k2; //!< Pointer to alpha*w. Shared address with all the BasisFunction subclasses.

            double *exp_factor; //!< Pointer to a factor precalculated by set_qnum_indie_terms(). Shared address with all the BasisFunction subclasses.

            /*!
             * Calculates the quantum numbers of the oscillator and stores them in the matrix qnums.
             */
            void get_qnums();
            void setup_basis();

            double get_sp_energy(int qnum) const;

            //! Method for calculating Hermite polynomials.
            /*!
             * @param n The degree of the Hermite polynomial.
             * @param x The argument for evaluating the polynomial.
             */
            double H(int n, double x) const;

            /*!
             * @return The variational parameter alpha.
             */
            double get_parameter(int n) {
                (void) n;

                return *alpha;
            }

            /*!
             * Sets a new value for the alpha and updates all the pointer values.
             */
            void set_parameter(double parameter, int n) {
                (void) n;

                *alpha = parameter;
                *k2 = parameter*w;
                *k = sqrt(*k2);
        }

    };

}
