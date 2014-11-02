
using namespace arma;

namespace QMC2
{

    struct GeneralParams;
    struct VariationalParams;

    /*! \brief The class representing the Jastrow correlation functions
     * Holds all data concerning the Jastrow function and it's influence
     * on the QMC algorithm.
     */
    class Jastrow {
    protected:
        int n_p;
        int n2;
        int dim;

    public:
        Jastrow(int n_p, int dim);
        Jastrow();

        bool active; //!< Parameter false if No_Jastrow is loaded.

        //! Sets variational parameters
        /*!
         * @param n The index of the sought variational parameter
         * @param param The new value of parameter [n]
         */
        virtual void set_parameter(double param, int n) = 0;

        //! Returns variational parameters
        /*!
         * @param n The index of the sought variational parameter
         * @return Variational parameter with index [n]
         */
        virtual double get_parameter(int n) = 0;

        /*!
         * Initializes the non-variational parameters needed by the Jastrow Factor.
         */
        virtual void initialize() = 0;

        /*!
         * Calculates the value of the Jastrow Factor at the walker's position.
         */
        virtual double get_val(const mat& r) const = 0;

    };


    /*!
     * \brief The Pade Jastrow factor with a single variational parameter.
     */
    class Pade_Jastrow : public Jastrow {
    protected:
        double beta; //!< The variational parameter
        arma::mat a; //!< The spin-dependent variables taking care of the cusp condition.

        void set_parameter(double param, int n) {
            (void) n;
            beta = param;
        }

        double get_parameter(int n) {
            (void) n;

            return beta;
        }

    public:

        Pade_Jastrow(GeneralParams &, VariationalParams &);

        /*!
         * In case of Pade Jastrow, initializing means setting up the a matrix.
         */
        void initialize();

        double get_val(const mat& r) const;

    };


    /*!
     * \brief Class loaded when no correlation factor is used.
     */
    class No_Jastrow : public Jastrow {
    protected:

        double get_parameter(int n) {
            (void) n;

            return 0;
        }

        void set_parameter(double param, int n) {
            (void) param;
            (void) n;
        }

    public:

        No_Jastrow();

        void initialize() {}

        double get_val(const mat& r) const {
            (void) r;

            return 1;
        }

    };


}
