#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

/* WARNING! THE CLASS "HARMONICOSCILLATOR" ONLY WORKS
 * IN THE 2-DIM CASE! DO NOT USE IT FOR THE 3-DIM
 * CASE!
 */

using namespace arma;

namespace QMC2
{

    // Basis functions superclass
    class BasisFunctions {
    public:
        BasisFunctions() {}
        virtual double eval(const mat& r, int i) = 0;
    };


    // Harmonic oscillator specfic basis functions
    class HarmonicOscillator : public BasisFunctions {
    protected:
        double* k;
        double* k2;
        double* exp_factor;
        double H;
        double x;
        double y;
        double x2;
        double y2;
    public:
        HarmonicOscillator(double* k, double* k2, double* exp_factor);
    };


    /* Now we need an HarmonicOscillator function (packed into a class
     * for convenience) for each one of the f***ing particles.
     * IMPORTANT: Each BasisFunction contains 2 electrons! So, we only
     * need 6 BasisFunctions to deal with 12 electrons.
     */
    class HarmonicOscillator_0 : public HarmonicOscillator {
    public:
        HarmonicOscillator_0(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat &r, int i);
    };

    class HarmonicOscillator_1 : public HarmonicOscillator {
    public:
        HarmonicOscillator_1(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat& r, int i);
    };

    class HarmonicOscillator_2 : public HarmonicOscillator {
    public:
        HarmonicOscillator_2(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat& r, int i);
    };

    class HarmonicOscillator_3 : public HarmonicOscillator {
    public:
        HarmonicOscillator_3(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat& r, int i);
    };

    class HarmonicOscillator_4 : public HarmonicOscillator {
    public:
        HarmonicOscillator_4(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat& r, int i);
    };

    class HarmonicOscillator_5 : public HarmonicOscillator {
    public:
        HarmonicOscillator_5(double* k, double* k2, double* exp_factor);
        virtual double eval(const mat& r, int i);
    };

}

#endif // BASISFUNCTIONS_H
