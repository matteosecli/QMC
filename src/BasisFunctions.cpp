#include <armadillo>
#include "BasisFunctions.h"

using namespace QMC2;

// Superclass constructor
HarmonicOscillator::HarmonicOscillator(double* k, double* k2, double* exp_factor) {
    this->k = k;
    this->k2 = k2;
    this->exp_factor = exp_factor;
}


// Subclass constructors
HarmonicOscillator_0::HarmonicOscillator_0(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
HarmonicOscillator_1::HarmonicOscillator_1(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
HarmonicOscillator_2::HarmonicOscillator_2(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
HarmonicOscillator_3::HarmonicOscillator_3(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
HarmonicOscillator_4::HarmonicOscillator_4(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
HarmonicOscillator_5::HarmonicOscillator_5(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}


// Evals    // H_nx,ny // CHECK THIS STUFF!!!
double HarmonicOscillator_0::eval(const mat& r, int i) {    // H_0,0
    (void) r;
    (void) i;
    //exp(-k^2*r^2/2)
    H = 1;
    return H*(*exp_factor);
}
double HarmonicOscillator_1::eval(const mat& r, int i) {    // H_0,1
    y = r(i, 1);
    //y*exp(-k^2*r^2/2)
    H = y;
    return H*(*exp_factor);
}
double HarmonicOscillator_2::eval(const mat& r, int i) {    // H_1,0
    x = r(i, 0);
    //x*exp(-k^2*r^2/2)
    H = x;
    return H*(*exp_factor);
}
double HarmonicOscillator_3::eval(const mat& r, int i) {    // H_0,2
    y = r(i, 1);
    y2 = y*y;
    //(2*k^2*y^2 - 1)*exp(-k^2*r^2/2)
    H = 2*(*k2)*y2 - 1;
    return H*(*exp_factor);
}
double HarmonicOscillator_4::eval(const mat& r, int i) {    // H_1,1
    x = r(i, 0);
    y = r(i, 1);
    //x*y*exp(-k^2*r^2/2)
    H = x*y;
    return H*(*exp_factor);
}
double HarmonicOscillator_5::eval(const mat& r, int i) {    // H_2,0
    x = r(i, 0);
    x2 = x*x;
    //(2*k^2*x^2 - 1)*exp(-k^2*r^2/2)
    H = 2*(*k2)*x2 - 1;
    return H*(*exp_factor);
}
