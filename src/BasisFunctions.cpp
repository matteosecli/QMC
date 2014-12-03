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
dell_HarmonicOscillator_0_x::dell_HarmonicOscillator_0_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_0_y::dell_HarmonicOscillator_0_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_0::lapl_HarmonicOscillator_0(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 0 ---//
HarmonicOscillator_1::HarmonicOscillator_1(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_1_x::dell_HarmonicOscillator_1_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_1_y::dell_HarmonicOscillator_1_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_1::lapl_HarmonicOscillator_1(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 1 ---//
HarmonicOscillator_2::HarmonicOscillator_2(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_2_x::dell_HarmonicOscillator_2_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_2_y::dell_HarmonicOscillator_2_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_2::lapl_HarmonicOscillator_2(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 2 ---//
HarmonicOscillator_3::HarmonicOscillator_3(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_3_x::dell_HarmonicOscillator_3_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_3_y::dell_HarmonicOscillator_3_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_3::lapl_HarmonicOscillator_3(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 3 ---//
HarmonicOscillator_4::HarmonicOscillator_4(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_4_x::dell_HarmonicOscillator_4_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_4_y::dell_HarmonicOscillator_4_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_4::lapl_HarmonicOscillator_4(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 4 ---//
HarmonicOscillator_5::HarmonicOscillator_5(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_5_x::dell_HarmonicOscillator_5_x(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
dell_HarmonicOscillator_5_y::dell_HarmonicOscillator_5_y(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
lapl_HarmonicOscillator_5::lapl_HarmonicOscillator_5(double* k, double* k2, double* exp_factor)
    : HarmonicOscillator(k, k2, exp_factor) {}
//--- END 5 ---//



// Evals    // H_nx,ny // CHECK THIS STUFF!!!
double HarmonicOscillator_0::eval(const mat& r, int i) {    // H_0,0
    (void) r;
    (void) i;
    //exp(-k^2*r^2/2)
    H = 1;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_0_x::eval(const mat& r, int i) {
    x = r(i, 0);
    //-k^2*x*exp(-k^2*r^2/2)
    H = -(*k2)*x;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_0_y::eval(const mat& r, int i) {
    y = r(i, 1);
    //-k^2*y*exp(-k^2*r^2/2)
    H = -(*k2)*y;
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_0::eval(const mat& r, int i) {
    //k^2*(k^2*r^2 - 2)*exp(-k^2*r^2/2)
    H = (*k2)*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 2);
    return H*(*exp_factor);
}
//--- END 0 ---//

double HarmonicOscillator_1::eval(const mat& r, int i) {    // H_0,1
    y = r(i, 1);
    //y*exp(-k^2*r^2/2)
    H = y;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_1_x::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    //-k^2*x*y*exp(-k^2*r^2/2)
    H = -(*k2)*x*y;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_1_y::eval(const mat& r, int i) {
    y = r(i, 1);
    //-(k*y - 1)*(k*y + 1)*exp(-k^2*r^2/2)
    H = -((*k)*y - 1)*((*k)*y + 1);
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_1::eval(const mat& r, int i) {
    y = r(i, 1);
    //k^2*y*(k^2*r^2 - 4)*exp(-k^2*r^2/2)
    H = (*k2)*y*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 4);
    return H*(*exp_factor);
}
//--- END 1 ---//

double HarmonicOscillator_2::eval(const mat& r, int i) {    // H_1,0
    x = r(i, 0);
    //x*exp(-k^2*r^2/2)
    H = x;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_2_x::eval(const mat& r, int i) {
    x = r(i, 0);
    //-(k*x - 1)*(k*x + 1)*exp(-k^2*r^2/2)
    H = -((*k)*x - 1)*((*k)*x + 1);
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_2_y::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    //-k^2*x*y*exp(-k^2*r^2/2)
    H = -(*k2)*x*y;
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_2::eval(const mat& r, int i) {
    x = r(i, 0);
    //k^2*x*(k^2*r^2 - 4)*exp(-k^2*r^2/2)
    H = (*k2)*x*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 4);
    return H*(*exp_factor);
}
//--- END 2 ---//

double HarmonicOscillator_3::eval(const mat& r, int i) {    // H_0,2
    y = r(i, 1);
    y2 = y*y;
    //(2*k^2*y^2 - 1)*exp(-k^2*r^2/2)
    H = 2*(*k2)*y2 - 1;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_3_x::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    y2 = y*y;
    //-k^2*x*(2*k^2*y^2 - 1)*exp(-k^2*r^2/2)
    H = -(*k2)*x*(2*(*k2)*y2 - 1);
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_3_y::eval(const mat& r, int i) {
    y = r(i, 1);
    y2 = y*y;
    //-k^2*y*(2*k^2*y^2 - 5)*exp(-k^2*r^2/2)
    H = -(*k2)*y*(2*(*k2)*y2 - 5);
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_3::eval(const mat& r, int i) {
    y = r(i, 1);
    y2 = y*y;
    //k^2*(k^2*r^2 - 6)*(2*k^2*y^2 - 1)*exp(-k^2*r^2/2)
    H = (*k2)*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 6)*(2*(*k2)*y2 - 1);
    return H*(*exp_factor);
}
//--- END 3 ---//

double HarmonicOscillator_4::eval(const mat& r, int i) {    // H_1,1
    x = r(i, 0);
    y = r(i, 1);
    //x*y*exp(-k^2*r^2/2)
    H = x*y;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_4_x::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    //-y*(k*x - 1)*(k*x + 1)*exp(-k^2*r^2/2)
    H = -y*((*k)*x - 1)*((*k)*x + 1);
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_4_y::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    //-x*(k*y - 1)*(k*y + 1)*exp(-k^2*r^2/2)
    H = -x*((*k)*y - 1)*((*k)*y + 1);
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_4::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    //k^2*x*y*(k^2*r^2 - 6)*exp(-k^2*r^2/2)
    H = (*k2)*x*y*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 6);
    return H*(*exp_factor);
}
//--- END 4 ---//

double HarmonicOscillator_5::eval(const mat& r, int i) {    // H_2,0
    x = r(i, 0);
    x2 = x*x;
    //(2*k^2*x^2 - 1)*exp(-k^2*r^2/2)
    H = 2*(*k2)*x2 - 1;
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_5_x::eval(const mat& r, int i) {
    x = r(i, 0);
    x2 = x*x;
    //-k^2*x*(2*k^2*x^2 - 5)*exp(-k^2*r^2/2)
    H = -(*k2)*x*(2*(*k2)*x2 - 5);
    return H*(*exp_factor);
}
double dell_HarmonicOscillator_5_y::eval(const mat& r, int i) {
    x = r(i, 0);
    y = r(i, 1);
    x2 = x*x;
    //-k^2*y*(2*k^2*x^2 - 1)*exp(-k^2*r^2/2)
    H = -(*k2)*y*(2*(*k2)*x2 - 1);
    return H*(*exp_factor);
}
double lapl_HarmonicOscillator_5::eval(const mat& r, int i) {
    x = r(i, 0);
    x2 = x*x;
    //k^2*(k^2*r^2 - 6)*(2*k^2*x^2 - 1)*exp(-k^2*r^2/2)
    H = (*k2)*((*k2)*(r(i,0)*r(i,0)+r(i,1)*r(i,1)) - 6)*(2*(*k2)*x2 - 1);
    return H*(*exp_factor);
}
//--- END 5 ---//
