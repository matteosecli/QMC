/* 
 * File:   AlphaHarmonicOscillator.cpp
 * Author: jorgehog
 * 
 * Created on 26. juni 2012, 17:41
 */

/*
 * This is a simplified version of
 * Jorgen's original file, suited
 * to this project
 */

#include "AlphaHarmonicOscillator.h"

#include "Structs.h"
#include "BasisFunctions.h"

using namespace QMC2;


Orbitals::Orbitals(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = ceil(n_p / 2.0);
    this->n_half = n_p/2;
    this->dim = dim;

    max_implemented = 6; //for 12 particles
    basis_functions = new BasisFunctions*[max_implemented];
    dell_basis_functions = new BasisFunctions**[dim];
    for (int i = 0; i < dim; i++) {
        dell_basis_functions[i] = new BasisFunctions*[max_implemented];
    }
    lapl_basis_functions = new BasisFunctions*[max_implemented];

    nCap = n_p;
}


Orbitals::Orbitals() {}


double Orbitals::phi(const mat& r, int particle, int q_num) {
    return basis_functions[q_num]->eval(r, particle);
}

double Orbitals::del_phi(const mat& r, int particle, int q_num, int d) {
    return dell_basis_functions[d][q_num]->eval(r, particle);
}

double Orbitals::lapl_phi(const mat& r, int particle, int q_num) {
    return lapl_basis_functions[q_num]->eval(r, particle);
}

double Orbitals::SlaterD(mat& r) {
    mat D_up(n_half,n_half), D_down(n_half,n_half);

    for (int i = 0; i < n_half; i++) {
        this->set_qnum_indie_terms(r, i);
        for (int j = 0; j < n_half; j++) {
            D_up(i,j) = this->phi(r, i, j);
        }
    }

    for (int i = 0; i < n_half; i++) {
        this->set_qnum_indie_terms(r, i+n_half);
        for (int j = 0; j < n_half; j++) {
            D_down(i,j) = this->phi(r, i+n_half, j);
        }
    }

    double slater = det(D_up)*det(D_down);

    D_up.reset();
    D_down.reset();

    return slater;
}

double Orbitals::SlaterD_grad(mat r, int particle, int dimension){
    int offset;
    double grad = 0;
    mat S(n_half,n_half);
    mat S_inv(n_half,n_half);

    if ( particle < n_half ) {
        offset = 0;
    } else {
        offset = n_half;
    }

    for (int k = 0; k < n_half; k++) {
        this->set_qnum_indie_terms(r, k+offset);
        for (int j = 0; j < n_half; j++) {
            S(k,j) = this->phi(r, k+offset, j);
        }
    }

    S_inv = inv(S);

    this->set_qnum_indie_terms(r, particle);
    for (int k = 0; k < n_half; k++) {
        grad += this->del_phi(r, particle, k, dimension)*S_inv(k,particle);
    }

    S.reset();
    S_inv.reset();

    return grad;
}

double Orbitals::SlaterD_lapl(mat r, int particle) {
    int offset;
    mat S(n_half,n_half);
    mat S_inv(n_half,n_half);
    double lapl = 0;

    if ( particle < n_half ) {
        offset = 0;
    } else {
        offset = n_half;
    }

    for (int k = 0; k < n_half; k++) {
        this->set_qnum_indie_terms(r, k+offset);
        for (int j = 0; j < n_half; j++) {
            S(k,j) = this->phi(r, k+offset, j);
        }
    }

    S_inv = inv(S);

    this->set_qnum_indie_terms(r, particle);
    for (int k = 0; k < n_half; k++) {
        lapl += this->lapl_phi(r, particle, k)*S_inv(k,particle);
    }

    S.reset();
    S_inv.reset();

    return lapl;
}


AlphaHarmonicOscillator::AlphaHarmonicOscillator(GeneralParams & gP, VariationalParams & vP)
: Orbitals(gP.number_particles, gP.dimension) {

    this->alpha = new double();
    this->k = new double();
    this->k2 = new double();
    this->exp_factor = new double();

    this->w = gP.frequency;
    set_parameter(vP.alpha, 0);

    qnums = arma::zeros<arma::imat > (n2, dim);

    name = "QDots";
    get_qnums();
    setup_basis();

}


void AlphaHarmonicOscillator::setup_basis() {

    basis_functions[0] = new HarmonicOscillator_0(k, k2, exp_factor);
    basis_functions[1] = new HarmonicOscillator_1(k, k2, exp_factor);
    basis_functions[2] = new HarmonicOscillator_2(k, k2, exp_factor);
    basis_functions[3] = new HarmonicOscillator_3(k, k2, exp_factor);
    basis_functions[4] = new HarmonicOscillator_4(k, k2, exp_factor);
    basis_functions[5] = new HarmonicOscillator_5(k, k2, exp_factor);

    dell_basis_functions[0][0] = new dell_HarmonicOscillator_0_x(k, k2, exp_factor);
    dell_basis_functions[1][0] = new dell_HarmonicOscillator_0_y(k, k2, exp_factor);
    dell_basis_functions[0][1] = new dell_HarmonicOscillator_1_x(k, k2, exp_factor);
    dell_basis_functions[1][1] = new dell_HarmonicOscillator_1_y(k, k2, exp_factor);
    dell_basis_functions[0][2] = new dell_HarmonicOscillator_2_x(k, k2, exp_factor);
    dell_basis_functions[1][2] = new dell_HarmonicOscillator_2_y(k, k2, exp_factor);
    dell_basis_functions[0][3] = new dell_HarmonicOscillator_3_x(k, k2, exp_factor);
    dell_basis_functions[1][3] = new dell_HarmonicOscillator_3_y(k, k2, exp_factor);
    dell_basis_functions[0][4] = new dell_HarmonicOscillator_4_x(k, k2, exp_factor);
    dell_basis_functions[1][4] = new dell_HarmonicOscillator_4_y(k, k2, exp_factor);
    dell_basis_functions[0][5] = new dell_HarmonicOscillator_5_x(k, k2, exp_factor);
    dell_basis_functions[1][5] = new dell_HarmonicOscillator_5_y(k, k2, exp_factor);

    lapl_basis_functions[0] = new lapl_HarmonicOscillator_0(k, k2, exp_factor);
    lapl_basis_functions[1] = new lapl_HarmonicOscillator_1(k, k2, exp_factor);
    lapl_basis_functions[2] = new lapl_HarmonicOscillator_2(k, k2, exp_factor);
    lapl_basis_functions[3] = new lapl_HarmonicOscillator_3(k, k2, exp_factor);
    lapl_basis_functions[4] = new lapl_HarmonicOscillator_4(k, k2, exp_factor);
    lapl_basis_functions[5] = new lapl_HarmonicOscillator_5(k, k2, exp_factor);

}


void AlphaHarmonicOscillator::get_qnums() {
    double n_x, n_y;

    int n_shells = (int) (0.5 * (sqrt(1 + 4 * n_p) - 1));

    int q = 0;
    for (int shell = 0; shell < n_shells; shell++) {

        n_x = 0;
        n_y = shell;

        for (int subshell_i = 0; subshell_i <= shell; subshell_i++) {
            qnums(q, 0) = n_x;
            qnums(q, 1) = n_y;

            n_x++;
            n_y--;

            q++;
        }
    }
}


double AlphaHarmonicOscillator::H(int n, double x) const {
    if (n < 0) {
        return 0;
    } else if (n == 0) {
        return 1;
    } else if (n == 1) {
        return 2 * x;
    } else if (n == 2) {
        return 4 * x * x - 2;
    } else {
        std::cout << "Unsopported Hermite polynomial level: " << n << std::endl;
        exit(1);
    }
}


void AlphaHarmonicOscillator::set_qnum_indie_terms(mat &r, int i) {
    double r_i2 = pow(r(i,0),2) + pow(r(i,1),2);
    *exp_factor = exp(-0.5 * (*k2) * r_i2);
}


double AlphaHarmonicOscillator::get_sp_energy(int qnum) const {
    int n_x = qnums(qnum, 0);
    int n_y = qnums(qnum, 1);

    return w * (n_x + n_y + 1);

}
