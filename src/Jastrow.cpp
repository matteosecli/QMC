#include <armadillo>
#include "Jastrow.h"
#include "Structs.h"

using namespace arma;
using namespace QMC2;

Jastrow::Jastrow(int n_p, int dim) {
    this->n_p = n_p;
    this->n2 = ceil(n_p / 2.0);
    this->dim = dim;
}

Jastrow::Jastrow() {}


Pade_Jastrow::Pade_Jastrow(GeneralParams & gP, VariationalParams & vP)
: Jastrow(gP.number_particles, gP.dimension) {

    active = true;

    set_parameter(vP.beta, 0);

    a = arma::zeros<arma::mat > (n_p, n_p);

    initialize();
}

void Pade_Jastrow::initialize() {
    int i, j;
    double a_sym, a_asym;
    a_sym = 1. / 3;
    a_asym = 1.0;

    for (i = 0; i < n_p; i++) {
        for (j = 0; j < n_p; j++) {
            if ((j < n2 && i < n2) || (j >= n2 && i >= n2)) {
                a(i, j) = a_sym;
            } else {
                a(i, j) = a_asym;
            }
        }
    }


}

double Pade_Jastrow::get_val(const mat& r) const {
    int i, j, k;
    double r_ij, r_ij2, arg;

    arg = 0;
    for (i = 0; i < n_p - 1; i++) {
        for (j = i + 1; j < n_p; j++) {
            r_ij2 = 0;
            for (k = 0; k < dim; k++) {
                r_ij2 += pow((r(i,k)-r(j,k)),2);
            }
            r_ij = sqrt(r_ij2);

            arg += a(i, j) * r_ij / (1.0 + beta * r_ij);
        }
    }

    return exp(arg);
}

double Pade_Jastrow::get_grad(const mat r, int particle, int dimension) const {
    double r_ij, r_ij2, arg=0;

    for (int j = 0; j < n_p; j++) {
        if ( j != particle ) {

            r_ij2 = 0;
            for (int k = 0; k < dim; k++) {
                r_ij2 += pow((r(particle,k)-r(j,k)),2);
            }
            r_ij = sqrt(r_ij2);

            arg += ( r(particle, dimension) - r(j,dimension) ) * a(particle, j) / r_ij / pow((1.0 + beta * r_ij),2);

        }
    }

    return arg;
}

double Pade_Jastrow::get_lapl(const mat r, int particle) const {
    double r_ij, r_ij2, arg;
    double grad = 0;

    arg = 0;

    for (int j = 0; j < n_p; j++) {
        if ( j != particle ) {

            r_ij2 = 0;
            for (int k = 0; k < dim; k++) {
                r_ij2 += pow((r(particle,k)-r(j,k)),2);
            }
            r_ij = sqrt(r_ij2);

            arg += a(particle, j) * ( (dim-3)*(1.0 + beta * r_ij) + 2 ) / r_ij / pow((1.0 + beta * r_ij),3);

        }
    }

    for (int d = 0; d < dim; d++) {
        grad = this->get_grad(r, particle, d);
        arg += grad*grad;
    }

    return arg;
}

No_Jastrow::No_Jastrow() {
    n_p = 0;
    active = false;
}
