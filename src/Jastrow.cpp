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

No_Jastrow::No_Jastrow() {
    n_p = 0;
    active = false;
}
