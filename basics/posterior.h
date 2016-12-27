#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <armadillo>

#include "phi.h"

using arma::vec;
using arma::mat;

double kappaj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                          const Phij& phij, double nuj, double sj, double theta,
                          double kappa_var);

double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       const Phij& phij, double nuj, double sj, double theta);

#endif
