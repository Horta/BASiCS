#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <armadillo>

using arma::vec;

double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       double phij, double nuj, double sj, double theta);

#endif
