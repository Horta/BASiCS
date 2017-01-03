#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <armadillo>

#include "gig/gig.h"

#include "phi.h"
#include "s.h"

using arma::vec;
using arma::mat;

double kappaj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                          const Phij &phij, double nuj, const Sj &sj,
                          double theta, double kappa_var);

double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       const Phij &phij, double nuj, const Sj &sj,
                       double theta);

double sj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                      const Phij &phij, double nuj, const Sj &sj, double theta,
                      Random &random);

#endif
