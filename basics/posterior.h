#ifndef POSTERIOR_H
#define POSTERIOR_H

#include <armadillo>

#include "gig/gig.h"

#include "phi.h"
#include "s.h"

using arma::mat;

double mui_post_loglik(const mat &xi, double mui, double deltai, const Phi &phi,
                       const mat &nu);

double theta_post_loglik(size_t n, double nuj, const Sj &sj, double theta);

double deltai_post_loglik(const mat &xi, double mui, double deltai,
                          const Phi &phi, const mat &nu);

double kappaj_post_loglik(const mat &xj, const mat &mu, const mat &delta,
                          const Phij &phij, double nuj, const Sj &sj,
                          double theta, double kappa_var);

double nuj_post_loglik(const mat &xj, const mat &mu, const mat &delta,
                       const Phij &phij, double nuj, const Sj &sj,
                       double theta);

double sj_post_loglik(const mat &xj, const mat &mu, const mat &delta,
                      const Phij &phij, double nuj, const Sj &sj, double theta,
                      Random &random);

#endif
