#include "posterior.h"

using arma::log;
using arma::sum;
using arma::find;
using arma::uvec;

double mui_post_loglik(const vec &xi, double mui, double deltai, const Phi &phi,
                       const vec &nu) {
  double left = sum(xi - 1) * mui;
  double right =
      sum((xi + 1 / deltai) % log(phi.get() % nu * mui + 1 / deltai));

  return left - right;
}

double kappaj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                          const Phij &phij, double nuj, const Sj &sj,
                          double theta, double kappa_var) {
  uvec pos = find(delta > 0);

  vec id = 1 / delta(pos);

  double phij_ = phij.get();
  double a = sum(xj(pos)) * log(phij_);
  double b = -sum((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id));

  double kappaj = phij.get_kappa();

  return a + b - kappaj * kappaj / (2 * kappa_var);
}

double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       const Phij &phij, double nuj, const Sj &sj,
                       double theta) {
  uvec pos = find(delta > 0);
  uvec zer = find(delta == 0);

  size_t q = xj.size();
  vec id = 1 / delta(pos);

  double a = (sum(xj) + q / theta - q) * log(nuj);

  double b0 = -sum((xj(pos) + id) % log(mu(pos) * phij.get() * nuj + id));
  double b1 = -nuj * sum(mu(zer) + 1 / (theta * sj.get()));
  return a + b0 + b1;
}

double sj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                      const Phij &phij, double nuj, const Sj &sj, double theta,
                      Random &random) {
  double p = sj.gamma_shape() - 1 / theta;
  double b = 2 * sj.gamma_rate();
  return random.gig(p, 2 * nuj / theta, b);
}
