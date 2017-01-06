#include "posterior.h"

#include <cmath>

using arma::log;
using arma::accu;
using arma::find;
using arma::uvec;

using std::lgamma;

double mui_post_loglik(const vec &xi, double mui, const Deltai &deltai, const Phi &phi,
                       const vec &nu) {
  double left = accu(xi - 1) * mui;
  double right =
      accu((xi + 1 / deltai.get()) % log(phi.get() % nu * mui + 1 / deltai.get()));

  return left - right;
}

double deltai_post_loglik(const vec &xi, double mui, const Deltai &deltai,
                          const Phi &phi, const vec &nu) {
  double a = xi.n_elem * std::lgamma(1 / deltai.get());

  vec b0 = arma::lgamma(xi + 1 / deltai.get());
  vec b1 = (xi + 1 / deltai.get()) % log(phi.get() % nu * mui + 1 / deltai.get());

  return -a + accu(b0 - b1);
}

double theta_post_loglik(size_t n, double nuj, const Sj &sj, double theta) {
  double a = accu(log(nuj) - log(sj.get())) / theta;
  double b = n * std::lgamma(1 / theta);
  return a - b;
}

double kappaj_post_loglik(const vec &xj, const vec &mu, const Delta &delta,
                          const Phij &phij, double nuj, const Sj &sj,
                          double theta, double kappa_var) {
  uvec pos = find(delta.get() > 0);

  vec id = 1 / delta.get()(pos);

  double phij_ = phij.get();
  double a = accu(xj(pos)) * log(phij_);
  double b = -accu((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id));

  double kappaj = phij.get_kappa();

  return a + b - kappaj * kappaj / (2 * kappa_var);
}

double nuj_post_loglik(const vec &xj, const vec &mu, const Delta &delta,
                       const Phij &phij, double nuj, const Sj &sj,
                       double theta) {
  uvec pos = find(delta.get() > 0);
  uvec zer = find(delta.get() == 0);

  size_t q = xj.size();
  vec id = 1 / delta.get()(pos);

  double a = (accu(xj) + q / theta - q) * log(nuj);

  double b0 = -accu((xj(pos) + id) % log(mu(pos) * phij.get() * nuj + id));
  double b1 = -nuj * accu(mu(zer) + 1 / (theta * sj.get()));
  return a + b0 + b1;
}

double sj_post_loglik(const vec &xj, const vec &mu, const Delta &delta,
                      const Phij &phij, double nuj, const Sj &sj, double theta,
                      Random &random) {
  double p = sj.gamma_shape() - 1 / theta;
  double b = 2 * sj.gamma_rate();
  return random.gig(p, 2 * nuj / theta, b);
}
