#include "posterior.h"

using arma::log;
using arma::sum;
using arma::find;

double kappaj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                          const Phij& phij, double nuj, double sj, double theta,
                          double kappa_var) {
  auto pos = find(delta > 0);

  auto id = 1 / delta(pos);

  double phij_ = phij.get();

  double a = sum(xj(pos)) * log(phij_);
  double b = -sum((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id));

  double kappaj = phij.get_kappa();
  
  return a + b - kappaj * kappaj / (2 * kappa_var);
}


double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       const Phij& phij, double nuj, double sj, double theta) {
  auto pos = find(delta > 0);
  auto zer = find(delta == 0);

  auto q = xj.size();
  auto id = 1 / delta(pos);

  double a = (sum(xj) + q / theta - q) * log(nuj);

  double b0 = -sum((xj(pos) + id) % log(mu(pos) * phij.get() * nuj + id));
  double b1 = -nuj * sum(mu(zer) + 1 / (theta * sj));

  return a + b0 + b1;
}
