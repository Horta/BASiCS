#include "nu.h"
#include <armadillo>

using arma::vec;

void Nu::sample() {
  // for (size_t j = 0; j < n; ++j) {
  //   posterior_loglikelihood(data.X(j), param.delta, param.phi(j), param.y(j),
  //                           param.mu, param.nu0(j));
  // }
}

void posterior_loglikelihood(const vec& x, const vec& delta, double phi,
                             double y, const vec& mu, double nu0) {
  auto deltai = 1 / delta;
  auto thetai = 1 / theta;
  auto a = phi * y * mu + deltai;
  auto b = phi * nu0 * mu + deltai;
  double r0 = -sum((x + deltai) * log(a / b));
  double r1 = (log(y) - log(nu0)) * (sum(x) + thetai);
  double r2 = - (y - nu0) * (SumSpikeInput + (thetai/sj);
  return r0 + r1 + r2;
}
