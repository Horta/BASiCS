#include "nu.h"

void Nu::sample() {
  // for (size_t j = 0; j < n; ++j) {
  //   posterior_loglikelihood(data.X(j), param.delta, param.phi(j), param.y(j),
  //                           param.mu, param.nu0(j));
  // }
}

// void posterior_loglikelihood(arma::vec Xj, arma::vec delta, double phij,
//                              double yj, arma::vec mu, double nu0j) {
//   auto deltai = 1 / delta;
//   auto thetai = 1 / theta;
//   auto a = phij * yj * mu + deltai;
//   auto b = phij * nu0j * mu + deltai;
//   double r0 = -sum((Xj + deltai) * log(a / b));
//   double r1 = (log(yj) - log(nu0j)) * (sum(Xj) + thetai);
//   double r2 = - (yj - nu0j) * (SumSpikeInput + (thetai/sj);
//   return r0 + r1 + r2;
// }
