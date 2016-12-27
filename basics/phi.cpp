#include "phi.h"

using arma::exp;
using arma::sum;

vec Phi::get(void) const
{
  vec phi(kappa.size() + 1);
  phi(0) = phi.size();
  phi.tail(kappa.size()) = phi(0) * exp(kappa) / sum(exp(kappa));
  return phi;
}


Phij Phi::phij(size_t j) const
{
  double phi0 = kappa.size() + 1;
  return Phij(kappa(j-1), phi0, sum(exp(kappa)));
}
