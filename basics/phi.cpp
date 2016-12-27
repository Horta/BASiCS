#include "phi.h"

using arma::exp;
using arma::sum;

vec Phi::get(void)
{
  vec phi(kappa.size() + 1);
  phi(0) = phi.size();
  phi.tail(kappa.size()) = phi(0) * exp(kappa) / sum(exp(kappa));
  return phi;
}
