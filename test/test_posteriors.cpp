#include "basics/basics.h"
#include "basics/posterior.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

using std::abs;

bool isclose(double x, double y) {
  double rtol = 1e-7;
  double atol = 0;
  return abs(x - y) <= atol + rtol * abs(y);
}

void test_nu()
{
  vec xj({5.0, 3.0});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});
  double phij = 1.1;
  double nuj = 0.9;
  double sj = 2.1;
  double theta = 0.64;

  assert(isclose(nuj_post_loglik(xj, mu, delta, phij, nuj, sj, theta),
                 -18.209099180638268933));
}

int main() {
  test_nu();

  return 0;
}
