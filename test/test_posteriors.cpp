#include "basics/basics.h"
#include "basics/phi.h"
#include "basics/posterior.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

using std::abs;
using arma::span;

bool isclose(double x, double y) {
  double rtol = 1e-7;
  double atol = 0;
  return abs(x - y) <= atol + rtol * abs(y);
}

void test_nu() {
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

void test_kappa() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});
  vec s({2.1, 2.2, 1.9});
  double theta = 0.64;

  vec kappa({1.3, 1.36});
  double kappa_var = 1;

  Phi phi(3);
  phi.set(kappa);

  assert(isclose(kappaj_post_loglik(X(span::all, 0), mu, delta, phi.get()(0),
                                    nu(0), s(0), theta, kappa(0), kappa_var),
                 -11.031259249539308698));
}

int main() {
  test_nu();
  test_kappa();

  return 0;
}
