#include "basics/basics.h"
#include "basics/s.h"
#include "basics/phi.h"
#include "basics/posterior.h"

#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>


using std::cout;
using std::endl;
using std::flush;
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
  double nuj = 0.9;
  double sj = 2.1;
  double theta = 0.64;

  Phi phi(3);
  phi.set({1.6008, 1.05425});

  assert(isclose(nuj_post_loglik(xj, mu, delta, phi.phij(2), nuj, sj, theta),
                 -18.2090911836940598789169598604));
}

void test_kappa() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  cout << "X" << X << endl << flush;
  vec mu({-0.2, 1.1});
  cout << "mu" << mu << endl << flush;
  vec delta({0.2, 2.1});
  cout << "delta" << delta << endl << flush;

  vec nu({0.9, 0.7, 1.1});
  cout << "nu" << nu << endl << flush;
  vec s({2.1, 2.2, 1.9});
  cout << "s" << s << endl << flush;
  double theta = 0.64;

  vec kappa({1.3, 1.36});
  cout << "kappa" << kappa << endl << flush;
  double kappa_var = 1;

  Phi phi(3);
  phi.set(kappa);

  cout << "passou" << endl << flush;
  assert(isclose(kappaj_post_loglik(X(span::all, 1), mu, delta, phi.phij(1),
                                    nu(1), s(1), theta, kappa_var),
                 -12.671154074905397025));
}

void test_s()
{
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});
  S s(3, 1.1, 2.1);
  s.get() = {2.1, 2.2, 1.9};

  std::cout << s.get() << std::endl;

  double theta = 0.64;

  vec kappa({1.3, 1.36});
  double kappa_var = 1;

  Phi phi(3);
  phi.set(kappa);
}

int main() {
  cout << "Testing nu." << endl << flush;
  test_nu();

  cout << "Testing kappa." << endl << flush;
  test_kappa();

  cout << "Testing s." << endl << flush;
  test_s();

  return 0;
}
