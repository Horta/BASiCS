#include "basics/basics.h"
#include "basics/phi.h"
#include "basics/posterior.h"
#include "basics/s.h"

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

void test_mu() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});

  vec kappa({1.3, 1.36});

  Phi phi(3);
  phi.set(kappa);

  std::mt19937_64 generator(0);
  Random random(generator);

  double mui = mui_post_loglik(X(1, span::all).t(), mu(1), delta(1), phi, nu);
  assert(isclose(mui, -2.8960000990179226577));
}

void test_delta() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});

  vec kappa({1.3, 1.36});

  Phi phi(3);
  phi.set(kappa);

  std::mt19937_64 generator(0);
  Random random(generator);

  double r = deltai_post_loglik(X(1, span::all).t(), mu(1), delta(1), phi, nu);
  assert(isclose(r, -3.3039472650673218368));
}

void test_theta() {
  vec nu({0.9, 0.7, 1.1});

  S s(3, 1.0, 1.0);
  s.get() = {1.1, 1.2, 2.1};
  double theta = 0.64;

  double r = theta_post_loglik(3, nu(1), s.sj(1), theta);

  assert(isclose(r, -0.49205598225087332498));
}

void test_nu() {
  vec xj({5.0, 3.0});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});
  double nuj = 0.9;
  S s(3, 1.0, 1.0);
  s.get() = {1.1, 1.2, 2.1};
  double theta = 0.64;

  Phi phi(3);
  phi.set({1.6008, 1.05425});

  double r = nuj_post_loglik(xj, mu, delta, phi.phij(2), nuj, s.sj(2), theta);
  assert(isclose(r, -18.2090911836940598789169598604));

  delta(1) = 0;
  r = nuj_post_loglik(xj, mu, delta, phi.phij(2), nuj, s.sj(2), theta);
  assert(isclose(r, -18.31138418403616086));

  s.get()(2) = 3.1;
  r = nuj_post_loglik(xj, mu, delta, phi.phij(2), nuj, s.sj(2), theta);
  assert(isclose(r, -18.09537035915137082));
}

void test_kappa() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});

  S s(3, 1.0, 1.0);
  s.get() = {2.1, 2.2, 1.9};
  double theta = 0.64;

  vec kappa({1.3, 1.36});
  double kappa_var = 1;

  Phi phi(3);
  phi.set(kappa);

  double r = kappaj_post_loglik(X(span::all, 1), mu, delta, phi.phij(1), nu(1),
                                s.sj(1), theta, kappa_var);
  assert(isclose(r, -12.671154074905397025));
}

void test_s() {
  mat X({{5.0, 3.0, 2.0}, {1.0, 2.0, 2.0}});
  vec mu({-0.2, 1.1});
  vec delta({0.2, 2.1});

  vec nu({0.9, 0.7, 1.1});
  S s(3, 1.1, 2.1);
  s.get() = {2.1, 2.2, 1.9};

  double theta = 0.64;

  vec kappa({1.3, 1.36});
  double kappa_var = 1;

  Phi phi(3);
  phi.set(kappa);

  std::mt19937_64 generator(0);
  Random random(generator);

  double r = sj_post_loglik(X(span::all, 1), mu, delta, phi.phij(1), nu(1),
                            s.sj(1), theta, random);
  assert(isclose(r, 1.6883913147833378154));
}

int main() {
  cout << "Testing mu." << endl << flush;
  test_mu();

  cout << "Testing nu." << endl << flush;
  test_nu();

  cout << "Testing kappa." << endl << flush;
  test_kappa();

  cout << "Testing s." << endl << flush;
  test_s();

  cout << "Testing delta." << endl << flush;
  test_delta();

  cout << "Testing theta." << endl << flush;
  test_theta();

  return 0;
}
