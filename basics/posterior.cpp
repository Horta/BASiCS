#include "posterior.h"

using arma::log;
using arma::sum;
using arma::find;
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

double kappaj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                          const Phij &phij, double nuj, double sj, double theta,
                          double kappa_var) {
  cout << "ponto 1" << endl << flush;
  auto pos = find(delta > 0);

  cout << "ponto 2" << endl << flush;
  auto id = 1 / delta(pos);

  cout << "ponto 3" << endl << flush;
  double phij_ = phij.get();

  cout << "ponto 4" << endl << flush;

  cout << "pos" << pos << endl << flush;
  cout << "xj(pos)" << xj(pos) << endl << flush;
  cout << "log(phij_)" << log(phij_) << endl << flush;

  double a = sum(xj(pos)) * log(phij_);
  cout << "a" << a << endl << flush;

  cout << "id" << id << endl << flush;
  cout << "(xj(pos) + id)" << (xj(pos) + id) << endl << flush;
  cout << "mu(pos) * phij_ * nuj + id" << mu(pos) * phij_ * nuj + id << endl << flush;
  cout << "log(mu(pos) * phij_ * nuj + id)" << log(mu(pos) * phij_ * nuj + id) << endl << flush;
  cout << "antes da sum:" << endl << flush;
  cout << (xj(pos) + id) % log(mu(pos) * phij_ * nuj + id) << endl << flush;
  cout << "depois da sum:" << endl << flush;
  cout << sum((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id)) << endl << flush;
  cout << "separacao" << endl << flush;
  cout << -sum((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id)) << endl << flush;
  cout << "-------------" << endl << flush;
  double b = -sum((xj(pos) + id) % log(mu(pos) * phij_ * nuj + id));
  cout << "-------------" << endl << flush;
  cout << "b" << b << endl << flush;

  cout << "ponto 5" << endl << flush;
  double kappaj = phij.get_kappa();

  cout << "ponto 6" << endl << flush;
  return a + b - kappaj * kappaj / (2 * kappa_var);
}

double nuj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
                       const Phij &phij, double nuj, double sj, double theta) {
  auto pos = find(delta > 0);
  auto zer = find(delta == 0);

  auto q = xj.size();
  auto id = 1 / delta(pos);

  double a = (sum(xj) + q / theta - q) * log(nuj);

  double b0 = -sum((xj(pos) + id) % log(mu(pos) * phij.get() * nuj + id));
  double b1 = -nuj * sum(mu(zer) + 1 / (theta * sj));

  return a + b0 + b1;
}

// double sj_post_loglik(const vec &xj, const vec &mu, const vec &delta,
//                       const Phij &phij, double nuj, const Sj &sj, double theta) {
//   double p = sj.gamma_shape() - 1 / theta;
//   double b = 2 * sj.gamma_rate();
//
//   return RgigDouble(p, 2 * nuj / theta, b);
// }
