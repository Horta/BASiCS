#ifndef PHI_H
#define PHI_H

#include <armadillo>
#include <cmath>
#include <list>

using arma::vec;
using std::initializer_list;

class Phij {
private:
  double kappa;
  double scale;
  double denominator;

public:
  Phij(double kappa, double scale, double denominator)
      : kappa(kappa), scale(scale), denominator(denominator) {}

  double get(void) const { return scale * std::exp(kappa) / denominator; }
  double get_kappa(void) const { return kappa; }
};

class Phi {
private:
  vec kappa;

public:
  Phi(size_t n) : kappa(n - 1) {}
  Phi(initializer_list<double> kappa_list) : kappa(kappa_list) {}

  vec get(void) const;
  void set(const vec &kappa) { this->kappa = kappa; }

  Phij phij(size_t j) const;
};

#endif
