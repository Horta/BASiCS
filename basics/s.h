#ifndef S_H
#define S_H

#include <armadillo>
#include <cmath>

using arma::vec;

class Sj {
private:
  double sj;
  double _gamma_shape;
  double _gamma_rate;

public:
  Sj(double sj, double gamma_shape, double gamma_rate)
      : sj(sj), _gamma_shape(gamma_shape), _gamma_rate(gamma_rate) {}

  double get(void) const { return sj; }
  double gamma_shape(void) const { return _gamma_shape; }
  double gamma_rate(void) const { return _gamma_rate; }
};

class S {
private:
  vec s;
  double _gamma_shape;
  double _gamma_rate;

public:
  S(size_t n, double gamma_shape, double gamma_rate)
      : s(n), _gamma_shape(gamma_shape), _gamma_rate(gamma_rate) {}

  vec &get(void) { return s; }
  const vec &get(void) const { return s; }

  Sj sj(size_t j) const;

  double gamma_shape(void) const { return _gamma_shape; }
  double gamma_rate(void) const { return _gamma_rate; }
};

#endif
