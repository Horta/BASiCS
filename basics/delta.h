#ifndef DELTA_H
#define DELTA_H

#include "prior.h"

#include <armadillo>
#include <cmath>
#include <list>

using arma::vec;
using std::initializer_list;

class Delta {
private:
  vec delta;
  Prior prior;

public:
  Delta(size_t n) : delta(n), prior(Gamma(1, 1)) {}
  Delta(initializer_list<double> delta_list)
      : delta(delta_list), prior(Gamma(1, 1)) {}

  const vec& get(void) const { return delta; }
  vec& get(void) { return delta; }
};

class Deltai {
private:
  const Delta& delta;
  size_t index;

public:
  Deltai(const Delta& delta, size_t index) : delta(delta), index(index) {}

  double get(void) const { return delta.get()(index); }
  // double get_kappa(void) const { return kappa; }
};

#endif
