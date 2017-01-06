#ifndef DELTA_H
#define DELTA_H

#include "prior.h"

#include <armadillo>
#include <cmath>
#include <list>

using arma::vec;
using std::initializer_list;

class Deltai;

class Delta {
private:
  vec delta;
  Prior prior;

public:
  Delta(size_t n) : delta(n), prior(Gamma(1, 1)) {}
  Delta(initializer_list<double> delta_list)
      : delta(delta_list), prior(Gamma(1, 1)) {}

  const vec &get(void) const { return delta; }
  vec &get(void) { return delta; }

  Deltai deltai(size_t i) const;
};

class Deltai {
private:
  const Delta &delta;
  size_t index;

public:
  Deltai(const Delta &delta, size_t index) : delta(delta), index(index) {}

  double get(void) const { return delta.get()(index); }
};

Deltai Delta::deltai(size_t i) const { return Deltai(*this, i); }

#endif
