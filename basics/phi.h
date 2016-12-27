#ifndef PHI_H
#define PHI_H

#include <armadillo>

using arma::vec;

class Phi {
private:
  vec kappa;
public:
  Phi(size_t n) : kappa(n-1) {}

  vec get(void);
  void set(const vec& kappa) { this->kappa = kappa; }
};

#endif
