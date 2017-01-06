#ifndef PRIOR_H
#define PRIOR_H

class Prior {};

class Gamma : public Prior {
public:
  double shape;
  double rate;

  Gamma(double shape, double rate) : shape(shape), rate(rate) {}
};

#endif
