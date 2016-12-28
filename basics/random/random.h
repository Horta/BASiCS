#ifndef RANDOM_H
#define RANDOM_H

#include <random>

class Random
{
private:
  std::default_random_engine generator;

public:
  Random(unsigned int seed) : generator(seed) {}

  double normal(void);
  double uniform(void);
  double gig(double lambda, double chi, double psi);
  double gamma(double shape, double scale);
};

#endif
