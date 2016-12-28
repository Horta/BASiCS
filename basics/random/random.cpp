#include "random.h"

double Random::normal(void)
{
  return std::normal_distribution<double>()(generator);
}

double Random::uniform(void)
{
  return std::uniform_real_distribution<double>()(generator);
}

double Random::gig(double lambda, double chi, double psi)
{
  
}

double Random::gamma(double shape, double scale)
{
  return std::gamma_distribution<double>(shape, scale)(generator);
}
