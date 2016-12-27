#include "s.h"

Sj S::sj(size_t j) const { return Sj(s(j), _gamma_shape, _gamma_rate); }
