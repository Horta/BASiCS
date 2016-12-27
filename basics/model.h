#include <armadillo>

using arma::vec;

class Model {
public:
  vec &mu();
  vec &delta();
  vec &phi();
  vec &nu();
  vec &theta();
  vec &s();
};
