#include <cmath>
#include <limits>
#include <random>

static double Z0TOL = std::numeric_limits<double>::epsilon();
static double Z1TOL = std::sqrt(Z0TOL);

double rgamma(double shape, double scale) {
  using std::gamma_distribution;
  static std::default_random_engine generator;

  return gamma_distribution<double>(shape, scale)(generator);
}

double rnormal(void) {
  using std::normal_distribution;
  static std::default_random_engine generator;

  return normal_distribution<double>()(generator);
}

double runiform(void) {
  using std::uniform_real_distribution;
  static std::default_random_engine generator;

  return uniform_real_distribution<double>()(generator);
}

/* Auxiliary function taken from rgig.c */
double gig_y_gfn(double y, double m, double beta, double lambda)
{
/*
 * gig_y_gfn:
 *
 * evaluate the function that we need to find the root
 * of in order to construct the optimal rejection
 * sampler to obtain GIG samples
 */
  double y2, g;
  y2 = y * y;
  g = 0.5 * beta * y2 * y;
  g -= y2 * (0.5 * beta * m + lambda + 1.0);
  g += y * ((lambda - 1.0) * m - 0.5 * beta) + 0.5 * beta * m;
  return(g);
}

/* Auxiliary function taken from rgig.c */
double rinvgauss(const double mu, const double lambda) {
  /*
   * rinvgauss:
   *
   * Michael/Schucany/Haas Method for generating Inverse Gaussian
   * random variable with mean mu and shape parameter lambda,
   * as given in Gentle's book on page 193
   */
  double u, y, x1, mu2, l2;

  y = rnormal(); /// ******
  y *= y;
  mu2 = mu * mu;
  l2 = 2.0 * lambda;
  x1 =
      mu + mu2 * y / l2 - (mu / l2) * sqrt(4.0 * mu * lambda * y + mu2 * y * y);

  u = runiform(); // ******
  if (u <= mu / (mu + x1))
    return x1;
  else
    return mu2 / x1;
}

/* Auxiliary function taken from rgig.c */
double zeroin_gig(double ax, double bx,
                  double f(double x, double m, double beta, double lambda),
                  double tol, double m, double beta,
                  double lambda) /* An estimate to the root  */
// double ax;				/* Left border | of the range	*/
// double bx;  				/* Right border| the root is seeked*/
/* Function under investigation	*/
// double (*f)(double x, double m, double beta, double lambda);
// double tol;				/* Acceptable tolerance	*/
// double m;                               /* specific to gig_y_gfn */
// double beta;                            /* specific to gig_y_gfn */
// double lambda;                          /* specific to gig_y_gfn */
{

  double a, b, c; /* Abscissae, descr. see above	*/
  double fa;      /* f(a)				*/
  double fb;      /* f(b)				*/
  double fc;      /* f(c)				*/

  a = ax;
  b = bx;
  fa = (f)(a, m, beta, lambda);
  fb = (f)(b, m, beta, lambda);
  c = a;
  fc = fa;

  for (;;) /* Main iteration loop	*/
  {
    double prev_step = b - a; /* Distance from the last but one*/
                              /* to the last approximation	*/
    double tol_act;           /* Actual tolerance		*/
    double p;                 /* Interpolation step is calcu- */
    double q;                 /* lated in the form p/q; divi- */
                              /* sion operations is delayed   */
    /* until the last moment	*/
    double new_step; /* Step at this iteration       */

    if (fabs(fc) < fabs(fb)) { /* Swap data for b to be the 	*/
      a = b;
      b = c;
      c = a; /* best approximation		*/
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol_act = 2.0 * Z0TOL * fabs(b) + tol / 2.0;
    new_step = (c - b) / 2.0;

    if (fabs(new_step) <= tol_act || fb == (double)0)
      return b; /* Acceptable approx. is found	*/

    /* Decide if the interpolation can be tried	*/
    if (fabs(prev_step) >= tol_act /* If prev_step was large enough*/
        && fabs(fa) > fabs(fb))    /* and was in true direction,	*/
    {                              /* Interpolatiom may be tried	*/
      double t1, cb, t2;
      cb = c - b;
      if (a == c)     /* If we have only two distinct	*/
      {               /* points linear interpolation 	*/
        t1 = fb / fa; /* can only be applied		*/
        p = cb * t1;
        q = 1.0 - t1;
      } else /* Quadric inverse interpolation*/
      {
        q = fa / fc;
        t1 = fb / fc;
        t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > (double)0) /* p was calculated with the op-*/
        q = -q;          /* posite sign; make p positive	*/
      else               /* and assign possible minus to	*/
        p = -p;          /* q				*/

      if (p < (0.75 * cb * q -
               fabs(tol_act * q) / 2.0)     /* If b+p/q falls in [b,c]*/
          && p < fabs(prev_step * q / 2.0)) /* and isn't too large	*/
        new_step = p / q;                   /* it is accepted	*/
                                            /* If p/q is too large then the	*/
                                            /* bissection procedure can 	*/
                                            /* reduce [b,c] range to more	*/
                                            /* extent			*/
    }

    if (fabs(new_step) < tol_act) { /* Adjust the step to be not less*/
      if (new_step > (double)0)     /* than tolerance		*/
        new_step = tol_act;
      else
        new_step = -tol_act;
    }

    a = b;
    fa = fb; /* Save the previous approx.	*/
    b += new_step;
    fb = (*f)(b, m, beta, lambda); /* Do step to a new approxim. */
    if ((fb > 0.0 && fc > 0.0) ||
        (fb < 0.0 && fc < 0.0)) { /* Adjust c for it to have a sign*/
      c = a;
      fc = fa; /* opposite to that of b	*/
    }
  }
}

/* Random sample generator from a Generalized Inverse Gaussian distribution
   (modified version of rgig.c using Rcpp classes)
*/
double RgigDouble(double lambda, double chi, double psi) {
  double samps;
  /* special case which is basically a gamma distribution */
  if ((chi < Z1TOL) & (lambda > 0.0)) {
    // TODO: check if Im using the right parameters
    samps = rgamma(lambda, 2.0 / psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if ((psi < Z1TOL) & (lambda < 0.0)) {
    // TODO: check if Im using the right parameters
    samps = 1.0 / rgamma(0.0 - lambda, 2.0 / chi);
    return samps;
  }

  /* special case which is basically an inverse gaussian distribution */
  if (lambda == -0.5) {
    double alpha;
    alpha = sqrt(chi / psi);
    samps = rinvgauss(alpha, chi);
    return samps;
  }

  /*
   * begin general purpose rgig code, which was basically
   * translated from the R function rgig in the ghyp package v_1.5.2
   */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2,
      Y;
  int need;

  alpha = sqrt(chi / psi);
  beta2 = psi * chi;
  beta = sqrt(psi * chi);
  lm1 = lambda - 1.0;
  lm12 = lm1 * lm1;
  m = (lm1 + sqrt(lm12 + beta2)) / beta;
  m1 = m + 1.0 / m;

  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) {
    upper *= 2.0;
  }

  yM = zeroin_gig(0.0, m, gig_y_gfn, Z1TOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, Z1TOL, m, beta, lambda);

  a = (yP - m) * pow(yP / m, 0.5 * lm1);
  a *= exp(-0.25 * beta * (yP + 1.0 / yP - m1));
  b = (yM - m) * pow(yM / m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1 / yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);

  need = 1;
  while (need) {
    R1 = runiform();
    R2 = runiform();

    Y = m + a * R2 / R1 + b * (1.0 - R2) / R1;
    if (Y > 0.0) {
      if (-log(R1) >= -0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0 / Y) + c) {
        need = 0;
      }
    }
  }
  samps = Y * alpha;
  return (samps);
}
