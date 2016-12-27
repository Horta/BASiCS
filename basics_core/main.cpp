/* Metropolis-Hastings updates of mu
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdate(
  arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& sum_bycell_bio, /* Sum of expression counts by cell (biological genes only) */
  double const& s2_mu,
  int const& q,
  int const& q_bio,
  int const& n,
  arma::vec & mu,
  arma::vec & ind)
{
}

/* Metropolis-Hastings updates of mu
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat deltaUpdate(
  arma::vec const& delta0, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& prop_var,  /* Current value of the proposal variances for $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_delta, /* Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  double const& b_delta, /* Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  int const& q_bio,
  int const& n,
  double const& s2_delta,
  double const& prior_delta)
{
}

/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::vec sUpdate(
  arma::vec const& s0, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& theta, /* Current value of $\theta$ */
  double const& as, /* Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double const& bs, /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  int const& n)
{
}

/* Auxiliary function required for some of the Metropolis-Hastings updates of nu */
arma::mat UpdateAux_nuTrick(
  arma::vec const& nu, /* Auxiliary variable (function of the current and proposed value of nu) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi) /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
{
}

/* Metropolis-Hastings updates of nu
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat nuUpdate(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  double const& SumSpikeInput,
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  double const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q,
  int const& q_bio,
  int const& n)
{
}

/* Metropolis-Hastings updates of theta
 */
arma::vec thetaUpdate(
  double const& theta0, /* Current value of $\theta$ */
  double const& prop_var, /* Current value of the proposal variances for $\theta$ */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_theta, /* Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  double const& b_theta, /* Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  int const& n)
  {
}

/* Metropolis-Hastings updates of phi
 * Joint updates using Dirichlet proposals
 */
Rcpp::List phiUpdate(
  arma::vec const& phi0, // Current value of $\phi=(\phi_1,...,\phi_n)'$
  double const& prop_var, // Current value of the proposal precision
  arma::mat const& Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom)
  arma::vec const& mu, // Current value of $\mu=(\mu_1,...,\mu_q)'$
  arma::vec const& delta, // Current value of $\delta=(\delta_1,...,\delta_{q_0})'$
  arma::vec const& nu, // Current value of $\nu=(\nu_1,...,\nu_n)'$
  arma::vec const& p_phi, // Dirichlet hyper-parameter of the prior for $\phi / n$
  arma::vec const& sum_bygene_bio, // Sum of expression counts by gene (biological genes only)
  int const& q_bio, // Number of biological genes
  int const& n) // Total number of cells
{
}






/* Draws for cell-specific normalising constants s[j]
 * Metropolis-Hastings updates are not required as full conditionals have a closed form (Generalized Inverse Gaussian)
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::vec sUpdateBatch(
  arma::vec const& s0_arma, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu_arma, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& theta, /* Current value of $\theta$ */
  double const& as, /* Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double const& bs, /* Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  arma::mat const& BatchDesign,
  int const& n)
{

}

/* Metropolis-Hastings updates of nu
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat nuUpdateBatch(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  double const& SumSpikeInput,
  arma::mat const& BatchDesign,
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q,
  int const& q_bio,
  int const& n)
{

}

/* Metropolis-Hastings updates of theta
 */
arma::mat thetaUpdateBatch(
  arma::vec const& theta0, /* Current value of $\theta$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\theta$ */
  arma::mat const& BatchDesign,
  arma::vec const& s, /* Current value of $s=(s_1,...,s_n)$' */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_theta, /* Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  double const& b_theta, /* Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$ */
  int const& n,
  int const& nBatch)
  {

}



/* Metropolis-Hastings updates of mu
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat muUpdateNoSpikes(
    arma::vec const& mu0, /* Current value of $\mu=(\mu_1,...,\mu_q_0)'$ */
    arma::vec const& prop_var, /* Current value of the proposal variances for $\mu_{bio}=(\mu_1,...,\mu_{q_0})'$ */
    double const& Constrain,
    arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
    arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
    arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
    arma::vec const& sum_bycell_all, /* Sum of expression counts by cell (biological genes only) */
    double const& s2_mu,
    int const& q_bio,
    int const& n,
    arma::vec & mu,
    arma::vec & ind,
    int const& RefGene,
    arma::uvec const& ConstrainGene,
    arma::uvec const& NotConstrainGene,
    int const& ConstrainType)
{

}

/* Metropolis-Hastings updates of delta
 * Updates are implemented simulateaneously for all biological genes, significantly reducing the computational burden.
 */
arma::mat deltaUpdateNoSpikes(
  arma::vec const& delta0, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& prop_var,  /* Current value of the proposal variances for $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& nu, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  double const& a_delta, /* Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  double const& b_delta, /* Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$ */
  int const& q_bio,
  int const& n,
  double const& s2_delta,
  double const& prior_delta)
{

}

/* Metropolis-Hastings updates of nu
 * Updates are implemented simulateaneously for all cells, significantly reducing the computational burden.
 */
arma::mat nuUpdateNoSpikes(
  arma::vec const& nu0, /* Current value of $\nu=(\nu_1,...,\nu_n)'$ */
  arma::vec const& prop_var, /* Current value of the proposal variances for $\nu=(\nu_1,...,\nu_n)'$ */
  arma::mat const& Counts, /* $q \times n$ matrix of expression counts (technical genes at the bottom) */
  arma::mat const& BatchDesign,
  arma::vec const& mu, /* Current value of $\mu=(\mu_1,...,\mu_q)'$ */
  arma::vec const& delta, /* Current value of $\delta=(\delta_1,...,\delta_{q_0})'$ */
  arma::vec const& phi, /* Current value of $\phi=(\phi_1,...,\phi_n)$' */
  arma::vec const& theta, /* Current value of $\theta$' */
  arma::vec const& sum_bygene_all, /* Sum of expression counts by gene (all genes) */
  int const& q_bio,
  int const& n,
  arma::vec const& BatchInfo,
  arma::vec const& BatchIds,
  int const& nBatch,
  arma::vec const& BatchSizes,
  arma::vec const& BatchOffSet)
{

}
