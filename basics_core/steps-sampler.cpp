/* MCMC sampler
 */
// [[Rcpp::export]]
Rcpp::List HiddenBASiCS_MCMCcpp(
  int N, // Total number of MCMC draws
  int thin, // Thinning period for MCMC chain
  int burn, // Burning period for MCMC chain
  NumericMatrix Counts, // $q \times n$ matrix of expression counts (technical genes at the bottom)
  NumericVector mu0, // Starting value of $\mu=(\mu_1,...,\mu_q)'$ (true mRNA content for technical genes)
  NumericVector delta0, // Starting value of $\delta=(\delta_1,...,\delta_{q_0})'$
  NumericVector phi0, // Starting value of $\phi=(\phi_1,...,\phi_n)$'$
  NumericVector s0, // Starting value of $s=(s_1,...,s_n)$'$
  NumericVector nu0, // Starting value of $\nu=(\nu_1,...,\nu_n)$'$
  double theta0, // Starting value of $\theta$
  double s2mu,
  double adelta, // Shape hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$
  double bdelta, // Rate hyper-parameter of the Gamma($a_{\delta}$,$b_{\delta}$) prior assigned to each $\delta_i$
  NumericVector p_Phi, // Dirichlet hyper-parameter for $\phi / n$
  double as, // Shape hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$ */
  double bs, // Rate hyper-parameter of the Gamma($a_s$,$b_s$) prior assigned to each $s_j$
  double atheta, // Shape hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double btheta, // Rate hyper-parameter of the Gamma($a_{\theta}$,$b_{\theta}$) prior assigned to $\theta$
  double ar, // Optimal acceptance rate for adaptive Metropolis-Hastings updates
  NumericVector LSmu0, // Starting value of adaptive proposal variance of $\mu=(\mu_1,...,\mu_q)'$ (log-scale)
  NumericVector LSdelta0, // Starting value of adaptive proposal variance of $\delta=(\delta_1,...,\delta_{q_0})'$ (log-scale)
  double LSphi0, // Starting value of adaptive proposal precision of $\phi=(\phi_1,...,\phi_n)'$ (log-scale)
  NumericVector LSnu0, // Starting value of adaptive proposal variance of $\nu=(\nu_1,...,\nu_n)'$ (log-scale)
  double LStheta0, // Starting value of adaptive proposal variance of $\theta$ (log-scale)
  NumericVector sumByCellAll, // Sum of expression counts by cell (all genes)
  NumericVector sumByCellBio, // Sum of expression counts by cell (biological genes only)
  NumericVector sumByGeneAll, // Sum of expression counts by gene (all genes)
  NumericVector sumByGeneBio,
  int StoreAdapt,
  int EndAdapt,
  int PrintProgress,
  double s2_delta,
  double prior_delta)
{


// NUMBER OF CELLS, GENES AND STORED DRAWS
  int n = Counts.ncol(); int q = Counts.nrow(); int qbio = delta0.size(); int Naux = N/thin - burn/thin;
  // TRANSFORMATION TO ARMA ELEMENTS
  arma::vec sumByCellAll_arma = as_arma(sumByCellAll); arma::vec sumByCellBio_arma = as_arma(sumByCellBio);
  arma::vec sumByGeneAll_arma = as_arma(sumByGeneAll); arma::vec sumByGeneBio_arma = as_arma(sumByGeneBio);
  arma::mat Counts_arma = as_arma(Counts);
  arma::vec mu0_arma = as_arma(mu0);
  arma::vec phi0_arma = as_arma(phi0);
  // OBJECTS WHERE DRAWS WILL BE STORED
  arma::mat mu = arma::zeros(Naux,q);
  arma::mat delta = arma::zeros(Naux,qbio);
  arma::mat phi = arma::ones(Naux,n);
  arma::mat s = arma::zeros(Naux,n);
  arma::mat nu = arma::zeros(Naux,n);
  arma::vec theta = arma::zeros(Naux);
  arma::mat LSmu;
  arma::mat LSdelta;
  arma::vec LSphi;
  arma::mat LSnu;
  arma::vec LStheta;

  // LOG-PROPOSAL VARIANCES
  if(StoreAdapt == 1)
  {
    LSmu = arma::zeros(Naux,q);
    LSdelta = arma::zeros(Naux,qbio);
    LSphi = arma::ones(Naux);
    LSnu = arma::zeros(Naux,n);
    LStheta = arma::zeros(Naux);
  }
  // ACCEPTANCE RATES FOR ADAPTIVE METROPOLIS-HASTINGS UPDATES
  arma::vec muAccept = arma::zeros(q); arma::vec PmuAux = arma::zeros(q);
  arma::vec deltaAccept = arma::zeros(qbio); arma::vec PdeltaAux = arma::zeros(qbio);
  double phiAccept = 0; double PphiAux = 0;
  arma::vec nuAccept = arma::zeros(n); arma::vec PnuAux = arma::zeros(n);
  double thetaAccept=0; double PthetaAux=0;
  // INITIALIZATION OF VALUES FOR MCMC RUN
  arma::mat muAux = arma::zeros(q,2); muAux.col(0)=mu0_arma; arma::vec LSmuAux = as_arma(LSmu0);
  arma::mat deltaAux = arma::zeros(qbio,2); deltaAux.col(0)=as_arma(delta0); arma::vec LSdeltaAux = as_arma(LSdelta0);
  arma::vec phiAux = phi0_arma; double LSphiAux = LSphi0; Rcpp::List phiAuxList;
  arma::vec sAux = as_arma(s0);
  arma::mat nuAux = arma::zeros(n,2); nuAux.col(0)=as_arma(nu0); arma::vec LSnuAux = as_arma(LSnu0);
  arma::vec thetaAux = arma::zeros(2); thetaAux(0) = theta0; double LSthetaAux = LStheta0;
  // OTHER AUXILIARY QUANTITIES FOR ADAPTIVE METROPOLIS UPDATES
  arma::vec PmuAux0 = arma::zeros(q); arma::vec PdeltaAux0 = arma::zeros(qbio);
  double PphiAux0 = 0;
  arma::vec PnuAux0 = arma::zeros(n); double PthetaAux0=0;

  double SumSpikeInput = sum(mu0_arma(arma::span(qbio,q -1)));



  NumericVector LSphiValuesNV = NumericVector::create(1, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000);
  arma::vec LSphiValues = as_arma(LSphiValuesNV);
  NumericVector LSthetaValuesNV = NumericVector::create(-7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2);
  arma::vec LSthetaValues = as_arma(LSthetaValuesNV);





  // BATCH INITIALIZATION FOR ADAPTIVE METROPOLIS UPDATES (RE-INITIALIZE EVERY 50 ITERATIONS)
  int Ibatch = 0; int i;

  // INITIALIZATION OF PARAMETERS TO RETURN IN UPDATE FUNCTIONS
  arma::vec muUpdateAux = arma::ones(q); muUpdateAux(arma::span(qbio, q-1)) = mu0_arma(arma::span(qbio, q-1));
  arma::vec indQ = arma::zeros(q);


  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "MCMC sampler has been started: " << N << " iterations to go." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;

  // START OF MCMC LOOP
  for (i=0; i<N; i++) {

//    Rcpp::checkUserInterrupt();

    if(i==burn)
    {
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
      Rcpp::Rcout << "End of burn-in period."<< std::endl;
      Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
    }

    Ibatch++;

    // UPDATE OF PHI: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
//    LSphiAux = rUnifDisc(LSphiValues, LSphiValuesNV.size());
    phiAuxList = phiUpdate(phiAux, exp(LSphiAux), Counts_arma, muAux.col(0), deltaAux.col(0),
                           nuAux.col(0), p_Phi, sumByGeneBio_arma, qbio, n);
    phiAux = Rcpp::as<arma::vec>(phiAuxList["phi"]);
    PphiAux += Rcpp::as<double>(phiAuxList["ind"]); if(i>=burn) {phiAccept += Rcpp::as<double>(phiAuxList["ind"]);}

    // UPDATE OF THETA: 1st ELEMENT IS THE UPDATE, 2nd ELEMENT IS THE ACCEPTANCE INDICATOR
//    LSthetaAux = rUnifDisc(LSthetaValues, LSthetaValuesNV.size());
    thetaAux = thetaUpdate(thetaAux(0), exp(LSthetaAux), sAux, nuAux.col(0), atheta, btheta, n);
    PthetaAux += thetaAux(1); if(i>=burn) {thetaAccept += thetaAux(1);}

    // UPDATE OF MU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
//    LSmuAux = rUnifDiscMult(q, LSthetaValues, LSthetaValuesNV.size());
    muAux = muUpdate(muAux.col(0), exp(LSmuAux), Counts_arma, deltaAux.col(0),
                     phiAux, nuAux.col(0), sumByCellBio_arma, s2mu, q, qbio, n,
                     muUpdateAux, indQ);
    PmuAux += muAux.col(1); if(i>=burn) {muAccept += muAux.col(1);}




    // UPDATE OF S
    sAux = sUpdate(sAux, nuAux.col(0), thetaAux(0), as, bs, n);

    // UPDATE OF DELTA: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
//    LSdeltaAux = rUnifDiscMult(qbio, LSthetaValues, LSthetaValuesNV.size());
    deltaAux = deltaUpdate(deltaAux.col(0), exp(LSdeltaAux), Counts_arma, muAux.col(0),
                           nuAux.col(0), phiAux, adelta, bdelta, qbio, n, s2_delta, prior_delta);
    PdeltaAux += deltaAux.col(1); if(i>=burn) {deltaAccept += deltaAux.col(1);}

    // UPDATE OF NU: 1st COLUMN IS THE UPDATE, 2nd COLUMN IS THE ACCEPTANCE INDICATOR
    nuAux = nuUpdate(nuAux.col(0), exp(LSnuAux), Counts_arma, SumSpikeInput,
                     muAux.col(0), deltaAux.col(0), phiAux, sAux, thetaAux(0), sumByGeneAll_arma, q, qbio, n);
    PnuAux += nuAux.col(1); if(i>=burn) {nuAccept += nuAux.col(1);}

    // STOP ADAPTING THE PROPOSAL VARIANCES AFTER EndAdapt ITERATIONS
    if(i < EndAdapt)
    {
      // UPDATE OF PROPOSAL VARIANCES (ONLY EVERY 50 ITERATIONS)
      if(Ibatch==50)
      {
        PmuAux=PmuAux/50; PmuAux = -1+2*arma::conv_to<arma::mat>::from(PmuAux>ar);
//        LSmuAux=LSmuAux+PmuAux*std::min(0.03,1/sqrt(i));
        LSmuAux=LSmuAux+PmuAux*0.1;
        PdeltaAux=PdeltaAux/50; PdeltaAux = -1+2*arma::conv_to<arma::mat>::from(PdeltaAux>ar);
//        LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.03,1/sqrt(i));
        LSdeltaAux=LSdeltaAux+PdeltaAux*0.1;
//        if(i < 1000) {LSdeltaAux=LSdeltaAux+PdeltaAux*0.5;}
//        else {LSdeltaAux=LSdeltaAux+PdeltaAux*std::min(0.01,1/sqrt(i));}
        PphiAux=PphiAux/50; PphiAux = -1+2*(PphiAux>ar);
        LSphiAux=LSphiAux - PphiAux*0.1;
//        LSphiAux=LSphiAux - PphiAux*std::min(0.03,1/sqrt(i));
        PnuAux=PnuAux/50; PnuAux = -1+2*arma::conv_to<arma::mat>::from(PnuAux>ar);
//        LSnuAux=LSnuAux+PnuAux*std::min(0.03,1/sqrt(i));
        LSnuAux=LSnuAux+PnuAux*0.1;
        PthetaAux=PthetaAux/50; PthetaAux = -1+2*(PthetaAux>ar);
//        LSthetaAux=LSthetaAux+PthetaAux*std::min(0.03,1/sqrt(i));
        LSthetaAux=LSthetaAux+PthetaAux*0.1;

        Ibatch = 0;
        PmuAux = PmuAux0; PdeltaAux = PdeltaAux0;
        PphiAux = PphiAux0;
        PnuAux = PnuAux0; PthetaAux = PthetaAux0;
      }

    }

    // STORAGE OF DRAWS
    if((i%thin==0) & (i>=burn))
    {
      mu.row(i/thin - burn/thin) = muAux.col(0).t();
      delta.row(i/thin - burn/thin) = deltaAux.col(0).t();
      phi.row(i/thin - burn/thin) = phiAux.t();
      s.row(i/thin - burn/thin) = sAux.t();
      nu.row(i/thin - burn/thin) = nuAux.col(0).t();
      theta(i/thin - burn/thin) = thetaAux(0);

      if(StoreAdapt == 1)
      {
        LSmu.row(i/thin - burn/thin) = LSmuAux.t();
        LSdelta.row(i/thin - burn/thin) = LSdeltaAux.t();
        LSphi(i/thin - burn/thin) = LSphiAux;
        LSnu.row(i/thin - burn/thin) = LSnuAux.t();
        LStheta(i/thin - burn/thin) = LSthetaAux;
      }
    }

    // PRINT IN CONSOLE SAMPLED VALUES FOR FEW SELECTED PARAMETERS
    if((i%(2*thin) == 0) & (PrintProgress == 1))
    {
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "MCMC iteration " << i << " out of " << N << " has been completed." << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current draws of some selected parameters are displayed below." << std::endl;
        Rcpp::Rcout << "mu (gene 1): " << muAux(0,0) << std::endl;
        Rcpp::Rcout << "delta (gene 1): " << deltaAux(0,0) << std::endl;
        Rcpp::Rcout << "phi (cell 1): " << phiAux(0) << std::endl;
        Rcpp::Rcout << "s (cell 1): " << sAux(0) << std::endl;
        Rcpp::Rcout << "nu (cell 1): " << nuAux(0,0) << std::endl;
        Rcpp::Rcout << "theta: " << thetaAux(0) << std::endl;
        Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
        Rcpp::Rcout << "Current proposal variances for Metropolis Hastings updates (log-scale)." << std::endl;
        Rcpp::Rcout << "LSmu (gene 1): " << LSmuAux(0) << std::endl;
        Rcpp::Rcout << "LSdelta (gene 1): " << LSdeltaAux(0) << std::endl;
        Rcpp::Rcout << "LSphi: " << LSphiAux << std::endl;
        Rcpp::Rcout << "LSnu (cell 1): " << LSnuAux(0) << std::endl;
        Rcpp::Rcout << "LStheta: " << LSthetaAux << std::endl;
    }
  }

//  ProfilerStop();

  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "All " << N << " MCMC iterations have been completed." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;
  // ACCEPTANCE RATE CONSOLE OUTPUT
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << "Please see below a summary of the overall acceptance rates." << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;

  Rcpp::Rcout << "Minimum acceptance rate among mu[i]'s: " << min(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among mu[i]'s: " << mean(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among mu[i]'s: " << max(muAccept(arma::span(0, qbio - 1))/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among delta[i]'s: " << min(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among delta[i]'s: " << mean(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among delta[i]'s: " << max(deltaAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for phi (joint): " << phiAccept/(N-burn) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Minimum acceptance rate among nu[j]'s: " << min(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Average acceptance rate among nu[j]'s: " << mean(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << "Maximum acceptance rate among nu[j]'s: " << max(nuAccept/(N-burn)) << std::endl;
  Rcpp::Rcout << " " << std::endl;
  Rcpp::Rcout << "Acceptance rate for theta: " << thetaAccept/(N-burn) << std::endl;
  Rcpp::Rcout << "--------------------------------------------------------------------" << std::endl;
  Rcpp::Rcout << " " << std::endl;


  if(StoreAdapt == 1)
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu")=mu,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("phi")=phi,
           Rcpp::Named("s")=s,
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta,
           Rcpp::Named("ls.mu")=LSmu,
           Rcpp::Named("ls.delta")=LSdelta,
           Rcpp::Named("ls.phi")=LSphi,
           Rcpp::Named("ls.nu")=LSnu,
           Rcpp::Named("ls.theta")=LStheta));
  }

  else
  {
    // OUTPUT (AS A LIST)
    return(Rcpp::List::create(
           Rcpp::Named("mu")=mu,
           Rcpp::Named("delta")=delta,
           Rcpp::Named("phi")=phi,
           Rcpp::Named("s")=s,
           Rcpp::Named("nu")=nu,
           Rcpp::Named("theta")=theta));
  }

}
