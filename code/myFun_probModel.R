#### Fig.S5: Probabilistic model of a transmission pair  ####

myFun_probModel <- function(C_i, K, theta, is_isolateInfector_beforeInfecteeOnset, D_i) {
  ## shape, rate of gamma distribution of infectiousness -------
  gpar = c(K, 1/theta);  
  
  # incubation period (Li et al NEJM 2020) --------
  # lognormal mean = 5.2; 95% CI = c(4.1, 7.0)
  ln.par1 = 1.434065
  ln.par2 = 0.6612
  lnpar = c(ln.par1, ln.par2)
  
  #### compute PDF of serial interval  ####
  
  ## used for setting integration points
  dtar = 0.05  
  dtau = dtar
  dSI = dtar
  
  SI.ij = seq(-20, 30, by=dSI)  # range of serial interval to integrate
  
  df_SI = tibble(
    SI=rep(SI.ij, length(D_i)),
    Di=rep(D_i, each=length(SI.ij))
    )
  
  if ( !is_isolateInfector_beforeInfecteeOnset ) {   # infector can be isolated after onset time of infectee, allowing serial interval < infector's isolation delay
    pSI_byDi_vec = mapply(
      function(SI, Di, Ci, gpar, lnpar, dtau) {
        tau_vec = seq(0, Ci + Di, dtau);             # infectee's time of infection
        E_vec = SI + Ci - tau_vec;                   # infectee's incubation period
        f.beta = dgamma(tau_vec, gpar[1], gpar[2]);  # gamma distribution of infector's infectiousness
        f.pE = dlnorm(E_vec, lnpar[1], lnpar[2]);    # lognormal distribution of infectee's incubation period
        pSI = sum( f.beta * f.pE * dtau );           # distribution of serial interval
        return(pSI)
      },
      df_SI$SI,
      df_SI$Di,
      MoreArgs = list(
        "Ci" = C_i,
        "gpar"  = gpar,
        "lnpar" = lnpar,
        "dtau" = dtau     # unit of integral 
      )
    )
  } else {                                            
    pSI_byDi_vec = mapply(
      function(SI, Di, Ci, gpar, lnpar, dtau) {
        if (SI < Di) {
          pSI = 0;                                     # becoz infector is assumed to isolate before infectee's onset time
        } else {
          tau_vec = seq(0, Ci + Di, dtau);             # infectee's time of infection
          E_vec = SI + Ci - tau_vec;                   # infectee's incubation period
          f.beta = dgamma(tau_vec, gpar[1], gpar[2]);  # gamma distribution of infector's infectiousness
          f.pE = dlnorm(E_vec, lnpar[1], lnpar[2]);    # lognormal distribution of infectee's incubation period
          pSI = sum( f.beta * f.pE * dtau );           # distribution of serial interval
        }
        return(pSI)
      },
      df_SI$SI,
      df_SI$Di,
      MoreArgs = list(
        "Ci" = C_i,
        "gpar"  = gpar,
        "lnpar" = lnpar,
        "dtau" = dtau     # unit of integral 
      )
    )
  }
  
  df_SI <- df_SI %>% mutate( pdf = pSI_byDi_vec )
  
  ## normalize pdf ------------------
  B = sapply(
    D_i, 
    function(tarDi, df_SI, dSI) {
      pdfvec = df_SI %>% filter(Di == tarDi) %>% pull("pdf");
      return(pdfvec / sum(pdfvec * dSI))
    },
    df_SI = df_SI,
    dSI = dSI,
    simplify = "vector"
    ) 
  
  df_SI$pdf <- as.vector(B)
  
  ## compute CDF of serial intervals ------------------
  cdf_vec <- mapply(
    function(tarSI, tarDi, df_SI, dSI) 
      sum( df_SI %>% filter(Di == tarDi & SI <= tarSI) %>% pull("pdf") ) * dSI,
    df_SI$SI,
    df_SI$Di,
    MoreArgs = list(
      "df_SI" = df_SI,
      "dSI" = dSI
      )
    )
  
  df_SI <- df_SI %>% mutate( cdf = cdf_vec )
  
  ## mean, median, 95% CrI of serial interval given isolation delay ---------------
  SI_stat_byDi = sapply(
    D_i, 
    function (tarDi, df_SI, dSI) {
      A = df_SI %>% filter(Di == tarDi);
      meanSI   = sum(A$SI * A$pdf * dSI);
      cdfvec = A$cdf;
      medianSI = A$SI[ max(which(cdfvec <= 0.5)) ];
      SI_LB = A$SI[ max(which(cdfvec <= 0.25)) ];
      SI_UB = A$SI[ max(which(cdfvec <= 0.75)) ];
      return( c(meanSI, medianSI, SI_LB, SI_UB) )
    },
    df_SI = df_SI,
    dSI = dSI
  )
  
  df_stat = tibble(
    Di = D_i,
    SI_mean = SI_stat_byDi[1, ],
    SI_median = SI_stat_byDi[2, ],
    LB = SI_stat_byDi[3, ],
    UB = SI_stat_byDi[4, ]
    )
  
  return(df_stat)
}


