## mcmc posterior analyser

MCMC_posterior_analyser <- function(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data) {
  
  quantile_vec = c(0.025, 0.5, 0.975);
  
  is_Gumbel_fit = likelihood_options$is_Gumbel;
  is_Normal_fit = likelihood_options$is_Normal;
  is_Logist_fit = likelihood_options$is_Logistic
  
  chainRecord = chainRecord[(0.4 * numStepsPerParameter):nrow(chainRecord), ];
  pars_posterior_sampler = chainRecord[seq(1, nrow(chainRecord), 10), ];
  
  if (is_Gumbel_fit) {
    posterior_mean = pars_posterior_sampler[, 1] + pars_posterior_sampler[, 2] * 0.5772156649;
    posterior_sd = sqrt( ( pi * pars_posterior_sampler[, 2] )^2 / 6 );
    meanVal_fit = quantile(posterior_mean, quantile_vec);
    sdVal_fit   = quantile(posterior_sd, quantile_vec );
  } else if (is_Normal_fit) {
    meanVal_fit = quantile( pars_posterior_sampler[, 1], quantile_vec );
    sdVal_fit   = quantile( pars_posterior_sampler[, 2], quantile_vec );
  } else if (is_Logist_fit) {
    posterior_mean = pars_posterior_sampler[, 1];
    posterior_sd   = pars_posterior_sampler[, 2] * pi / sqrt(3);
    meanVal_fit = quantile(posterior_mean, quantile_vec );
    sdVal_fit   = quantile(posterior_sd, quantile_vec );
  }
  
  df_pars_estim$meanSI_median <- meanVal_fit[2];
  df_pars_estim$meanSI_LB <- meanVal_fit[1];
  df_pars_estim$meanSI_UB <- meanVal_fit[3];
  df_pars_estim$sdSI_median <- sdVal_fit[2];
  df_pars_estim$sdSI_LB   <- sdVal_fit[1];
  df_pars_estim$sdSI_UB   <- sdVal_fit[3];
  df_pars_estim$fitMedian <- qnorm(0.5,  mean=df_pars_estim$meanSI_median, sd=df_pars_estim$sdSI_median);
  df_pars_estim$fitIQR_LB <- qnorm(0.25, mean=df_pars_estim$meanSI_median, sd=df_pars_estim$sdSI_median);
  df_pars_estim$fitIQR_UB <- qnorm(0.75, mean=df_pars_estim$meanSI_median, sd=df_pars_estim$sdSI_median);
  df_pars_estim$DIC <- Get_DIC(pars_posterior_sampler, data, likelihood_options);
  
  return(df_pars_estim)

}