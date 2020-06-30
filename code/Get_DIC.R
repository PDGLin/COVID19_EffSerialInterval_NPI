
Get_DIC <- function( pars_posterior_sampler, data, likelihood_options) {
  n_sampler = dim(pars_posterior_sampler)[1];
  loglikelihood_posterior = rep(0, n_sampler)
  
  for (idx_like in 1:n_sampler) {
    loglikelihood_posterior[idx_like] = logLikelihood( pars_posterior_sampler[idx_like, ], data, likelihood_options);
  }
  
  # Posterior Mean
  posteriorMean = colMeans(pars_posterior_sampler); 
  D_expect_theta = -2 * logLikelihood(posteriorMean, data, likelihood_options);
  
  D_theta_vec = -2 * loglikelihood_posterior;
  p_D = mean(D_theta_vec) - D_expect_theta;
  
  DIC = p_D + mean(D_theta_vec);
  
  return(DIC)
}