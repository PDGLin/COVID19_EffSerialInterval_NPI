# logLikelihood function

logLikelihood <- function(parms, data, likelihood_options) {   # fixedVariables
  x = data;
  
  is_Normal_fit = likelihood_options$is_Normal;
  is_Gumbel_fit  = likelihood_options$is_Gumbel;
  is_Logistic_fit = likelihood_options$is_Logistic;
  is_LogNormal_fit = likelihood_options$is_LogNormal;
  
  is_transPair_fit = likelihood_options$is_transPair_fit;

  if (is_Normal_fit) { ## Normal distribution ---------------
    miu   = parms[1];
    sigma = parms[2];
    LogL = -log(sigma * sqrt(2*pi)) - 0.5 * ( (x - miu) / sigma )^2;
  
  } else if (is_Gumbel_fit) { ## Gumbel distribution ---------
    miu = parms[1];
    beta = parms[2];
    z = (x - miu) / beta;
    LogL = log(1/beta) - z - exp(- z);
  
  } else if (is_Logistic_fit) {
    miu = parms[1];
    s = parms[2];
    y = exp(- (x - miu) / s);
    LogL = log(y) - log(s) - 2 * log(1 + y);
  
  } else if (is_LogNormal_fit) {
    miu   = parms[1];
    sigma = parms[2];
    LogL = dlnorm(x, miu, sigma, log = T);
  
  
  } else if (is_transPair_fit) {
    
    P_SI_vec = mapply(
      function (D_i, SI, eta, gamma)
        1.0 / D_i * exp(- eta * D_i) * (exp(gamma * D_i) - 1) * exp(- gamma * SI),
      data$lag_isol_onset_primaryCase,                        # D_i
      data[["p_Date of onset"]] - data[["s_Date of onset"]],  # SI
      MoreArgs = list(
        "eta"   = parms[1],
        "gamma" = parms[2]
        ),
      SIMPLIFY = "vector"
    )
    
    LogL = log(P_SI_vec);
    
  } else {
    # Log normal likelihood ------------------------------
    # miu   = parms[1];
    # sigma = parms[2];
    # logL =  log(1/sigma/x/sqrt(2*pi))-0.5*((log(x)-miu)/sigma)^2;
  }
   
  return(sum(LogL));
}


# logLikelihood <- function(parameters,epiCurveData)
# {
#   PopSize = 100000;
#   durInf = 2.5;   # Mean duration of infectiousness
#   
#   R0 = parameters[1];
#   p = parameters[2];
#   I0 = parameters[3];
#   
#   S0 = PopSize - I0;
#   tmax = length(epiCurveData);
#   initialConditions = c(S0,I0);
#   
#   dailyIncidence = SIR(R0,durInf,tmax,PopSize,initialConditions);
#   
#   X = epiCurveData;
#   N = dailyIncidence;
# 
#   if (any( N<X ))
#   {
#     funval = -1e100;
#   } else
#   {
#     funval = lgamma(N+1)-lgamma(X+1)-lgamma(N-X+1)+ X*log(p)+(N-X)*log(1-p);
#     funval = sum(funval);
#   }
#   
# 	return(funval);
# }