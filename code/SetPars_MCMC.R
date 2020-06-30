
SetPars_MCMC <- function( numStepsPerParameter ) {
  # To Do: Prepare pars 
  startingPoint = c(1, 1);  # Initial guess values for pars
  LB = c(0, 0) + 1e-6;      # Lower bound for pars
  UB = c(100, 100) - 1e-6;  # Upper bound for pars
  minProbAccept = rep(0.3, length(LB)); # Min acceptance probability of candidates
  maxProbAccept = rep(0.7, length(LB)); # Max acceptance probability of candidates
  
  pars_ls = list(
    "startingPoint" = startingPoint,  
    "LB" = LB,
    "UB" = UB,
    "minProbAccept" = minProbAccept,
    "maxProbAccept" = maxProbAccept
    )
  
  df_pars_estim_full = tibble(
    phase = integer(),
    stratification = integer(),
    onsetDay_LB = integer(), onsetDay_UB = integer(),
    sample_size = integer(),
    meanSI_median = double(), meanSI_LB = double(), meanSI_UB = double(),
    sdSI_median   = double(),   sdSI_LB = double(),   sdSI_UB = double(),
    fitMedian = double(), fitIQR_LB = double(), fitIQR_UB = double(),
    DIC = double()
  )
  
  df_pars_estim = tibble(
    phase = 0,
    stratification = 0,
    onsetDay_LB = 0,
    onsetDay_UB = 0,
    sample_size = 0,
    meanSI_median = 0,
    meanSI_LB = 0,
    meanSI_UB = 0,
    sdSI_median = 0,
    sdSI_LB = 0,
    sdSI_UB = 0,
    fitMedian = 0,
    fitIQR_LB = 0,
    fitIQR_UB = 0,
    DIC = 0
  )
  
  pars_ls[["df_pars_estim_full"]] <- df_pars_estim_full
  pars_ls[["df_pars_estim"]] <- df_pars_estim
  
  return(pars_ls)
  
}


