
myFun_choose_likelihood <- function( tar_likelihood_model ){
  
  likelihood_options = tibble(
    is_Normal = F,
    is_Gumbel = F,
    is_Logistic = F,
    is_LogNormal = F,
    is_transPair_fit = F
    )
  
  if ( strcmp(tar_likelihood_model, "Normal_fit") ) {
    likelihood_options$is_Normal = T;
  } else if ( strcmp(tar_likelihood_model, "Gumbel_fit") ) {
    likelihood_options$is_Gumbel = T;
  } else if ( strcmp(tar_likelihood_model, "Logist_fit") ) {
    likelihood_options$is_Logistic = T;
  } else if ( strcmp(tar_likelihood_model, "LogNormal_fit") ) {
    likelihood_options$is_LogNormal = T;
  } else if ( strcmp(tar_likelihood_model, "transPair_fit") ) {
    likelihood_options$is_transPair_fit = T;
  } else {  }
  
  return( likelihood_options )
}
