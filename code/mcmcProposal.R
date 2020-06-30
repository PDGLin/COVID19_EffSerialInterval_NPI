# Proposed step for MCMC
# Proposed new step based on the current status

mcmcProposal <- function(currParm, parmStepSize, parLB, parUB) {
  newParm = rep(NA, length(currParm));
  # While-loop: generate new proposed step 
  # until the new step meets the parameter range constraints 
  while(is.na(newParm) || newParm < parLB || newParm > parUB) {
    # Random walk 
    newParm = currParm + parmStepSize * runif(length(currParm), min=-0.5, max=0.5);
    if(newParm < parLB || newParm > parUB) {
      while(newParm < parLB || newParm>parUB) {
        newParm[newParm < parLB] = parLB[newParm < parLB] + (parLB[newParm < parLB] - newParm[newParm < parLB]);
        newParm[newParm > parUB] = parUB[newParm > parUB] - (newParm[newParm > parUB] - parUB[newParm > parUB]);
      }
    }
  }
  return(newParm);
}