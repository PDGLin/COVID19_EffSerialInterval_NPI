#######################################################################
## Function of Markov chain Monte Carlo (MCMC)                       ##
## Ver 1.0 (2017/05) by Prof. Joseph T. Wu, Mr. CK Lam, Dr. Lin Wang ##
## Ver 2.0 (2020/05) by Dr. Lin Wang                                 ##
## Email: lw660@cam.ac.uk                                            ##
## Web: http://pdg.gen.cam.ac.uk/                                    ##
#######################################################################

MCMC <- function(data_to_fit, numStepsPerParameter, likelihood_options)
{
  MCMC_pars_ls <- SetPars_MCMC( numStepsPerParameter )
  startingPoint = MCMC_pars_ls$startingPoint;
  LB = MCMC_pars_ls$LB;
  UB = MCMC_pars_ls$UB;
  minProbAccept = MCMC_pars_ls$minProbAccept;
  maxProbAccept = MCMC_pars_ls$maxProbAccept;

  minIterBeforeRestart = 1000 * length(LB);
	
  numStepsBetweenDisplay = 500;
  storageSize = min(numStepsPerParameter, 1000);
  probAccept = (minProbAccept + maxProbAccept) / 2;
  oldProbAccept = array(-1, c(1, length(LB)));

	stepSTD = rep(0.1, length(LB)); # kronecker(matrix(1,1,length(LB)),0.1);
	delta = 30;

	chainRecord = array(0, c(numStepsPerParameter * length(LB), length(LB) + 2));
	chainGood = -1;

	bestParameters = startingPoint;
	bestLogLikelihood = logLikelihood(bestParameters, data_to_fit, likelihood_options);
	
	while (chainGood == -1){
	  x = bestParameters;
	  xLogLikelihood = logLikelihood(x, data_to_fit, likelihood_options);
	  xTarget = xLogLikelihood + logPrior(x);
    
	  chainRecord[1, ] = c(x, xLogLikelihood, 1);
		numAccept = array(1, c(1, length(x)));
		numSteps  = array(1, c(1, length(x)));
		
		stepSizeStorage = array(0, c(storageSize, length(x)));
		for (para in 1:length(LB)) {
		  if (probAccept[para] < minProbAccept[para]) {
		    p_delta = minProbAccept[para] - probAccept[para];
			  scalingFactor = 1/(1 + p_delta*delta);
			  stepSTD[para] = scalingFactor*stepSTD[para]; # need revision according to matlab code
			}
			if (probAccept[para] > maxProbAccept[para]) {
			  p_delta = probAccept[para] - maxProbAccept[para];
			  scalingFactor = 1 + p_delta*delta;
				stepSTD[para] = scalingFactor*stepSTD[para]; # need revision
			}
  	}
    	
    paraUpdate = 1;
		tic();
		totalSecElapsed = 0;
		probStorage = runif(storageSize*length(x), 0, 1);
		for (para in 1:length(x)) {
			stepSizeStorage[, para] = rnorm(storageSize, 0, stepSTD[para]);
		}
		numTotalSteps = 1;
	
		for (iter in 2:(numStepsPerParameter*length(x))) {
			if (numSteps[paraUpdate]%%storageSize==0) {
				stepSizeStorage[,paraUpdate] = rnorm(storageSize, 0, stepSTD[paraUpdate]);
			}
			if (iter%%storageSize==0) {
				probStorage = runif(storageSize, 0, 1);
			}
      
		  # next step in the random walk
			y = x;
			
			y[paraUpdate] = x[paraUpdate] + stepSizeStorage[1+(numSteps[paraUpdate]-1)%%storageSize,paraUpdate];
			
			while (y[paraUpdate] < LB[paraUpdate] || y[paraUpdate] > UB[paraUpdate]) {
				if (y[paraUpdate]<LB[paraUpdate]) {
					y[paraUpdate] = LB[paraUpdate] + (LB[paraUpdate]-y[paraUpdate]);
				} else {
					y[paraUpdate] = UB[paraUpdate] - (y[paraUpdate]-UB[paraUpdate]);
				}
			}
			
			yLogLikelihood = logLikelihood(y, data_to_fit, likelihood_options);
			yTarget = yLogLikelihood + logPrior(y);
			
			if (probStorage[1+((iter-1)%%storageSize)] < min(1,exp(yTarget-xTarget))) {
				x = y;
				xLogLikelihood = yLogLikelihood;
				xTarget = yTarget;
				numAccept[paraUpdate] = numAccept[paraUpdate]+1;
			}
			
			if (all(xLogLikelihood>bestLogLikelihood)) {
				bestLogLikelihood = xLogLikelihood;
				bestParameters = x;
			}
			
			numSteps[paraUpdate] = numSteps[paraUpdate]+1;
			
			chainRecord[iter,] = c(x,xLogLikelihood,probAccept[paraUpdate]);
        	
			probAccept[paraUpdate] = numAccept[paraUpdate]/numSteps[paraUpdate];
			
			if (probAccept[paraUpdate]>minProbAccept[paraUpdate] && probAccept[paraUpdate]<=maxProbAccept[paraUpdate])
			{
				chainGood = 1;
			}else
			{
				chainGood = -1;
			
				if (iter>minIterBeforeRestart)
				{
					print('');
					print('Restarting MCMC because P(Accept) is out of range...');
					print(cat('Parameter ', toString(paraUpdate), ': no. of steps = ', toString(numSteps[paraUpdate]), ', P(Accept) = ', toString(probAccept[paraUpdate]), ', stepSTD = ', toString(stepSTD[paraUpdate])));
				
					plot.new();
					
					# plot graph
					rsize = dim(chainRecord)[1];
					csize = dim(chainRecord)[2];
					numFigRow = floor(sqrt(length(LB)));
					numFigCol = min(length(LB), ceiling(length(LB)/numFigRow));
					numFigRow = numFigRow+1;
					
					par(mfrow=c(numFigRow,numFigCol));
					
					for (ii in 1:length(LB))
					{                 
						plot(chainRecord[1:iter,ii],type="l",xlab="",ylab=cat('Parameter ',toString(ii)));
					}
					
					plot(chainRecord[1:iter,1+length(LB)],type="l",xlab="",ylab="logLikelihood");
					
					barplot(probAccept,names.arg=c(1:length(probAccept)),ylim=c(0,1),xlab="Parameter",ylab="P(Accept)");
					
					Sys.sleep(1);
					break;
				}
			}
			
			paraUpdate = paraUpdate+1;
			if (paraUpdate>length(LB))
			{
				paraUpdate = 1;
			}
			
			totalSecElapsed = totalSecElapsed + toc(echo=FALSE);
			
			tic();
			
			if (paraUpdate==1 && numSteps[paraUpdate]%%numStepsBetweenDisplay==0)
			{
				tt = cat('MCMC ',toString(round(numSteps[paraUpdate]/numStepsPerParameter*100,digits=2)),'% complete. Time elapsed = ',					toString(round(totalSecElapsed/60,digits=2)),' min. Remaining time = ',toString(round(totalSecElapsed/numSteps[paraUpdate]*(numStepsPerParameter-numSteps[paraUpdate])/60,digits=2)),' min.');
				print(tt);
			}
        }
    }
	
	#return(chainRecord,probAccept);
	return(list("chainRecord"=chainRecord,"probAccept"=probAccept));
}
