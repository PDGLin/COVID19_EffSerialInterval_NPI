#####################################################################################################
## Fig.S6, Simulation of serial intervals with isolation delay from 10 to 0 days                 ####
## Wrote by Dr. Eric H.Y. Lau                                                                    ####
## Email: ehylau@hku.hk                                                                          ####          
## Web: https://sph.hku.hk/en/about-us/faculty-and-staff/academic-related-staff/lau,-ho-yin-eric ####
#####################################################################################################

rm(list=ls())

# Incubation period with mean 5.2 and SD 3.9 days (lognormal) (Li et al., NEJM 2020)
ln.par1 <- 1.434065
ln.par2 <- 0.6612

## Estimated serial intervals in the pre-peak period in our analysis (Ali et al.)
## Normal distribution: mean 7.8, 95% CI=7.0-8.6
si.norm.se <- (8.6-7.0)/(2*1.96)
si.norm.par <- c(7.8, si.norm.se)

## Estimated serial intervals from Zhang et al. Lancet Inf Dis 2020
## Gamma distribution: mean 5.1, 95% CI=1.3-11.6
## Gamma parameters: 3.63, 0.71
si.gamma.par <- c(3.63, 0.71)

sim.si <- function(R0, onset.iso){
# Under condition with no population immunity
n = 50000
data = data.frame(id=1:n, inf.t=NA, infector=NA, ons.t=NA)
n.inf0 = 5000 # initial no. infectious cases allowing for low R0
data$inf.t[1:n.inf0] = 0

## Approximate the generation time by assumed serial interval

counter=n.inf0
# Cases with isolation: cannot infect others after isolation, 
# allowing for a delay from onset to isolation

for (i in 1:10000){ # simulate more than 500 transmission events
  data$ons.t[i] <- data$inf.t[i] + rlnorm(1, ln.par1, ln.par2)
  n.infect <- rpois(1,R0)
  if (n.infect>0){

# Two scenarios (modify as needed):
# Scenario (1) based on our estimate at the pre-peak period (mean = 7.8d)
    t.infect <- data$inf.t[i] + sort(rnorm(n.infect,si.norm.par[1],si.norm.par[2]))
# Scenario (2) based on Zhang et al. (mean = 5.1d)
#    t.infect <- data$inf.t[i] + sort(rgamma(n.infect,si.gamma.par[1],si.gamma.par[2]))

    n.infect <- n.infect - length(which(t.infect > data$ons.t[i] + onset.iso))
    t.infect <- t.infect[which(t.infect <= data$ons.t[i] + onset.iso)]
  }
  if (n.infect>0){
    data$inf.t[(counter+1):(counter+n.infect)] <- t.infect
    data$infector[(counter+1):(counter+n.infect)] <- i
    counter=counter+n.infect
  }
  if (counter==i) break
}

# si: serial intervals
n.si <- 500
data$si <- NA
data$si[(n.inf0+1):(n.inf0+n.si)] <- data$ons.t[(n.inf0+1):(n.inf0+n.si)] - data$ons.t[data$infector[(n.inf0+1):(n.inf0+n.si)]] 

print(mean(data$si, na.rm=T))
}
# function end


R0.range <- c(3.0,2.5,2,1.5,1,0.5,0.3)
iso.delay.range <- 0:10

# Simulation for all scenarios
N.sim = 200

out.sim.si <- array(NA, dim=c(length(iso.delay.range),length(R0.range),N.sim))
for (k in 1:N.sim){
  for (j in 1:length(R0.range)){
    for (i in 1:length(iso.delay.range)){
      out.sim.si[i,j,k] <- sim.si(R0=R0.range[j], onset.iso=iso.delay.range[i])
    }
  }
}

mean.out.sim.si <- data.frame(matrix(NA, nrow=length(iso.delay.range), ncol=length(R0.range)+1))
LB.out.sim.si <- data.frame(matrix(NA, nrow=length(iso.delay.range), ncol=length(R0.range)+1))
UB.out.sim.si <- data.frame(matrix(NA, nrow=length(iso.delay.range), ncol=length(R0.range)+1))

colnames(mean.out.sim.si) <- colnames(LB.out.sim.si) <- colnames(UB.out.sim.si) <- c('delay',paste0('R0_',R0.range))
mean.out.sim.si$delay <- LB.out.sim.si$delay <- UB.out.sim.si$delay <- iso.delay.range

for (j in 1:length(R0.range)){
  for (i in 1:length(iso.delay.range)){
    mean.out.sim.si[i,j+1] <- mean(out.sim.si[i,j,])
    LB.out.sim.si[i,j+1] <- quantile(out.sim.si[i,j,],0.025)
    UB.out.sim.si[i,j+1] <- quantile(out.sim.si[i,j,],0.975)
  }
}

# mean, 95% CI for the simulation scenarios
mean.out.sim.si
LB.out.sim.si
UB.out.sim.si



