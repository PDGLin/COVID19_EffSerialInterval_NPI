##############################################################
## Fig.S5, Serial interval distribution                     ##
## estimated by probabilistic model of a transmission pair  ##
## Wrote by Dr. Lin Wang, 2020/06/30                        ##                                                                 
## Email: lw660@cam.ac.uk                                   ##
## Web: https://www.pdg.gen.cam.ac.uk/                      ##
##############################################################


rm(list=ls())

setwd("")  # Must Do: Specify the working directory !!
getwd()

####  Install packages if necessary  ####
installed_package_details <- installed.packages()
installed_package_names    <- as.vector(installed_package_details[, 1])
installed_package_versions <- as.vector(installed_package_details[, 3])

## List of package dependencies ----
dependencies <- c(
  "ggplot2",
  "readr",
  "dplyr",
  "pracma"
)
uninstalled_dependencies <- dependencies[!dependencies %in% installed_package_names]
lapply(uninstalled_dependencies, install.packages)

## load packages -----------------
library(ggplot2)
library(readr)
library(dplyr)
library(pracma)

source("myFun_probModel.R");

linesize = 0.5;
pointsize = 1.5

####  Compute the probabilistic model of a transmission pair  ####
## Infectiousness: Probability to infect the infectee \tau days after infector i starts to be infectious at time T^{EA}_i 
## parameters estimated from (He et al Nat Med 2020, https://www.nature.com/articles/s41591-020-0869-5)
## C_i   -> time duration from the start of infectiousness to symptom onset
## K     -> shape of gamma distribution for infectiousness
## theta -> scale of gamma distribution for infectiousness

pars_estimated <- tibble(
  C_i = seq(1, 7, 1),
  shape = c(4.070587, 2.146901, 2.542802, 2.735434, 1.9622, 1.528505, 1.248171),
  scale = c(0.463752, 1.33708, 1.496252, 1.711471, 2.907554, 4.380512, 6.155327)
  )

D_i = seq(0, 10, by=1)  # Isolation delay

is_isolateInfector_beforeInfecteeOnset = F;  # To Do: whether infector is assumed to be isolated before infectee's onset

df_stat = tibble(
  Di = double(),
  SI_mean = double(),
  SI_median = double(),
  LB = double(),
  UB = double(),
  C_i = double()
  )

for (idx_Ci in c(2, 4, 6)) {
  A = myFun_probModel(
    pars_estimated$C_i[idx_Ci],
    pars_estimated$shape[idx_Ci],
    pars_estimated$scale[idx_Ci],
    is_isolateInfector_beforeInfecteeOnset,
    D_i
    )
  A <- A %>% mutate( C_i = rep(pars_estimated$C_i[idx_Ci], length(D_i)) );
  
  df_stat = bind_rows( df_stat, A )
}


pd = position_dodge(0.6);

df_stat_p1 <- df_stat # %>% filter(C_i == 2 | C_i == 4 | C_i == 6)

plot_figS5 <- ggplot(df_stat_p1, aes(x=Di, y=SI_mean, ymin=LB, ymax=UB, group=C_i, color=factor(C_i))) +
  geom_point(position=pd, size=pointsize) +
  geom_errorbar(position=pd, width=0, size=linesize) +
  scale_x_reverse(breaks=rev(unique(df_stat_p1$Di))) +
  scale_color_manual( values=c("#009E73", "#E69F00", "#D55E00") ) +
  coord_cartesian( xlim=c(10, 0), ylim=c(0, 8) ) +
  labs(x="Isolation delay (days)", y="Serial interval distribution (probabilistic model)") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border= element_blank(),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    axis.line.x = element_line(color="black", size=linesize),
    axis.line.y = element_line(color="black", size=linesize),
    axis.ticks.y =element_line(size=linesize, color="black"),
    axis.ticks.x =element_line(size=linesize, color="black"),
    # axis.title =  element_blank()
  )

X11()
plot.new()
print(plot_figS5)

