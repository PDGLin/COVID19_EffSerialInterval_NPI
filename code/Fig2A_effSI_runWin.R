
rm(list=ls())

setwd("") # Must Do: Specify the working directory !!
getwd()

source("myFun_readTransPairs.R");
source("myFun_choose_likelihood.R")
source("mcmcProposal.R");
source("MCMC.R");
source("Get_DIC.R");
source("logLikelihood.R");
source("logPrior.R");
source("MCMC_posterior_analyser.R")
source("SetPars_MCMC.R")

####  Install packages if necessary  ####

installed_package_details <- installed.packages()
installed_package_names    <- as.vector(installed_package_details[, 1])
installed_package_versions <- as.vector(installed_package_details[, 3])

# List of package dependencies
dependencies <- c(
  "readr",
  "dplyr",
  "stringr",
  "ggplot2",
  "cowplot",
  "pracma",
  "ggpubr",
  "matrixStats",
  "extrafont"
  )
uninstalled_dependencies <- dependencies[!dependencies %in% installed_package_names]
lapply(uninstalled_dependencies, install.packages)

## load packages ----
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(pracma)
library(ggpubr)
library(matrixStats)
library(extrafont)
loadfonts(device="win")

myBlack = "black";
myRed = "#d7191c";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

####  To Do: Choose likelihood func  ####
# To Do: Choose Likelihood Model: "Normal_fit", "Gumbel_fit", "Logist_fit", "LogNormal_fit", or "transPair_fit"
tar_likelihood_model = "Normal_fit"; 

likelihood_options = myFun_choose_likelihood( tar_likelihood_model );

####  read "linelist" and "transmission pair" data  ####

folder_nm = "InputData";
file_name = paste(getwd(), folder_nm, "TableS1_1407TransPairs.csv", sep = .Platform$file.sep);

Covid19.China <- Read_TransPairs( file_name )  # To Do: Check if first index case of each cluster needs to be excluded


####  MCMC fittings: by Isolation  ####

MCMC_FitSI_Isol <- function( likelihood_options, df_plot, onset_LB_UB_Mat, contacts ) {
  n_stratification = length(df_plot$is_plot_isol);
  win_LBvec = onset_LB_UB_Mat[, 1];
  win_UBvec = onset_LB_UB_Mat[, 2];
  nwin = length(win_UBvec);
  len_win = length(seq(win_UBvec[1], win_LBvec[1]));
  
  ## MCMC Settings 
  numStepsPerParameter = 100000;  # To Do: Choose Chain Length
  
  MCMC_pars_ls <- SetPars_MCMC( numStepsPerParameter );
  df_pars_estim_full = MCMC_pars_ls$df_pars_estim_full;
  df_pars_estim = MCMC_pars_ls$df_pars_estim;
  
  for (idx_window in 1:nwin) {
    LB_win = win_LBvec[idx_window]; 
    UB_win = win_UBvec[idx_window];
    for (idx_stra in 1:n_stratification) {
      df_pars_estim$onsetDay_LB <- LB_win;
      df_pars_estim$onsetDay_UB <- UB_win;
      print(c(LB_win, UB_win));
      
      sub_contacts <- contacts %>% filter( date_onset_infector >= LB_win & date_onset_infector <= UB_win );
      
      is_plot_isol  = df_plot$is_plot_isol[idx_stra];
      is_byShortIsol = df_plot$is_byShortIsol[idx_stra];
      
      if (!is_plot_isol) { 
        sub_contacts_p2 = sub_contacts; # plot all cases
      } else {
        tar_lag_isol_onset = sub_contacts$delayIsol_infector;
        split_lag = median(tar_lag_isol_onset[ tar_lag_isol_onset < 500]);  # To Do: Choose criterion to split the lag
        
        if (is_byShortIsol) {
          sub_contacts_p2 <- sub_contacts %>% filter( delayIsol_infector < split_lag );
        } else {
          sub_contacts_p2 <- sub_contacts %>% filter( delayIsol_infector >= split_lag & delayIsol_infector < 500 );
        }
      }
      
      data <- sub_contacts_p2$serialInterval;
      df_pars_estim$sample_size <- length(data);

      result = MCMC(data, numStepsPerParameter, likelihood_options);
      
      chainRecord = result$chainRecord;
      df_pars_estim <- MCMC_posterior_analyser(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data);
      df_pars_estim_full = bind_rows( df_pars_estim_full, df_pars_estim );
    }
  }
  
  df_pars_estim_full$stratification <- rep(seq(1, length(df_plot$is_plot_isol)), nwin)
  
  return(df_pars_estim_full)
}


## Fitting ------------------

df_plot = tibble(
  is_plot_isol   = c(F, T, T),
  is_byShortIsol = c(F, T, F)
  )

len_window = 13;  # To Do: Choose time window

onset_LB_vec = seq( 
  min(Covid19.China$contacts$date_onset_infector), 
  max(Covid19.China$contacts$date_onset_infector) - len_window, 
  by=1
  ) 
onset_UB_vec = onset_LB_vec + len_window
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

df_pars_estim_full <- MCMC_FitSI_Isol( likelihood_options, df_plot, onset_LB_UB_Mat, Covid19.China$contacts )

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}
file_nm_write = paste0("Fig2A_effSI_", file_nm_write, "_", toString(len_window+1), "day", "RunWin")
save(
  df_pars_estim_full, 
  file = paste(getwd(), paste0(file_nm_write, ".Rdata"), sep = .Platform$file.sep)
  )
write_csv(
  df_pars_estim_full,
  paste( getwd(), paste0(file_nm_write, ".csv"), sep = .Platform$file.sep )
  )

## Plotting ----

linesize = 0.5;

myBlack = "black";
myRed = "#d7191c";
myDarkRed = "#800026";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

pd = position_dodge(0.6);

plot.Fig2A <- ggplot(
  df_pars_estim_full, 
  aes(x=onsetDay_LB, y=fitMedian, ymin=fitIQR_LB, ymax=fitIQR_UB, group=stratification, color=factor(stratification))
  ) + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(position=pd, width=0, size=linesize+0.25) +
  scale_y_continuous(breaks=seq(0, 12, 3) ) +
  scale_x_continuous(
    breaks= unique(df_pars_estim_full$onsetDay_LB),
    labels=c(
      "9" = "Jan 9 - Jan 22",
      "10" = "Jan 10 - Jan 23",
      "11" = "Jan 11 - Jan 24",
      "12" = "Jan 12 - Jan 25",
      "13" = "Jan 13 - Jan 26",
      "14" = "Jan 14 - Jan 27",
      "15" = "Jan 15 - Jan 28",
      "16" = "Jan 16 - Jan 29",
      "17" = "Jan 17 - Jan 30",
      "18" = "Jan 18 - Jan 31",
      "19" = "Jan 19 - Feb 01",
      "20" = "Jan 20 - Feb 02",
      "21" = "Jan 21 - Feb 03",
      "22" = "Jan 22 - Feb 04",
      "23" = "Jan 23 - Feb 05",
      "24" = "Jan 24 - Feb 06",
      "25" = "Jan 25 - Feb 07",
      "26" = "Jan 26 - Feb 08",
      "27" = "Jan 27 - Feb 09",
      "28" = "Jan 28 - Feb 10",
      "29" = "Jan 29 - Feb 11",
      "30" = "Jan 30 - Feb 12",
      "31" = "Jan 31 - Feb 13"
      )
    ) +
  scale_color_manual(
    values=c("grey40", "#009E73", "#D55E00")
  ) +
  coord_cartesian(
    xlim=c(9.5, 30.5),
    ylim=c(-1, 12.5)
  ) +
  labs(x="14-day running windows", y="Fitted serial intervals (days)") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle=45, hjust=1),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank(),
    # axis.title.x = element_blank(),
    axis.text = element_text(size=12),
    axis.line = element_line(color="black", size=linesize)
  )


X11()
plot.new()
print(plot.Fig2A)

ggsave(
  paste0(file_nm_write, ".eps"), 
  plot=plot.Fig2A, 
  width=20, height=8, units="cm", device="eps", dpi=300
  )





  










# scale_x_continuous(
#   breaks= unique(df_pars_estim_full$onsetDay_LB),
#   labels=c(
#     "7" = "Jan 7 - Jan 16",
#     "8" = "Jan 8 - Jan 17",
#     "9" = "Jan 9 - Jan 18",
#     "10" = "Jan 10 - Jan 19",
#     "11" = "Jan 11 - Jan 20",
#     "12" = "Jan 12 - Jan 21",
#     "13" = "Jan 13 - Jan 22",
#     "14" = "Jan 14 - Jan 23",
#     "15" = "Jan 15 - Jan 24",
#     "16" = "Jan 16 - Jan 25",
#     "17" = "Jan 17 - Jan 26",
#     "18" = "Jan 18 - Jan 27",
#     "19" = "Jan 19 - Jan 28",
#     "20" = "Jan 20 - Jan 29",
#     "21" = "Jan 21 - Jan 30",
#     "22" = "Jan 22 - Jan 31",
#     "23" = "Jan 23 - Feb 01",
#     "24" = "Jan 24 - Feb 02",
#     "25" = "Jan 25 - Feb 03",
#     "26" = "Jan 26 - Feb 04",
#     "27" = "Jan 27 - Feb 05"
#     )
#   ) +
# scale_x_continuous(
#   breaks= unique(df_pars_estim_full$onsetDay_LB),
#   labels=c(
#     "7" = "Jan 7 - Jan 16",
#     "8" = "Jan 8 - Jan 17",
#     "9" = "Jan 9 - Jan 18",
#     "10" = "Jan 10 - Jan 19",
#     "11" = "Jan 11 - Jan 20",
#     "12" = "Jan 12 - Jan 21",
#     "13" = "Jan 13 - Jan 22",
#     "14" = "Jan 14 - Jan 23",
#     "15" = "Jan 15 - Jan 24",
#     "16" = "Jan 16 - Jan 25",
#     "17" = "Jan 17 - Jan 26",
#     "18" = "Jan 18 - Jan 27",
#     "19" = "Jan 19 - Jan 28",
#     "20" = "Jan 20 - Jan 29",
#     "21" = "Jan 21 - Jan 30",
#     "22" = "Jan 22 - Jan 31",
#     "23" = "Jan 23 - Feb 01",
#     "24" = "Jan 24 - Feb 02",
#     "25" = "Jan 25 - Feb 03",
#     "26" = "Jan 26 - Feb 04",
#     "27" = "Jan 27 - Feb 05",
#     "28" = "Jan 28 - Feb 06",
#     "29" = "Jan 29 - Feb 07",
#     "30" = "Jan 30 - Feb 08",
#     "31" = "Jan 31 - Feb 09",
#     "32" = "Jan 31 - Feb 10",
#     "33" = "Jan 31 - Feb 11",
#     "34" = "Jan 31 - Feb 12",
#     "35" = "Jan 31 - Feb 13"
#     )
#   ) +
