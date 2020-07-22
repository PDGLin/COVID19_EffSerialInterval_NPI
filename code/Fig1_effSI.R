#########################################
## Fig.1 non-overlapping periods
## Ver 1 (2020/06/30)
## by Dr. Lin Wang & Dr. Sheikh Taslim Ali    
## Contact: lw660@cam.ac.uk              
## Website: https://www.pdg.gen.cam.ac.uk/ 
###############################################

rm(list=ls())

setwd(" ") # Must Do: Set path of main folder !! 
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

Covid19.China <- Read_TransPairs( file_name );


####  Fig. 1A: Histogram of empirical serial intervals (based on transmission pairs)  ####
## pre-peak:   9 - 22 (January 9 - 22)
## peak-week: 23 - 29 (January 23 - 29)
## post-peak: 30 - 44 (January 30 - February 13)

Plot_HistSI_nonOverlap <- function( panelTag_vec, onset_LB_UB_Mat, ymax_raw_vec, contacts ) {
  win_LBvec = onset_LB_UB_Mat[, 1];
  win_UBvec = onset_LB_UB_Mat[, 2];
  nwin = length(win_UBvec);
  
  df_hist_SI_full = tibble( idx_win=integer(), SI=integer() );
  df_sample_full  = tibble( LB=integer(), UB=integer(), size=integer() );
  
  df_SI_stat_full = tibble(
    idx_win = integer(),  
    y = double(),
    medianSI = double(), 
    SI_LB = double(), 
    SI_UB = double(),
    SI_mean = double(),
    SI_sd = double()
    );
  
  for (idx_window in 1:nwin) {
    panelTag = panelTag_vec[idx_window];
    LB_win = win_LBvec[idx_window]; 
    UB_win = win_UBvec[idx_window]; 
    
    sub_contacts <- contacts %>% filter( date_onset_infector >= LB_win & date_onset_infector <= UB_win );

    serialIntervalData <- sub_contacts %>% pull("serialInterval");
    print( c(LB_win, UB_win, length(serialIntervalData) ) );
    
    df_sample = tibble( LB=LB_win, UB=UB_win, size=length(serialIntervalData) );
    df_sample_full <- bind_rows(df_sample_full, df_sample);
    
    df_hist_SI = tibble( idx_win = rep(panelTag, length(serialIntervalData)), SI = serialIntervalData );
    
    ymax_vec = switch (
      idx_window,
      c(0, ymax_raw_vec[1]),
      c(0, ymax_raw_vec[2]),
      c(0, ymax_raw_vec[3]),
      c(0, ymax_raw_vec[4])
      )
    df_SI_stat = tibble(
      idx_win = rep(panelTag, 2),
      y = ymax_vec,
      medianSI = rep(quantile(serialIntervalData, 0.5), 2),
      SI_LB = rep(quantile(serialIntervalData, 0.25), 2),
      SI_UB = rep(quantile(serialIntervalData, 0.75), 2),
      SI_mean = rep(mean(serialIntervalData), 2),
      SI_sd = rep(sd(serialIntervalData), 2)
      )
    
    df_hist_SI_full <- bind_rows( df_hist_SI_full, df_hist_SI );
    df_SI_stat_full <- bind_rows( df_SI_stat_full, df_SI_stat );
  }
  
  write_csv( df_sample_full,  paste(getwd(), paste0("Fig1A_sampSize_nonOverlap.csv"), sep=.Platform$file.sep) );
  write_csv( df_SI_stat_full, paste(getwd(), paste0("Fig1A_histStat_nonOverlap.csv"), sep = .Platform$file.sep) );
  
  # plotting -------------
  linesize = 0.5;
  my_breaks <- seq(0, ymax_vec[2], 0.05)
  
  myBlack = "black";
  myRed = "#d7191c";
  myDarkRed = "#800026";
  myBlue = "#67a9cf";
  myDarkBlue = "#2b83ba";
  
  ps1 <- ggplot(df_hist_SI_full) +
    facet_grid(idx_win ~ ., scales = "free_y") +
    geom_histogram(aes(x=SI, fill=factor(idx_win), y = ..density..), colour="black", binwidth=1, size=linesize) +
    geom_line(aes(x=medianSI, y=y), df_SI_stat_full, color=myDarkRed,  size=linesize, linetype = "dashed") +
    geom_line(aes(x=SI_LB, y=y),    df_SI_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    geom_line(aes(x=SI_UB, y=y),    df_SI_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    scale_y_continuous(expand=c(0,0), breaks=my_breaks) +
    scale_x_continuous(breaks=seq(-10, 25, 5)) +
    scale_color_manual(values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    scale_fill_manual(
      breaks=panelTag_vec, 
      values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    labs(x ="Empirical serial intervals (days)") +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      # axis.text.x = element_blank(),
      # axis.text.y = element_blank(),
      # axis.title.y = element_blank(),
      # axis.title.x = element_blank(),
      axis.text = element_text(size=7),
      axis.line=element_line(color="black", size=linesize)
    )
  
  newlist = list(
    plot = ps1,
    df_sampleSize = df_sample_full,
    df_SI_stat_full = df_SI_stat_full
    )
  return( newlist )
}

# Plot SI for 3 non-overlap phase & whole 36-day period
panelTag_vec = seq(1, 4)

ymax_raw_vec = c(0.15, 0.15, 0.20, 0.15)
onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

hist3phase <- Plot_HistSI_nonOverlap(
  panelTag_vec,
  onset_LB_UB_Mat,
  ymax_raw_vec, 
  Covid19.China$contacts 
  )

# X11();
# plot.new();
# print(hist3phase$plot)


####  Fig.1B: MCMC Fittings: Full data, Stratification by isolation delay & household settings  ####

MCMC_FitSI_House_Isol <- function( likelihood_options, df_plot, onset_LB_UB_Mat, contacts ) {
  n_stratification = length(df_plot$is_plot_isol);
  win_LBvec = onset_LB_UB_Mat[, 1];
  win_UBvec = onset_LB_UB_Mat[, 2];
  nwin = length(win_UBvec);
  phaseTag_vec = rev( seq(1, nwin) )
  
  numStepsPerParameter = 100000;  # To Do: Choose Chain Length of MCMC
  
  MCMC_pars_ls <- SetPars_MCMC( numStepsPerParameter );
  df_pars_estim_full = MCMC_pars_ls$df_pars_estim_full;
  df_pars_estim = MCMC_pars_ls$df_pars_estim;
  
  for (idx_window in 1:nwin) {
    LB_win = win_LBvec[idx_window]; 
    UB_win = win_UBvec[idx_window];
    for (idx_stra in 1:n_stratification) {
      df_pars_estim$phase <- phaseTag_vec[idx_window]
      df_pars_estim$onsetDay_LB <- LB_win;
      df_pars_estim$onsetDay_UB <- UB_win;
      print(c(LB_win, UB_win));
      
      sub_contacts <- contacts %>% filter( date_onset_infector >= LB_win & date_onset_infector <= UB_win );
    
      is_plot_house = df_plot$is_plot_house[idx_stra]; 
      is_plot_isol  = df_plot$is_plot_isol[idx_stra];
      is_byHouse    = df_plot$is_byHouse[idx_stra];
      is_byShortIsol = df_plot$is_byShortIsol[idx_stra];
      
      if (is_plot_house==F & is_plot_isol==F) { # plot all cases
        sub_contacts_p2 = sub_contacts;
        
      } else if (is_plot_house==T & is_plot_isol==F) { # plot household or non-household
        if (is_byHouse) { 
          sub_contacts_p2 <- sub_contacts %>% filter( str_detect(sub_contacts$is_household, "yes") ); # household pairs
        } else { 
          sub_contacts_p2 <- sub_contacts %>% filter( str_detect(sub_contacts$is_household, "no") );  # non-household pairs
        }
      } else if (is_plot_house==F & is_plot_isol==T) { # plot shorter or longer isolation
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
      df_pars_estim <- MCMC_posterior_analyser(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data)
      
      df_pars_estim_full = bind_rows( df_pars_estim_full, df_pars_estim );
    }
  }
  
  df_pars_estim_full$stratification <- rep(seq(length(df_plot$is_plot_house), 1), nwin)
  
  return(df_pars_estim_full)
}


## Fitting ----------------
df_plot = tibble(
  is_plot_house = c(T, T, F, F, F),
  is_plot_isol = c(F, F, T, T, F),
  is_byHouse = c(T, F, F, F, F),
  is_byShortIsol = c(F, F, T, F, F) 
  )

onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

df_pars_estim_full <- MCMC_FitSI_House_Isol( likelihood_options, df_plot, onset_LB_UB_Mat, Covid19.China$contacts )

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}
save(
  df_pars_estim_full,
  file = paste(getwd(), paste0("Fig1B_effSI_", file_nm_write, "_nonOverlap.Rdata"), sep = .Platform$file.sep)
  )
write_csv(
  df_pars_estim_full,
  paste(getwd(), paste0("Fig1B_effSI_", file_nm_write, "_nonOverlap.csv"), sep = .Platform$file.sep )
)

## plotting ----------------------
myBlack = "black";
myRed = "#d7191c";
myDarkRed = "#800026";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

linesize = 0.5;
pd = position_dodge(0.5);

ps1 <- ggplot(
  df_pars_estim_full, 
  aes(x=fitMedian, xmin=fitIQR_LB, xmax=fitIQR_UB, y=phase, group=stratification, color=factor(stratification))
  ) + 
  geom_point(position=pd, size=2.5) +
  geom_errorbarh(position=pd, height=0, size=linesize+0.25) +
  scale_x_continuous(expand=c(0,0), breaks=seq(0, 12, 2)) +
  scale_color_manual(
    breaks = unique(df_pars_estim_full$stratification),
    values = c("#0072B2", "#E69F00", "#009E73", "#D55E00", "grey40")
  ) +
  scale_fill_manual(
    breaks = unique(df_pars_estim_full$stratification),
    values = c("#0072B2", "#E69F00", "#009E73", "#D55E00", "grey40")
  ) +
  coord_cartesian( xlim=c(-1.5, 13.5) ) +
  labs( x="Fitted serial intervals (days)" ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    # axis.text.x = element_blank(),
    # axis.text.y = element_blank(),
    # axis.title.y = element_blank(),
    # axis.title.x = element_blank(),
    axis.text = element_text(size=7),
    axis.line = element_line(color="black", size=linesize)
    )

plot.Fig1 <- plot_grid(
  hist3phase$plot,
  ps1,
  nrow = 1,
  ncol = 2
  )

X11()
plot.new()
print(plot.Fig1)

ggsave2(
  paste0("Fig1_effSI_", toString(file_nm_write), "_nonoverlap.eps"), 
  plot=plot.Fig1, 
  width=18, height=12, units="cm", device="eps", dpi=300
  )




