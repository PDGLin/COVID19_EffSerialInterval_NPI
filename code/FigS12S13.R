
rm(list=ls())

setwd("C:/Users/Lin/Dropbox/COVIE19/Collaborations/Shared HKU_IP_Camb/serial_interval_covid19/Manuscript/Science/Submission Science Report/Revision V1/Data Code to Upload") # Must Do: Specify the working directory !!
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
## "Normal_fit", "Gumbel_fit", "Logist_fit", "LogNormal_fit", or "transPair_fit"
tar_likelihood_model = "Normal_fit"; 

likelihood_options = myFun_choose_likelihood( tar_likelihood_model );

####  read "linelist" and "transmission pair" data  ####

folder_nm = "InputData";
file_name = paste(getwd(), folder_nm, "TableS1_1407TransPairs.csv", sep = .Platform$file.sep);

Covid19.China <- Read_TransPairs( file_name ) 

####  Histogram of empirical serial intervals for each non-overlapping period  ####
## pre-peak:   9 - 22 (January 9 - 22)
## peak-week: 23 - 29 (January 23 - 29)
## post-peak: 30 - 44 (January 30 - February 13)

Plot_HistSI_nonOverlap <- function( panelTag_vec, onset_LB_UB_Mat, ymax_raw_vec, contacts, len_win_samp ) {
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
    
    serialIntervalData <- sub_contacts$date_onset_infectee - sub_contacts$date_onset_infector;
    print( c(LB_win, UB_win, length(serialIntervalData) ) );
    
    df_sample = tibble( LB=LB_win, UB=UB_win, size=length(serialIntervalData) )
    df_sample_full <- bind_rows(df_sample_full, df_sample)
    df_hist_SI = tibble(
      idx_win = rep(panelTag, length(serialIntervalData)),
      SI = serialIntervalData
      )
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
    df_hist_SI_full <- bind_rows( df_hist_SI_full, df_hist_SI )
    df_SI_stat_full <- bind_rows( df_SI_stat_full, df_SI_stat )
  }
  
  write_csv( 
    df_sample_full, 
    paste(getwd(), paste0("sampSize_nonOverlap_lenWin", toString(len_win_samp), "day.csv"), sep=.Platform$file.sep) 
    )
  write_csv(
    df_SI_stat_full,
    paste(getwd(), paste0("histStat_nonOverlap_lenWin", toString(len_win_samp), "day.csv"), sep = .Platform$file.sep)
  )
  
  myBlack = "black";
  myRed = "#d7191c";
  myDarkRed = "#800026";
  myBlue = "#67a9cf";
  myDarkBlue = "#2b83ba";
  
  linesize = 0.5;
  my_breaks <- seq(0, ymax_vec[2], 0.05)
  
  ps1 <- ggplot(df_hist_SI_full) +
    facet_grid(idx_win ~ ., scales = "free_y") +
    geom_histogram(aes(x=SI, fill=factor(idx_win), y = ..density..), colour="black", binwidth=1, size=linesize) + # fill="white",
    geom_line(aes(x=medianSI, y=y), df_SI_stat_full, color=myDarkRed,  size=linesize, linetype = "dashed") +
    geom_line(aes(x=SI_LB, y=y),    df_SI_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    geom_line(aes(x=SI_UB, y=y),    df_SI_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    scale_y_continuous(expand=c(0,0), breaks=my_breaks) +
    scale_x_continuous(breaks=seq(-10, 25, 5)) +
    scale_color_manual(values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    scale_fill_manual(
      breaks=panelTag_vec, 
      values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    labs(
      y="Relative frequency", 
      x=paste0("Empirical serial interval (days) with \n", toString(len_win_samp), "-day uncertainty in recall bias" )
      ) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      axis.text.x = element_text(size=12),
      axis.text.y = element_text(size=12),
    #  axis.title.y = element_blank(),
    #  axis.title.x = element_blank(),
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

## set non-overlapping periods -----------------------
panelTag_vec = seq(1, 4)

ymax_raw_vec = c(0.15, 0.15, 0.20, 0.15) # c(0.10, 0.13, 0.20, 0.11)
onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

## plot transmission pair data without resampling of recall bias --------------

len_win_samp = 0;

hist3phase <- Plot_HistSI_nonOverlap(
  panelTag_vec,
  onset_LB_UB_Mat,
  ymax_raw_vec, 
  Covid19.China$contacts,
  len_win_samp
  )

####  Add uncertainty in recall bias for each case  ####

## function to add uncertainty -------------------------
Set_onset_in_pair <-function( id_vec_search, raw_linelist ) {
  tar_onset_vec = sapply(
    id_vec_search, 
    function (tar_id, case_linelist) case_linelist %>% filter( str_detect(case_linelist$id, tar_id) ) %>% pull("date_onset"),
    "case_linelist" = raw_linelist
    )
  return(tar_onset_vec);
}

Sample_recall_bias <- function(raw_data, len_win) {
  raw_linelist = raw_data$linelist;
  raw_contacts = raw_data$contacts;
  
  raw_linelist$date_onset = sapply(
    raw_linelist$date_onset, 
    function (tar_date, len_win) sample( seq(tar_date - len_win, tar_date + len_win, by=1), 1, replace=F ),
    "len_win" = len_win
    )
  raw_contacts$date_onset_infector = unname( Set_onset_in_pair( raw_contacts$from, raw_linelist ) );
  raw_contacts$date_onset_infectee = unname( Set_onset_in_pair( raw_contacts$to,   raw_linelist ) );
  raw_contacts$serialInterval = raw_contacts$date_onset_infectee - raw_contacts$date_onset_infector;
  
  contacts_resampled_onset = raw_contacts;
  return(contacts_resampled_onset)
}


## Sampling 2-day uncertainty in recall bias --------------
len_win_samp = 2;

contacts_resampled_onset <- Sample_recall_bias( Covid19.China, len_win_samp );

hist3phase_samp2day <- Plot_HistSI_nonOverlap(
  panelTag_vec,
  onset_LB_UB_Mat,
  ymax_raw_vec, 
  contacts_resampled_onset,
  len_win_samp
  )

## Sampling 5-day uncertainty in recall bias --------------
len_win_samp = 4;

contacts_resampled_onset <- Sample_recall_bias( Covid19.China, len_win_samp );

hist3phase_samp4day <- Plot_HistSI_nonOverlap(
  panelTag_vec,
  onset_LB_UB_Mat,
  ymax_raw_vec, 
  contacts_resampled_onset,
  len_win_samp
  )

## Plot: Comparing 0-day, 2-day, 5-day recall bias -----------
pHist_comp_uncertainty <- plot_grid(
  hist3phase$plot,
  hist3phase_samp2day$plot,
  hist3phase_samp4day$plot,
  nrow = 1,
  ncol = 3
  )

X11();
plot.new();
print(pHist_comp_uncertainty)

ggsave2(
  "FigS12_hist_nonOverlap_recallBias.eps", 
  plot = pHist_comp_uncertainty, 
  width=18, height=15, units="cm", device="eps", dpi=300
)



####  MCMC fitting: Isolation & Household  ####

MCMC_FitSI_House_Isol <- function( likelihood_options, onset_LB_UB_Mat, Covid19.data, len_win_vec ) {
  win_LBvec = onset_LB_UB_Mat[, 1];
  win_UBvec = onset_LB_UB_Mat[, 2];
  nwin = length(win_UBvec);
  phaseTag_vec = rev( seq(1, nwin) );
  
  numStepsPerParameter = 100000;  # To Do: Choose Chain Length of MCMC
  MCMC_pars_ls <- SetPars_MCMC( numStepsPerParameter );
  df_pars_estim_full = MCMC_pars_ls$df_pars_estim_full;
  df_pars_estim = MCMC_pars_ls$df_pars_estim;
  
  num_len_win = length(len_win_vec);
  for (idx_len_win in 1:num_len_win) {
    contacts = switch (
      toString(idx_len_win),
      "1" = Covid19.data$contacts,
      "2" = Sample_recall_bias( Covid19.data, len_win_vec[idx_len_win] ),
      "3" = Sample_recall_bias( Covid19.data, len_win_vec[idx_len_win] ),
      "4" = Sample_recall_bias( Covid19.data, len_win_vec[idx_len_win] )
      )
    for (idx_window in 1:nwin) {
      LB_win = win_LBvec[idx_window]; 
      UB_win = win_UBvec[idx_window];
    
      df_pars_estim$phase <- phaseTag_vec[idx_window]
      df_pars_estim$onsetDay_LB <- LB_win;
      df_pars_estim$onsetDay_UB <- UB_win;
      print(c(LB_win, UB_win));
      
      sub_contacts <- contacts %>% filter( date_onset_infector >= LB_win & date_onset_infector <= UB_win );
      
      data <- sub_contacts$serialInterval;
      df_pars_estim$sample_size <- length(data);
      
      result = MCMC(data, numStepsPerParameter, likelihood_options);
      chainRecord = result$chainRecord;
      
      df_pars_estim <- MCMC_posterior_analyser(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data);
      df_pars_estim_full = bind_rows( df_pars_estim_full, df_pars_estim );
    }
  }
  
  df_pars_estim_full$stratification <- rep( seq(1, num_len_win), each=nwin );
  
  return(df_pars_estim_full)
}


## Household & Short IsolDelay: Non-overlap Windows ------------------
len_win_vec = c(0, 1, 2, 4)

onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

##  fitting by 0-day, 2-day, 4-day uncertainty  --------------------

df_pars_estim_full <- MCMC_FitSI_House_Isol( likelihood_options, onset_LB_UB_Mat, Covid19.China, len_win_vec )

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}

save( 
  df_pars_estim_full,
  file = paste(
    getwd(), 
    paste0("mcmcStat_", file_nm_write, "_forwardSI_nonOverlap_uncertainty.Rdata"), 
    sep = .Platform$file.sep
    )
  )

file_nm_write = paste0("mcmcStat_", file_nm_write, "_forwardSI_nonOverlap_uncertainty.csv")
write_csv( 
  df_pars_estim_full, 
  paste(getwd(), file_nm_write, sep = .Platform$file.sep) 
  )

## Plotting -----------------------------
linesize = 0.5;

myBlack = "black";
myRed = "#d7191c";
myDarkRed = "#800026";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

pd = position_dodge(0.5);

ps1 <- ggplot(
  df_pars_estim_full, 
  aes(y=fitMedian, ymin=fitIQR_LB, ymax=fitIQR_UB, x=phase, group=stratification, color=factor(stratification))
  ) +
  geom_point(position=pd, size=2.5) +
  geom_errorbar(position=pd, width=0, size=linesize+0.25) +
  scale_y_continuous(expand=c(0,0), breaks=seq(-2, 12, 2)) +
  scale_x_continuous(
    breaks= sort(unique(df_pars_estim_full$phase)),
    labels=c(
      "1" = "Jan 9 - Feb 13 \n (whole period)",
      "2" = "Jan 30 - Feb 13 \n (post-peak)",
      "3" = "Jan 23 - Jan 29 \n (peak-week)",
      "4" = "Jan 9 - Jan 22 \n (pre-peak)"
      )
    ) +
  scale_color_manual(
    breaks = unique(df_pars_estim_full$stratification),
    values = c("grey40", "#D55E00", "#009E73", "#E69F00"),
    labels = c("Original data", "1-day uncertainty", "2-day uncertainty", "4-day uncertainty")
  ) +
  labs(y="Fitted serial intervals (days)") +
  coord_flip( ylim=c(-3, 13.5)  ) +
  theme(
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=12),
    axis.text = element_text(size=12),
    axis.line = element_line(color="black", size=linesize)
    )


X11()
plot.new()
print(ps1)

ggsave2(
  paste0("FigS13_NormalFitSI_nonoverlap_uncertainty.eps"), 
  plot=ps1, 
  width=10, height=12, units="cm", device="eps", dpi=300
)









































## MCMC settings: Age & Sex -------------------------------------------

# Plot_MCMC_Age_Sex_fit <- function( likelihood_options, df_plot, onset_LB_UB_Mat, tar_linelist ) {
#   n_stratification = length(df_plot$is_byAge);
#   win_LBvec = onset_LB_UB_Mat[, 1];
#   win_UBvec = onset_LB_UB_Mat[, 2];
#   nwin = length(win_UBvec);
#   phaseTag_vec = rev( seq(1, nwin) )
#   
#   ## MCMC Settings 
#   numStepsPerParameter = 100000;  # To Do: Choose Chain Length
#   
#   MCMC_pars_ls <- SetPars_MCMC( numStepsPerParameter )
#   startingPoint = MCMC_pars_ls$startingPoint
#   LB = MCMC_pars_ls$LB
#   UB = MCMC_pars_ls$UB
#   minProbAccept = MCMC_pars_ls$minProbAccept
#   maxProbAccept = MCMC_pars_ls$maxProbAccept
#   df_pars_estim_full = MCMC_pars_ls$df_pars_estim_full
#   df_pars_estim = MCMC_pars_ls$df_pars_estim
#   
#   age_primaryCase_vec = as.numeric( tar_linelist$s_age );
#   age_median = median( age_primaryCase_vec[ !is.na(age_primaryCase_vec) ] );
#   age_primaryCase_vec[ is.na(age_primaryCase_vec) ] <- 10000
#   
#   for (idx_window in 1:nwin) {
#     LB_win = win_LBvec[idx_window]; 
#     UB_win = win_UBvec[idx_window];
#     
#     for (idx_stra in 1:n_stratification) {
#       df_pars_estim$phase <- phaseTag_vec[idx_window]
#       df_pars_estim$onsetDay_LB <- LB_win;
#       df_pars_estim$onsetDay_UB <- UB_win;
#       print(c(LB_win, UB_win))
#       
#       is_byAge = df_plot$is_byAge[idx_stra]; 
#       is_bySex = df_plot$is_bySex[idx_stra];
#       is_young = df_plot$is_young[idx_stra];
#       is_male  = df_plot$is_male[idx_stra];
#       
#       tar_linelist_p1 <- tar_linelist %>% filter(
#         tar_linelist[["s_Date of onset"]] >= LB_win & 
#         tar_linelist[["s_Date of onset"]] <= UB_win
#         )
#       
#       if (is_byAge) { # stratification by age
#         age_primaryCase_vec_p1 = as.numeric( tar_linelist_p1$s_age );
#         age_primaryCase_vec_p1[ is.na(age_primaryCase_vec_p1) ] <- 10000
#         if (is_young) {
#           tar_linelist_p2 <- tar_linelist_p1 %>% filter( age_primaryCase_vec_p1 < age_median );
#         } else {
#           tar_linelist_p2 <- tar_linelist_p1 %>% filter( age_primaryCase_vec_p1 >= age_median & age_primaryCase_vec_p1 < 1000 );
#         }
#         
#       } else if (is_bySex) { # stratification by gender
#         if (is_male) {
#           tar_linelist_p2 <- tar_linelist_p1 %>% filter(
#             !( str_detect(tar_linelist_p1$s_gender, "female") | str_detect(tar_linelist_p1$s_gender, "unknown") )
#           );
#         } else {
#           tar_linelist_p2 <- tar_linelist_p1 %>% filter( str_detect(tar_linelist_p1$s_gender, "female") );
#         }
#       }
#       
#       serialIntervalData <- tar_linelist_p2[["p_Date of onset"]] - tar_linelist_p2[["s_Date of onset"]];
#       df_pars_estim$sample_size <- length(serialIntervalData);
#       
#       data = serialIntervalData;
#       result = MCMC(data, LB, UB, startingPoint, numStepsPerParameter, minProbAccept, maxProbAccept, likelihood_options);
#       
#       chainRecord = result$chainRecord;
#       # probAccept = result$probAccept;
#       df_pars_estim <- MCMC_posterior_analyser(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data)
#       
#       df_pars_estim_full = bind_rows( df_pars_estim_full, df_pars_estim );
#     }
#   }
#   
#   df_pars_estim_full$stratification <- rep(seq(n_stratification, 1), nwin)
#   
#   return(df_pars_estim_full)
# }
# 
# ## run TableS2: Age & Sex
# df_plot = tibble(
#   is_byAge = c(T, T, F, F),
#   is_bySex = c(F, F, T, T),
#   is_young = c(T, F, F, F),
#   is_male  = c(F, F, T, F)
#   )
# 
# onset_LB_vec = c(9, 23, 30, 9)
# onset_UB_vec = c(22, 29, 44, 44)
# onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)
# 
# df_pars_estim_full <- Plot_MCMC_Age_Sex_fit( likelihood_options, df_plot, onset_LB_UB_Mat, tar_linelist )
# 
# if (likelihood_options$is_Normal) {
#   file_nm_write = "NormalFit";
# } else if (likelihood_options$is_Gumbel) {
#   file_nm_write = "GumbelFit";
# } else if (likelihood_options$is_Logistic) {
#   file_nm_write = "LogisFit";
# }
# file_nm_write = paste0("TableS2_", file_nm_write, "_forwardSI_age_sex_nonOverlap")
# 
# save(
#   df_pars_estim_full,
#   file = paste(getwd(), paste0(file_nm_write, ".Rdata"), sep = .Platform$file.sep)
# )
# write_csv(
#   df_pars_estim_full,
#   paste(getwd(), paste0(file_nm_write, ".csv"), sep = .Platform$file.sep )
# )
# 
# 
# 

