###########################################################################################
## Fig.S2, Temporal change of the time delay in isolating infector (isolation delay D_i) 
## Wrote by Dr. Lin Wang, 2020/07/11                                                     
## Email: lw660@cam.ac.uk                                                                
## Web: https://www.pdg.gen.cam.ac.uk/                                                   
###########################################################################################

rm(list=ls())

setwd(" ") # Must Do: Specify the working directory !!
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


#### Fig.S2A: Histogram of empirical serial intervals (based on transmission pairs)  ####
## pre-peak:   9 - 22 (January 9 - 22)
## peak-week: 23 - 29 (January 23 - 29)
## post-peak: 30 - 44 (January 30 - February 13)

Plot_HistIsol_nonOverlap <- function(panelTag_vec, onset_LB_UB_Mat, ymax_raw_vec, contacts) {
  win_LBvec = onset_LB_UB_Mat[, 1];
  win_UBvec = onset_LB_UB_Mat[, 2];
  nwin = length(win_UBvec);
  df_hist_IsolDelay_full = tibble( idx_win=integer(), SI=double() );
  df_sample_full  = tibble( LB=integer(), UB=integer(), size=integer() );
  df_IsolDelay_stat_full = tibble(
    idx_win = integer(),  
    y = double(),
    median = double(), 
    IQR_LB = double(), 
    IQR_UB = double(),
    mean = double(),
    sd = double()
    )
  
  for (idx_window in 1:nwin) {
    panelTag = panelTag_vec[idx_window];
    LB_win = win_LBvec[idx_window]; 
    UB_win = win_UBvec[idx_window]; 
    
    sub_contacts <- contacts %>% filter(date_onset_infector >= LB_win & date_onset_infector <= UB_win);
    
    IsolDelayData <- sub_contacts %>% distinct(from, .keep_all=T) %>% 
      filter(delayIsol_infector < 25 & delayIsol_infector > -10) %>% pull("delayIsol_infector");
    print( c(LB_win, UB_win, length(IsolDelayData) ) );
    
    df_sample = tibble(LB=LB_win, UB=UB_win, size=length(IsolDelayData) );
    df_sample_full <- bind_rows(df_sample_full, df_sample);
    df_hist_IsolDelay = tibble( idx_win=rep(panelTag, length(IsolDelayData)), SI=IsolDelayData );
    
    ymax_vec = switch(idx_window, c(0, ymax_raw_vec[1]), c(0, ymax_raw_vec[2]), c(0, ymax_raw_vec[3]), c(0, ymax_raw_vec[4]));
    
    df_IsolDelay_stat = tibble(
      idx_win = rep(panelTag, 2),
      y = ymax_vec,
      median = rep(quantile(IsolDelayData, 0.5), 2),
      IQR_LB = rep(quantile(IsolDelayData, 0.25), 2),
      IQR_UB = rep(quantile(IsolDelayData, 0.75), 2),
      mean = rep(mean(IsolDelayData), 2),
      sd = rep(sd(IsolDelayData), 2)
      )
    df_hist_IsolDelay_full <- bind_rows( df_hist_IsolDelay_full, df_hist_IsolDelay )
    df_IsolDelay_stat_full <- bind_rows( df_IsolDelay_stat_full, df_IsolDelay_stat )
  }
  
  write_csv(df_sample_full, paste(getwd(), "FigS2A_sampSize_nonOverlap_Di.csv", sep=.Platform$file.sep ) );
  write_csv(df_IsolDelay_stat_full, paste(getwd(), "FigS2A_stat_nonOverlap_Di.csv", sep = .Platform$file.sep) );
  
  # plotting
  myBlack = "black";
  myRed = "#d7191c";
  myDarkRed = "#800026";
  myBlue = "#67a9cf";
  myDarkBlue = "#2b83ba";
  linesize = 0.5;
  
  ps1 <- ggplot(df_hist_IsolDelay_full) +
    facet_grid(idx_win ~ ., scales = "free_y") +
    geom_histogram(aes(x=SI, fill=factor(idx_win), y = ..density..), colour="black", binwidth=1, size=linesize) +
    geom_line(aes(x=median, y=y), df_IsolDelay_stat_full, color=myDarkRed,  size=linesize, linetype = "dashed") +
    geom_line(aes(x=IQR_LB, y=y), df_IsolDelay_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    geom_line(aes(x=IQR_UB, y=y), df_IsolDelay_stat_full, color=myDarkBlue, size=linesize, linetype = "longdash") +
    scale_y_continuous(expand=c(0,0), breaks=seq(0, ymax_vec[2], 0.05)) +
    scale_x_continuous(breaks=seq(-10, 25, 5)) +
    scale_color_manual(values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    scale_fill_manual( breaks=panelTag_vec, values=c("#f03b20", "#feb24c", "#ffeda0", "#FFFFFF")) +
    labs( x="Empirical isolation delay per infector (days)", y="Relative frequency") +
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
      axis.text = element_text(size=10),
      axis.line=element_line(color="black", size=linesize)
    )
  
  newlist = list(
    plot = ps1,
    df_sampleSize = df_sample_full,
    df_IsolDelay_stat_full = df_IsolDelay_stat_full
  )
  return( newlist )
}

# Plot Isolation Delay for 3 non-overlap phase & whole 36-day period
panelTag_vec = seq(1, 4)

ymax_raw_vec = c(0.15, 0.20, 0.15, 0.15)
onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

hist3phase_IsolDelay <- Plot_HistIsol_nonOverlap(
  panelTag_vec,
  onset_LB_UB_Mat,
  ymax_raw_vec, 
  Covid19.China$contacts 
  )

X11();
plot.new();
print(hist3phase_IsolDelay$plot)


####  Fig.S2B: MCMC Fittings: Full data, Stratification by household settings  ####

MCMC_IsolDelay <- function(likelihood_options, df_plot, onset_LB_UB_Mat, contacts ) {
  n_stratification = length(df_plot$is_plot_house);
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
      df_pars_estim$phase <- phaseTag_vec[idx_window];
      df_pars_estim$onsetDay_LB <- LB_win;
      df_pars_estim$onsetDay_UB <- UB_win;
      print(c(LB_win, UB_win));
      
      sub_contacts <- contacts %>% filter(date_onset_infector >= LB_win & date_onset_infector <= UB_win);
      
      is_plot_house = df_plot$is_plot_house[idx_stra]; 
      is_byHouse = df_plot$is_byHouse[idx_stra]
      if (is_plot_house) { # plot household or non-household
        if (is_byHouse) { 
          sub_contacts_p2 <- sub_contacts %>% filter( str_detect(sub_contacts$is_household, "yes") ); # household pairs
        } else { 
          sub_contacts_p2 <- sub_contacts %>% filter( str_detect(sub_contacts$is_household, "no") );  # non-household pairs
        }
      } else {
        sub_contacts_p2 = sub_contacts; # plot all pairs
      }
      
      data <- sub_contacts_p2 %>% distinct(from, .keep_all=T) %>% 
        filter(delayIsol_infector < 25 & delayIsol_infector > -10) %>% pull("delayIsol_infector");
      df_pars_estim$sample_size <- length(data);
      
      result = MCMC(data, numStepsPerParameter, likelihood_options);
      
      chainRecord = result$chainRecord;
      df_pars_estim <- MCMC_posterior_analyser(chainRecord, numStepsPerParameter, likelihood_options, df_pars_estim, data)
      
      df_pars_estim_full = bind_rows( df_pars_estim_full, df_pars_estim );
    }
  }
  
  df_pars_estim_full$stratification <- rep(seq(length(df_plot$is_plot_house), 1), nwin);
  
  return(df_pars_estim_full)
}

## MCMC Fitting ----------------
df_plot = tibble(
  is_plot_house = c(T, T, F),
  is_byHouse    = c(T, F, F)
  )
onset_LB_vec = c(9, 23, 30, 9)
onset_UB_vec = c(22, 29, 44, 44)
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

df_pars_estim_full <- MCMC_IsolDelay( likelihood_options, df_plot, onset_LB_UB_Mat, Covid19.China$contacts )

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}
file_nm_write = paste0("FigS2B_", file_nm_write, "_nonOverlap_Di.csv")
write_csv( df_pars_estim_full, paste(getwd(), file_nm_write, sep = .Platform$file.sep ) )

# plotting -----

linesize = 0.5;
pd = position_dodge(0.5);

myBlack = "black";
myRed = "#d7191c";
myDarkRed = "#800026";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

ps1 <- ggplot(
  df_pars_estim_full, 
  aes(y=fitMedian, ymin=fitIQR_LB, ymax=fitIQR_UB, x=phase, group=stratification, color=factor(stratification))
  ) + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(position=pd, width=0, size=linesize+0.25) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0, 12, 4)) +
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
    values = c("#0072B2", "#E69F00", "grey40"),
    labels = c("Household", "Non-household", "No stratification")
    ) +
  labs(y="Fitted isolation delay per infector (days)") +
  coord_flip( ylim=c(0, 12)  ) +
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

plot_Fig2SAB <- plot_grid( hist3phase_IsolDelay$plot, ps1, nrow = 1, ncol = 2 )
print(plot_Fig2SAB)

ggsave2(
  paste0("FigS1B_NormalFit_forwardDi_nonoverlap.eps"), 
  plot=plot_Fig2SAB, 
  width=10 , height=12, units="cm", device="eps", dpi=300
  )


####  Fig.S2C: 14-day running time windows  ####

len_window = 13;  # To Do: Choose time window

contacts = Covid19.China$contacts;
onset_LB_vec = seq( min(contacts$date_onset_infector), max(contacts$date_onset_infector) - len_window, by=1 ) 
onset_UB_vec = onset_LB_vec + len_window
onset_LB_UB_Mat = cbind(onset_LB_vec, onset_UB_vec)

df_pars_estim_runWin <- MCMC_IsolDelay( likelihood_options, df_plot, onset_LB_UB_Mat, contacts )

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}
file_nm_write = paste0("FigS2C_", file_nm_write, "_Di_", toString(len_window+1),"DayRunWin.csv");
write_csv( df_pars_estim_runWin, paste(getwd(), file_nm_write, sep = .Platform$file.sep) );

if (likelihood_options$is_Normal) {
  file_nm_write = "NormalFit";
} else if (likelihood_options$is_Gumbel) {
  file_nm_write = "GumbelFit";
} else if (likelihood_options$is_Logistic) {
  file_nm_write = "LogisFit";
}
file_nm_write = paste0("FigS2C_", file_nm_write, "_Di_", toString(len_window+1),"DayRunWin.Rdata");
save( df_pars_estim_runWin, file = paste(getwd(), file_nm_write, sep = .Platform$file.sep ) );

# plotting -----------------
linesize = 0.5;
pd = position_dodge(0.6);

myBlack = "black";
myRed = "#d7191c";
myDarkRed = "#800026";
myBlue = "#67a9cf";
myDarkBlue = "#2b83ba";

df_pars_estim_runWin <- df_pars_estim_runWin %>% filter(onsetDay_LB >= 9)

ps2 <- ggplot(
  df_pars_estim_runWin, 
  aes(x=onsetDay_LB, y=fitMedian, ymin=fitIQR_LB, ymax=fitIQR_UB, group=stratification, color=factor(stratification))
  ) + 
  geom_point(position=pd, size=2.5) +
  geom_errorbar(position=pd, width=0, size=linesize+0.25) +
  scale_y_continuous(expand=c(0,0), breaks=seq(0, 10, 2)) +
  scale_x_continuous(
    breaks= unique(df_pars_estim_runWin$onsetDay_LB),
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
  scale_color_manual( values=c("grey40", "#E69F00", "#0072B2") ) +
  labs(y="Fitted isolation delay per infector (days)") +
  coord_cartesian( 
    xlim=c(9.5, 30.5),
    ylim=c(0, 11)
    ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.background=element_blank(),
    panel.border=element_blank(),
    axis.text.x = element_text(angle=45, size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=12),
    axis.title.x = element_blank(),
    axis.text = element_text(size=12),
    axis.line = element_line(color="black", size=linesize)
  )

X11()
plot.new()
print(ps1)

ggsave(
  "FigS2C_NormalFit_Di_14DayRunWin.eps", 
  plot=ps2, 
  width=20, height=8, units="cm", device="eps", dpi=300
  )


