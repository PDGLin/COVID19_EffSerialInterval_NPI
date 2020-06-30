
## Prepare calendar day -> No. of days since 2020/01/01

Get_calndr_str <- function(days_vec, tar_mth_str, tar_yr_str) {
  tar_date_str = sapply(
    days_vec,
    function (tar_date, mth, yr) paste(yr, mth, toString(tar_date), sep = "/"),
    mth = tar_mth_str,
    yr = tar_yr_str
  )
  return( tar_date_str );
}

Get_ndays_from_20200101 <- function( tar_yr_str ) {
  # Jan 01 to Jan 31
  days_vec1 = seq(1, 31);
  date_mth1 = Get_calndr_str(days_vec1, "1", tar_yr_str);
  # Feb 01 to Feb 29
  days_vec2 = seq(1, 29);
  date_mth2 = Get_calndr_str(days_vec2, "2", tar_yr_str);
  days_vec2 = days_vec2 + 31;
  
  df_uniq_date = tibble(
    nday = c(days_vec1, days_vec2),
    str  = as.Date( c(date_mth1, date_mth2) )
    )
  return( df_uniq_date )
}

Get_onset_day <- function(onset_date, df_uniq_date) {
  tar_vec = sapply(
    onset_date,
    function (tar_date, df_uniq_date) {
      if (str_detect(tar_date, "unknown")) {
        tar_nday = 1000000;
      } else if (str_detect(tar_date, "local")) {
        tar_nday = 200000;
      } else {
        tar_nday = df_uniq_date$nday[ df_uniq_date$str == as.Date(tar_date) ];
      }
      return(tar_nday)
    },
    df_uniq_date = df_uniq_date
  )
  return(tar_vec)
}


Set_uniq_caseID <- function( tarID_vec ) {
  result_vec = unname( sapply( tarID_vec, function (tar_id) paste0("-", tar_id, "-") ) );
  return(result_vec)
}


####  Read raw transmission pair data  ####
##  Transform into "linelist" and "contacts" as package "epicontacts" (http://www.repidemicsconsortium.org/epicontacts/)

Read_TransPairs <- function( file_name ) {
  
  df_uniq_date = Get_ndays_from_20200101( "2020" ); # Calender dates -> No. of days since 2020/01/01 
  
  raw_linelist <- read_csv( file_name );
  raw_linelist$infector_id <- Set_uniq_caseID( raw_linelist$infector_id );
  raw_linelist$infectee_id <- Set_uniq_caseID( raw_linelist$infectee_id ); 
  
  tar_linelist <- raw_linelist %>% filter(
    !str_detect(raw_linelist$infector_onsetDate, "unknown" ) & 
    !str_detect(raw_linelist$infectee_onsetDate, "unknown" )
    )
  
  tar_linelist_p1 <- tar_linelist
  
  ## Get onset times (days) of infector and infectee -----------------------------
  
  tar_linelist_p1$infector_onsetDate <- unname( Get_onset_day(tar_linelist_p1$infector_onsetDate, df_uniq_date) );
  tar_linelist_p1$infectee_onsetDate <- unname( Get_onset_day(tar_linelist_p1$infectee_onsetDate, df_uniq_date) );
  
  # infectors' onset between 9 Jan 2020 and 13 Feb 2020
  tar_linelist = tar_linelist_p1 %>% filter( infector_onsetDate >= 9 & infector_onsetDate < 45 );

  ## Isolation delay: time delay from onset to isolation of infector --------------------------
  tar_linelist$infector_isolateDate_beforeSymptom <- unname( Get_onset_day( tar_linelist$infector_isolateDate_beforeSymptom, df_uniq_date ) ); 
  tar_linelist$infector_isolateDate_afterSymptom  <- unname( Get_onset_day( tar_linelist$infector_isolateDate_afterSymptom,  df_uniq_date ) );
  tar_linelist$infector_firstHospitalVisit <- unname( Get_onset_day( tar_linelist$infector_firstHospitalVisit, df_uniq_date ) );
  tar_linelist$infector_labConfirmDate <- unname( Get_onset_day( tar_linelist$infector_labConfirmDate, df_uniq_date ) );
  
  date_isol_vec = mapply(
    function (A, B) min(A, B),
    A = tar_linelist$infector_isolateDate_beforeSymptom, 
    B = tar_linelist$infector_isolateDate_afterSymptom
    )
  tar_linelist <- tar_linelist %>% mutate(
    infector_isolateDelay = date_isol_vec - tar_linelist$infector_onsetDate
    )
  
  ## Isolation delay: time delay from onset to isolation of infectee --------------------------
  
  tar_linelist$infectee_isolateDate_beforeSymptom <- unname( Get_onset_day( tar_linelist$infectee_isolateDate_beforeSymptom, df_uniq_date ) );
  tar_linelist$infectee_isolateDate_afterSymptom  <- unname( Get_onset_day( tar_linelist$infectee_isolateDate_afterSymptom,  df_uniq_date ) );
  tar_linelist$infectee_firstHospitalVisit <- unname( Get_onset_day( tar_linelist$infectee_firstHospitalVisit, df_uniq_date ) );
  
  date2_isol_vec = mapply(
    function (A, B) min(A, B),
    A = tar_linelist$infectee_isolateDate_beforeSymptom,
    B = tar_linelist$infectee_isolateDate_afterSymptom
  )
  tar_linelist <- tar_linelist %>% mutate(
    infectee_isolateDelay = date2_isol_vec - tar_linelist$infectee_onsetDate
  )
  
  #### Transform to linelist and contacts for "epicontact" ####
  
  linelist = data.frame(
    id  = c( tar_linelist$infector_id,  tar_linelist$infectee_id ),
    age = c( tar_linelist$infector_age, tar_linelist$infectee_age ),
    sex = c( tar_linelist$infector_sex, tar_linelist$infectee_sex ),
    city_report = c( tar_linelist$infector_reportCity, tar_linelist$infectee_reportCity ),
    city_infect = c( tar_linelist$infector_infectCity, tar_linelist$infectee_infectCity ),
    date_onset  = c( tar_linelist$infector_onsetDate,  tar_linelist$infectee_onsetDate ),
    delay_isol  = c( tar_linelist$infector_isolateDelay, tar_linelist$infectee_isolateDelay )
    )
  
  linelist_uniq = distinct(linelist, id, .keep_all = T)
  
  contacts = data.frame(
    from = tar_linelist$infector_id,
      to = tar_linelist$infectee_id,
    date_onset_infector = tar_linelist$infector_onsetDate,
    date_onset_infectee = tar_linelist$infectee_onsetDate,
    serialInterval = tar_linelist$infectee_onsetDate - tar_linelist$infector_onsetDate,
    is_household = tar_linelist$isHousehold,
    delayIsol_infector = tar_linelist$infector_isolateDelay,
    delayIsol_infectee = tar_linelist$infectee_isolateDelay
    )

  newlist = list(
    "tar_linelist" = tar_linelist,
    "linelist" = linelist_uniq,
    "contacts" = contacts
    )
  
  return( newlist )
}



