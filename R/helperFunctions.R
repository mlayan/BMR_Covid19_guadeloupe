################################################################################
# Data management
################################################################################
# Transform "OUI/NON/ " character variable to "oui/non/NA" factor variable------
factor_oui_non = function(var1) {
  out = tolower(var1)
  out[out == " "] = NA
  factor(out)
}

################################################################################
# Data management
################################################################################
# Unroll dates 
unroll_days = function(df) {
  out = data.frame(
    REASEC = df$REASEC,
    Date_day = seq.Date(df$REAENT, df$REASOR, by = 1)
  ) #%>%
    #filter(row_number() <= n()-1)
  return(out)
}

################################################################################
# Intubation dates
################################################################################
getWardInfo = function(df) {
  # Convert to dates
  df = df %>%
    mutate(
      PATSTART = as_date(PATSTART), 
      PATEND = as_date(PATEND), 
      REAENT = as_date(REAENT),
      REASOR1 = as_date(REASOR1),
      REASOR2 = as_date(REASOR2),
      REASOR3 = as_date(REASOR3),
      REASOR4 = as_date(REASOR4),
      REASOR5 = as_date(REASOR5),
      REAINT1D = as_date(REAINT1D),
      REAINT1F = as_date(REAINT1F),
      REAINT2D = as_date(REAINT2D),
      REAINT2F = as_date(REAINT2F),
      REAINT3D = as_date(REAINT3D),
      REAINT3F = as_date(REAINT3F),
      REAINT4D = as_date(REAINT4D),
      REAINT4F = as_date(REAINT4F)
    )
  
  # Output dataframe
  out = data.frame(
    Date_day = seq.Date(df$PATSTART, df$PATEND, 1),
    Secteur = NA,
    present = 0,
    covid = 0,
    intub = 0
  )  
  
  # Present in the ward with secteur
  for (i in 1:5) {
    s = df[[paste0("REASOR", i)]]
    l = df[[paste0("REASEC", i)]]
    if (i == 1) e = df$REAENT
    if (i > 1) e = df[[paste0("REASOR", i-1)]]
    
    if (!is.na(s)) {
      out$present[out$Date_day >= e & out$Date_day <= s] = 1
      out$Secteur[out$Date_day >= e & out$Date_day <= s] = l
    } else {
      break
    }
  }
  
  # Intubated
  if (df$REAINT == "OUI") {
    for (i in 1:4) {
      s = df[[paste0("REAINT", i, "D")]]
      e = df[[paste0("REAINT", i, "F")]]
      
      if (!is.na(s) & !is.na(e) & sum(out$Date_day >= s & out$Date_day <= e) > 0) {
        out$intub[out$Date_day >= s & out$Date_day <= e] = 1
      } 
    }
  }
  
  # Covid-19 status
  if (df$INFCOV == "OUI") out$covid[out$present == 1] = 1
  
  return(out) 
}

################################################################################
# Data analysis
################################################################################
# New algorithms to classify episodes
classifyEpisodes = function(df, list_pvt, algo, detailed = F) {
  out = df %>% 
    distinct(REAENT, REASOR, REABCAD, REABCAE) %>%
    summarise(REAENT = min(REAENT), REASOR = max(REASOR), REABCAD = unique(REABCAD), REABCAE = unique(REABCAE))
  df = df %>% 
    arrange(PVTDAT) %>% 
    filter(PVTNAT %in% list_pvt, PVTDAT >= REAENT, PVTDAT <= REASOR) %>% 
    mutate(PVTID = match(PVTDAT, unique(PVTDAT)),
           PVTRES = ifelse(PVTRES %in% "Négatif", "negative", "positive"))
  
  if (algo == 'algo 1') out = classifyEpisodes_algo1(df, out, detailed) 
  if (algo == 'algo 2') out = classifyEpisodes_algo2(df, out, detailed) 
  
  return(out)
}

classifyEpisodes_algo1 = function(df, out, detailed = F) {
  
  within48 = df %>% filter(PVTDAT >= REAENT, PVTDAT <= REAENT+2)
  after48 = df %>% filter(PVTDAT > REAENT+2)
  if (nrow(after48)>0) m = min(after48$PVTID)
  
  # Not tested
  if (nrow(df) == 0) {
    if (out$REABCAE == "Yes" & 
        !is.na(out$REABCAD) & 
        as.numeric(difftime(out$REAENT, out$REABCAD, units = "day")) <= 7) {
      out$episode = 'Positive before' 
    } else {
      out$episode = 'Not tested'
    }
  }
  
  # Tested upon admission (within 48h)
  if (nrow(within48) > 0) {
    if (any(within48$PVTRES == "positive")) out$episode = "Positive"
    if (all(within48$PVTRES == "negative") & any(after48$PVTRES == "positive")) out$episode = "Acquisition"
    if (all(df$PVTRES == "negative")) out$episode = "Negative"    
  }
  
  # Tested after 48 h
  if (nrow(df) > 0 & nrow(within48) == 0) {
    if (out$REABCAE == "Yes" &
        !is.na(out$REABCAD) & 
        as.numeric(difftime(out$REAENT, out$REABCAD, units = "day")) <= 7) {
      out$episode = "Positive before"
    } else {
      if (any(after48$PVTRES[after48$PVTID == m] == "positive")) out$episode = "Potential positive"
      if (all(after48$PVTRES[after48$PVTID == m] == "negative") & any(after48$PVTRES == "positive")) out$episode = "Potential acquisition"
      if (all(after48$PVTRES == "negative")) out$episode = "Potential negative"  
    }
  }
  
  return(out)
}

classifyEpisodes_algo2 = function(df, out, detailed = F) {
  
  # Not tested 
  not_tested = ifelse(nrow(df) == 0, 1, 0)
  
  # Status upon admission
  within48 = df %>% filter(PVTDAT >= REAENT, PVTDAT <= REAENT+2)
  
  status_admission = 0
  if (nrow(within48) > 0) status_admission = 1
  if (nrow(within48) == 0 & !is.na(out$REABCAD) & difftime(out$REAENT, out$REABCAD, units = "days") <= 7) {
    status_admission = 1
    df = df %>%
      add_row(REAENT = out$REAENT, REASOR = out$REASOR, PVTDAT = out$REAENT, PVTRES = "positive", .before = 1)
  }
  
  # Test sequence smoothing
  intern_samples = c("Aspiration bronchique", "Aspiration duodénale", "Aspiration endo-trachéale",
                     "LBA", "LCR", "Liquide d'ascite", "Liquide de ponction articulaire ou synoviale",
                     "Liquide péricardique", "Liquide péritonéal", "Liquide pleural", "PDP")
  other_samples = c("Cathéter", "Cornée", "Crachat", "Drain", "Gorge", "Nasopharynx", "Prélèvement ORL",
                    "Pus", "Sang cathéter", "Selles", "Sonde", "Urine", "Vaginal")
  
  if ("REASEC" %in% colnames(df)) df = df %>% group_by(REAENT, REASOR, PVTDAT, REASEC)
  if (!"REASEC" %in% colnames(df)) df = df %>% group_by(REAENT, REASOR, PVTDAT)
  df = df %>%
    # Multiple tests on the same day --> Positive tests dominate negative tests 
    summarise(
      PVTRES = case_when(
        all(PVTRES %in% "negative") ~ "negative",
        any(PVTRES %in% "positive") ~ "positive",
        all(is.na(PVTRES)) ~ NA,
      .default = "?"
      ), 
      PVTNAT = case_when(
        all(PVTNAT %in% "Pvt rectal") ~ "Pvt rectal",
        any(PVTNAT %in% "Sang périphérique") ~ "Sang périphérique",
        any(PVTNAT %in% intern_samples) & all(!PVTNAT %in% "Sang périphérique") ~ "intern",
        any(PVTNAT %in% other_samples) & all(!PVTNAT %in% c("Sang périphérique", intern_samples)) ~ "other",
        .default = "before study"
          ),
      .groups = "drop"
      ) %>%
    # Smooth sequence
    mutate(
      PVTRES = case_when(
        is.na(lag(PVTRES)) ~ PVTRES,
        (PVTRES == "negative" & lag(PVTRES, 1) == "positive" & is.na(lead(PVTRES, 1)) ) | 
          (PVTRES == "negative" & lag(PVTRES, 1) == "positive" & lead(PVTRES, 1) == "negative" & is.na(lead(PVTRES, 2)) ) |
          (PVTRES == "negative" & lag(PVTRES, 1) == "positive" & lead(PVTRES, 1) == "positive" ) | 
          (PVTRES == "negative" & lag(PVTRES, 1) == "positive" & lead(PVTRES, 1) == "negative" &  lead(PVTRES, 2) == "positive" ) ~ "positive",
        .default = PVTRES
      )
    )
  
  # Classification of the episode
  episode = NA
  split_res = split(df$PVTRES, rleid(df$PVTRES))
  if (status_admission) {
    # Clear episodes
    if (all(df$PVTRES == "positive")) episode = "Positive"
    if (all(df$PVTRES == "negative")) episode = "Negative"
    if (length(split_res) == 2 & all(split_res[[1]] == "negative")) episode = "Acquisition"
    if (length(split_res) == 2 & all(split_res[[1]] == "positive")) episode = "Decolonization"
    if (length(split_res) > 2) episode = "Complex"
  }
  if (!not_tested & !status_admission) {
    # Potential episodes
    if (all(df$PVTRES == "positive")) episode = "Potential positive"
    if (all(df$PVTRES == "negative")) episode = "Potential negative"
    if (length(split_res) == 2 & all(split_res[[1]] == "negative")) episode = "Potential acquisition"
    if (length(split_res) == 2 & all(split_res[[1]] == "positive")) episode = "Potential decolonization"
    if (length(split_res) > 2) episode = "Potential complex"
  }
  if (not_tested & !status_admission) episode = "Not tested"
  if (not_tested & status_admission) episode = "Positive"
  
  out$episode = episode
  if (detailed) {
    if (nrow(df) == 0) {
      df = df %>% add_row(REAENT = out$REAENT, REASOR = out$REASOR)
      return(df)
    } else {
      return(df)
    }
  }
  if (!detailed) return(out)
}

# Get episode type
# 0: no test
# 1: only negative tests
# 2: only positive tests
# 3: acquisition
# 4: clearance
# 5: more complicated 
summarizeTestResults = function(df) {
  if (nrow(df) == 0) {
    out = data.frame(test_results = 0, first_pos_before = NA, pos_inc = NA)
    
  } else {
    df_temp = df %>%
      arrange(PVTDAT) %>%
      mutate(result_lag = lag(result)) 
    result = df_temp$result
    result_lag = df_temp$result_lag
    
    out = data.frame()
    if(all(result %in% NA)) {
      out = data.frame(test_results = 0, first_pos_before = NA, pos_inc = NA)
    }
    if(all(result %in% 0)) {
      out = data.frame(test_results = 1, first_pos_before = NA, pos_inc = NA)
    }
    if(all(result %in% 1)) {
      out = data.frame(test_results = 2, first_pos_before = NA, pos_inc = NA)
    }
    if(1 %in% result & 0 %in% result & 
         (sum((result - result_lag) %in% 1) == 1 & sum((result - result_lag) %in% -1) == 0)
    )
       {
      fpb = df_temp$PVTDAT[result == 1][1] <= unique(df_temp$PATEND)
      out = data.frame(test_results = 3, first_pos_before = fpb, pos_inc = NA)
    }
    
    if (0 %in% result & 1 %in% result & (
      (sum((result - result_lag) %in% 1) == 1 & sum((result - result_lag) %in% -1) == 1) |
      (sum((result - result_lag) %in% 1) == 0 & sum((result - result_lag) %in% -1) == 1)
      )
      ){
      pi = any(unique(df_temp$PATSTART) <= df_temp$PVTDAT[result == 1] & df_temp$PVTDAT[result == 1] <= unique(df_temp$PATEND))
      out = data.frame(test_results = 4, first_pos_before = NA, pos_inc = pi)
    }
    
    if (nrow(out) == 0) {
      out = data.frame(test_results = 5, first_pos_before = NA, pos_inc = NA)
    }
  }
  return(out)
}

# Classify patients according to test result and timing of test
getResults = function(df) {
  out = data.frame()
  
  if (any(is.na(df$PVTDAT))) {
    out = data.frame(PVTRECTAL = NA, RESULT = NA, PVTDAT = NA)
    
  } else {
    pos_results = df %>% 
      arrange(PVTDAT) %>%
      filter(Result %in% 1)
    
    start_inc = unique(df$PATSTART)
    end_inc = unique(df$PATEND)
    rea_ent = unique(df$REAENT)
    
    if (nrow(pos_results) == 0) {
      # Individuals testing only negative
      pvtrectal = ifelse(any(df$PVTNAT %in% "Pvt rectal"), 1, 0)
      out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Negative", PVTDAT = NA)
      
    } else {
      if (any(pos_results$PVTDAT <= rea_ent + 2)) {
        
        # Individuals with positive test within 48 H after admission at the ICU
        within_48h = min(pos_results$PVTDAT[pos_results$PVTDAT <= rea_ent + 2])
        within_48h_nat = pos_results$PVTNAT[pos_results$PVTDAT <= rea_ent + 2]
        neg_results = df[df$PVTDAT < within_48h & df$PVTDAT >= start_inc & df$PVTNAT %in% within_48h_nat, ]
        
        if (nrow(neg_results) == 0) {
          pvtrectal = ifelse(any(pos_results$PVTNAT[pos_results$PVTDAT <= rea_ent + 2] %in% "Pvt rectal"), 1, 0)
          out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Positive", PVTDAT = NA)
        } else {
          pvtrectal = ifelse(any(pos_results$PVTNAT[pos_results$PVTDAT <= rea_ent + 2] %in% "Pvt rectal"), 1, 0)
          out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Acquisition", PVTDAT = within_48h)
        }
        
      } else {
        
        # Individuals with positive test 48 h after ICU admission
        first_positive = min(pos_results$PVTDAT)
        pvtrectal = ifelse(pos_results$PVTNAT[1] %in% "Pvt rectal", 1, 0)
        
        if (first_positive <= start_inc) {
          out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Positive", PVTDAT = NA)
        } else if (first_positive > end_inc) {
          out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Negative", PVTDAT = NA)
        } else {
          out = data.frame(PVTRECTAL = pvtrectal, RESULT = "Acquisition", PVTDAT = first_positive)
        }
        
      }
    }
  }
  
  return(out)
}

################################################################################
# Data simulation
################################################################################
# Classification of episodes----------------------------------------------------
getEpisodeCategory = function(df) {
  if (all(is.na(df$value)) | all(df$value %in% 0)) {
    return("negative")
  } else {
    # Smoothing
    # If two positive tests separated by up to 2 negative tests
    # the negative tests are considered positive
    df = df %>%
      arrange(PVTDAT) %>%
      replace_na(list(value = 0)) %>%
      mutate(
        value_smoothed = value,
        test1 = as.numeric(lag(value) == 1 & lead(value) == 1 & difftime(PVTDAT, lag(PVTDAT), units = "days") <= 7 & difftime(lead(PVTDAT), PVTDAT, units = "days") <= 7),
        test2 = as.numeric(lag(value) == 1 & lead(value,2) == 1 & difftime(PVTDAT, lag(PVTDAT), units = "days") <= 7 & difftime(lead(PVTDAT,2), PVTDAT, units = "days") <= 14), 
        test3 = as.numeric(lag(value,2) == 1 & lead(value) == 1 & difftime(PVTDAT, lag(PVTDAT,2), units = "days") <= 14 & difftime(lead(PVTDAT), PVTDAT, units = "days") <= 7)
        ) %>%
      replace_na(list(test1 = 0, test2 = 0, test3 = 0)) %>%
      mutate(value_smoothed = ifelse(
        value %in% 0 & (test1 | test2 | test3),
        1,
        value_smoothed
      )
      )
    
    # Classification
    if (nrow(df) == 1 ) {
      if (df$value_smoothed == 1) return("positive")
      if (df$value_smoothed == 0) return("negative")
    } else if (nrow(df) == 2) {
      if (all(df$value_smoothed %in% 1)) return("positive")
      if (all(df$value_smoothed %in% 0)) return("negative")
      if (!all(df$value_smoothed %in% 1)) return("clearance")
    }else {
      acquisition_clearance_events = df$value_smoothed - lag(df$value_smoothed)
      if (sum(acquisition_clearance_events %in% 1) == 1 & sum(acquisition_clearance_events %in% -1) == 0) return("positive")
      if (any(acquisition_clearance_events %in% 1) == 1 & sum(acquisition_clearance_events %in% -1) == 1) return("clearance")
      if (sum(acquisition_clearance_events %in% 1) > 1) return("successive acquisitions")      
    }
  }
}


# Draw random stay length in ICU from discretized lognormal distribution--------
rStayLength = function(n) {
  
  # Upper discretization: px = F(x+h) - F(x)
  discretized_rlnorm = plnorm((0:100)+1, meanlog = 3, sdlog = 1.6) - plnorm((0:100), meanlog = 3, sdlog = 1.6)
  discretized_rlnorm = discretized_rlnorm/sum(discretized_rlnorm)
  
  out = sample(x = 0:100, size = n, replace = T, prob = discretized_rlnorm)
  return(out)
}


# Draw prevalence data from binomial distributions------------------------------
rPrevalence = function(nDays, careOrganization, prevType, occupancy_data = NULL) {
  nSectors = length(careOrganization)
  prev = c()
  
  if (prevType == "covid") {
    prev = rep(0, nDays*length(careOrganization))
  }
  
  if (prevType == "intubation" | prevType == "dialysis") {
    
    if (prevType == "intubation") p = 0.7
    if (prevType == "dialysis") p = 0.4
    
    occupancy_data$event = sample(c(0,1), nrow(occupancy_data), replace = T, prob = c(1-p, p))
    occupancy_data$event_date = apply(occupancy_data, 1, function(x) {
      out = NA
      if (x[["event"]] == 1) {
        out = sample(x[["admission"]]:x[["discharge"]], 1, replace = F)
      }
      return(out)
    })
    
    prev = c()
    for (s in 0:(nSectors-1)) {
      for (d in 0:(nDays-1)) {
        prev_temp = sum(occupancy_data$subsector == s & 
                          occupancy_data$event == 1 & 
                          occupancy_data$event_date <= d & 
                          occupancy_data$discharge >= d &
                          occupancy_data$admission <= d
                        )
        prev_temp = round(prev_temp/careOrganization[s+1], 3)
        prev = c(prev, prev_temp) 
      }
    }
  }
  
  # data.frame(
  #   period = 0,
  #   subsector = rep(0:(nSectors-1), each = nDays),
  #   t = rep(0:(nDays-1), nSectors),
  #   val = prev
  # ) %>%
  #   write.table(., paste0("data/simulated/", prevType, "_data_sim.txt"), 
  #               quote = F, sep = " ", col.names = F, row.names = F)
  
  out = matrix(prev, nrow = nSectors, ncol = nDays, byrow = T)
  
  return(out)
}


# Transform simulation output into mcmc input-----------------------------------
getEpisode = function(df, nDays) {
  
  # Get colonization type
  if (df[["detection_date"]] < 0 | df[["detection_date"]] > df[["discharge"]]) colType = "neg"
  if (df[["detection_date"]] >= 0 & df[["detection_date"]] == df[["admission"]]) colType = "pos"
  if (df[["detection_date"]] >= 0 & 
      df[["detection_date"]] > df[["admission"]] & 
      df[["detection_date"]] <= df[["discharge"]]) colType = "acq"
  
  # Get status sequence
  status = rep(-1, nDays)
  status[(df[["admission"]]+1):(df[["discharge"]]+1)] = 0
  if (colType != "neg") {
    status[(df[["detection_date"]]+1):(df[["discharge"]]+1)] = 1
  }
  
  # Output in dataframe format
  out = data.frame(
    colType = colType,
    status = status
  )
  return(out)
}


#############################################
# Functions for correlation plots
#############################################
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0, 1, 1, alpha = 0.5), ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- cor(x, y)
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * abs(Cor)) # Resize the text by level of correlation
}