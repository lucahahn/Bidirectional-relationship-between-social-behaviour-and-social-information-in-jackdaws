
#Low quality feeders: 48 days, day 170 - 218, 18/06 to 05/08
#High quality feeders: 37 days, day 180 - 217, 28/06 to 04/08

# (1) IMPORT DATA  ----

library(dplyr)
library(stringi)
library(stringr)
library(svMisc)
library(tidyr)
library(tidyverse)

library(asnipe)
library(car)
library(carData)
library(chisq.posthoc.test)
library(ClusterR)
library(corrplot)
library(data.table)
library(DHARMa)
library(dplyr)
library(emmeans)
library(extrafont)
font_import()
loadfonts()
library(ggplot2)
library(glmmTMB)
library(gridGraphics)
library(grid)
library(gridBase)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(ggeffects)
library(ggdist)
library(hms)
library(igraph)
library(lme4)
library(lmerTest)
library(lubridate)
library(MASS)
library(multcomp)
library(NBDA)
library(performance)
library(psych)
library(RColorBrewer)
library(reshape2)
library(rptR)
library(scales)
library(stringi)
library(stringr)
library(svMisc)
library(tidyr)
library(tidyverse)
library(tidybayes)
library(assortnet)

library(brms)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(gdata)
library(dplyr)
library(parallel)
library(cowplot)

library(STbayes)
library(dplyr)

library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)
library(distributional)
install.packages("viridis")  #Install
library("viridis")

#Shows milliseconds
op <- options(digits.secs=1)

#Use directory where you want to look for RT files, concatenate paths
RT_files <- list.files("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Knowledge", pattern = "RT", recursive = TRUE)
RT_paths <- paste("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Knowledge",RT_files, sep = "/")

#Load saved life history csv file 
LH <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/LH20240716.csv", header = T, stringsAsFactors = F, fileEncoding="latin1")
LH$DATE <- strptime(LH$DATE,format = "%d/%m/%Y")
LH$DATE <- as.Date(LH$DATE, format = "%d/%m/%Y") # convert to date

#Load list with "knowledgeable" individuals
all_individuals_knowledgeable <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/all_individuals_knowledgeable.csv", header = T, stringsAsFactors = F, fileEncoding="latin1")
all_individuals_knowledgeable$age <- all_individuals$age[match(all_individuals_knowledgeable$JID, all_individuals$JID)]
all_individuals_knowledgeable$sex <- all_individuals$sex[match(all_individuals_knowledgeable$JID, all_individuals$JID)]
all_individuals_knowledgeable$knowledge <- "knowledgeable"

all_individuals$knowledgeable <- all_individuals_knowledgeable$knowledge[match(all_individuals$JID, all_individuals_knowledgeable$JID)]
all_individuals <- all_individuals %>%
  mutate(knowledge=replace_na(knowledge, "naive"))

#Create sub-strings that contain feeder ID (e.g. "Y1.1") and date
day_arrays <- unique(substr(RT_files,5,15))

#Create empty list to place data into as we go
collapsed_list <- list()

#Run through RT files in list 
for(i in 1:length(day_arrays)) {
  progress(i, max.value = length(day_arrays))
  
  day_array_list <- list()
  
  for (j in 1:length(RT_files[which(substr(RT_files,5,15) == day_arrays[i])])) {
    temp_day_file <- (read.delim(RT_paths[which(substr(RT_files,5,15) == day_arrays[i])][j], header = T, stringsAsFactors = F))[,1:13]
    temp_day_file$feeder <- substr(RT_files[which(substr(RT_files,5,15) == day_arrays[i])][j],5,8)
    day_array_list[[j]] <- temp_day_file
  }
  
  temp_RT <- do.call(rbind, day_array_list)
  
  temp_RT$Time <- strptime(paste(temp_RT$Date,stri_sub(temp_RT$Hmsec/1024,2,5), sep = ""), "%Y-%m-%d %H:%M:%OS")  # Add in miliseconds (1024 in a second) and format time
  
  temp_RT %>% filter(nchar(TagID_hex) == 10) -> temp_tags  #remove times when no tag
  
  
  if(dim(temp_tags)[1] >0){  
    
    visits <- data.frame(Event = temp_tags$Event, Start = temp_tags$Time, End = temp_tags$Time+(0.5*(temp_tags$Reps -1)), tag = temp_tags$TagID_hex, feeder = temp_tags$feeder)  # adds reps to visit length (0.5 seconds for every extra detection as that was resampling speed)
    visits <- arrange(visits, Start)
    visits <- arrange(visits, feeder)
    within_errors <- which(visits$End < lag(visits$End) & visits$tag == lag(visits$tag) & visits$feeder == lag(visits$feeder))
    if(length(within_errors) > 0){
      visits <- visits[-which(visits$End < lag(visits$End) & visits$tag == lag(visits$tag) & visits$feeder == lag(visits$feeder)),]  ## Get rid of reads within bouts - almost always erroneous single reads that are repeated later in the datastream
    }
    visits <- arrange(visits, Start)
    visits$count <- sapply(1:nrow(visits),function(x)sum(visits$tag[x]==visits$tag[1:x]))  ## add individual visit counter - if two bouts don't have sequential counts then bird seen elsewhere in between
    visits <- arrange(visits, feeder)
    
    visits$collapse <- "0"  # Temp column to tell me if this bout is to be collapsed
    
    visit_time <- 10 ## How long between detections before a new bout is classed
    
    # Is the last visit ending within 'visit_time' seconds of this one starting, and with the same tag & in sequence?
    visits$collapse[ which (visits$Start-lag(visits$End) < visit_time & lead(visits$Start) - visits$End < visit_time & visits$feeder == lag(visits$feeder) & visits$feeder == lead(visits$feeder) & visits$tag == lag(visits$tag) & visits$tag == lead(visits$tag) & (visits$count-lag(visits$count)) == 1 & (lead(visits$count)-visits$count) == 1 )] <- "1"  ## What about times they hop in between?!
    
    visits$collapse[which(visits$collapse == 0 & lag(visits$collapse == 1))] <- "End"  ## Can work out the end based on the 0s and 1s
    visits$collapse[which(visits$collapse == 0 & lead(visits$collapse == 1))] <- "Start" ## Can work out the start based on the 0s and 1s
    
    ## Mop up those that only have two potential detections (so no middle values to get assigned 1 above)
    visits$collapse[which(visits$collapse == 0 & lead(visits$Start) - visits$End < visit_time & lead(visits$Start) - visits$End > -2 & visits$tag == lead(visits$tag) & lead(visits$count) - visits$count == 1)] <- "Start"
    visits$collapse[which(visits$collapse == 0 & visits$Start - lag(visits$End) < visit_time & lag(visits$feeder) == visits$feeder & visits$tag == lag(visits$tag) & lag(visits$count) - visits$count == -1)] <- "End"
    
    # Add this file's data to the list
    collapsed_list[[i]] <- data.frame(start = visits$Start[which(visits$collapse %in% c("Start","0"))], end = visits$End[which(visits$collapse %in% c("End","0"))], tag = visits$tag[which(visits$collapse %in% c("Start","0"))], event = visits$Event[which(visits$collapse %in% c("Start","0"))], feeder =  visits$feeder[which(visits$collapse %in% c("Start","0"))], array = rep(stri_sub(day_arrays[i], 8,11), length(which(visits$collapse %in% c("Start","0")))))  }
}

#Turn the list into a data frame
visit_data <- do.call(rbind,collapsed_list)

#Find instances in which time is still NA
visit_data_NA <- subset(visit_data, is.na(visit_data$start))

#Filter instances in which time is not NA
visit_data <- subset(visit_data, !is.na(visit_data$start))

#Remove instances where JID = NA and test tags
visit_data$JID <- LH_RFID$ID[match(visit_data$tag,LH_RFID$RFID)]
visit_data <- subset(visit_data, JID != "NA")

#Read visit data ----

#Shows milliseconds
op <- options(digits.secs=1)

#Use directory where you want to look for RT files, concatenate paths
RT_files <- list.files("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Knowledge", pattern = "RT", recursive = TRUE)
RT_paths <- paste("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Knowledge",RT_files, sep = "/")

library(dplyr)
library(stringi)

# Unique days/positions
day_positions <- unique(substr(RT_files, 6, 16))
collapsed_list <- list()
visit_time <- 10  # seconds, gap to start new visit

for (i in seq_along(day_positions)) {
  
  progress(i, max.value = length(day_positions))
  
  day_files_idx <- which(substr(RT_files, 6, 16) == day_positions[i])
  
  # Read all files for this day
  day_position_list <- lapply(day_files_idx, function(j) {
    temp_day_file <- read.delim(RT_paths[j], header = TRUE, stringsAsFactors = FALSE)[,1:13]
    temp_day_file$feeder <- substr(RT_files[j], 6, 9)
    temp_day_file
  })
  
  temp_RT <- do.call(rbind, day_position_list)
  
  # Ensure Hmsec is numeric, replace NA with 0
  temp_RT$Hmsec <- as.numeric(temp_RT$Hmsec)
  temp_RT$Hmsec[is.na(temp_RT$Hmsec)] <- 0
  
  # Pad 2-digit Hmsec by multiplying by 10
  pad_idx <- which(temp_RT$Hmsec >= 10 & temp_RT$Hmsec < 100)
  temp_RT$Hmsec[pad_idx] <- temp_RT$Hmsec[pad_idx] * 10
  
  # Build exact timestamps including milliseconds
  temp_RT$Time_exact <- as.POSIXct(temp_RT$Date, format="%Y-%m-%d %H:%M:%S", tz="UTC") +
    (temp_RT$Hmsec / 1024)
  
  # Keep only valid tags
  temp_tags <- temp_RT %>% filter(nchar(TagID_hex) == 10)
  if (nrow(temp_tags) == 0) next
  
  # Create read intervals with 0.5*(Reps - 1) adjustment
  reads <- data.frame(
    Time  = temp_tags$Time_exact,
    End   = temp_tags$Time_exact + 0.5 * (temp_tags$Reps - 1),
    tag   = temp_tags$TagID_hex,
    feeder= temp_tags$feeder,
    event = temp_tags$Event
  )
  
  reads <- reads %>% arrange(feeder, Time)
  
  visits_all <- list()
  visit_counter <- 1
  
  # Loop over feeders
  for (f in unique(reads$feeder)) {
    feeder_reads <- reads %>% filter(feeder == f) %>% arrange(Time)
    if (nrow(feeder_reads) == 0) next
    
    current_tag   <- feeder_reads$tag[1]
    current_start <- feeder_reads$Time[1]
    current_end   <- feeder_reads$End[1]
    current_event <- feeder_reads$event[1]
    
    for (r in 2:nrow(feeder_reads)) {
      this_tag   <- feeder_reads$tag[r]
      this_time  <- feeder_reads$Time[r]
      this_end   <- feeder_reads$End[r]
      this_event <- feeder_reads$event[r]
      
      if (any(is.na(c(this_tag, this_time, this_end)))) next
      
      gap <- as.numeric(difftime(this_time, current_end, units = "secs"))
      
      # New visit condition
      if (this_tag != current_tag || gap > visit_time) {
        
        # Save previous visit
        visits_all[[visit_counter]] <- data.frame(
          feeder = f,
          tag    = current_tag,
          start  = current_start,
          end    = current_end,
          event  = current_event
        )
        visit_counter <- visit_counter + 1
        
        # Start new visit
        current_tag   <- this_tag
        current_start <- this_time
        current_end   <- this_end
        current_event <- this_event
        
      } else {
        # Continue same visit
        current_end <- max(current_end, this_end)
      }
    }
    
    # Save last visit for feeder
    visits_all[[visit_counter]] <- data.frame(
      feeder = f,
      tag    = current_tag,
      start  = current_start,
      end    = current_end,
      event  = current_event
    )
    visit_counter <- visit_counter + 1
  }
  
  collapsed_list[[i]] <- dplyr::bind_rows(visits_all)
}

# Combine all visits
visit_data <- dplyr::bind_rows(collapsed_list)

setDT(visit_data)

# Ensure timestamps are POSIXct with milliseconds
visit_data[, start := ymd_hms(start)]
visit_data[, end   := ymd_hms(end)]
options(digits.secs = 3)

# Assign a unique event ID per visit
visit_data[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), format(end, "%Y-%m-%d %H:%M:%OS3"), tag, sep = " ")]

#Identify duplicates
dup_rows <- visit_data %>%
  group_by(event_id) %>%
  filter(n() > 1) %>%
  ungroup()

#Remove duplicates
visit_data <- visit_data %>%
  distinct(event_id, .keep_all = TRUE)

# Ensure milliseconds are shown
options(digits.secs = 3)

# Check result
head(visit_data$start)
head(visit_data$end)

#Fix remaining instances of visit overlap
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(
    next_start = lead(start),
    end = dplyr::if_else(
      !is.na(next_start) & end > next_start,
      next_start,
      as.POSIXct(end, origin = "1970-01-01")
    )
  ) %>%
  dplyr::select(-next_start) %>%
  ungroup()

# Optional: compute time differences between consecutive visits per feeder
visit_data <- visit_data %>%
  arrange(feeder, start) %>%
  group_by(feeder) %>%
  mutate(time_diff_pre = as.numeric(start - lag(end))) %>%
  ungroup()

sum(visit_data$time_diff_pre < 0, na.rm = T)

#Find instances in which time is still NA
visit_data_NA <- subset(visit_data, is.na(visit_data$start))

#Filter instances in which time is not NA
visit_data <- subset(visit_data, !is.na(visit_data$start))

#Remove instances where JID = NA and test tags
visit_data$JID <- LH_RFID$ID[match(visit_data$tag, LH_RFID$RFID)]
visit_data <- subset(visit_data, JID != "NA")


# (2) ADD DATA ----

#Adding information about visit duration
visit_data$interval <- interval(visit_data$start,visit_data$end)
visit_data$visit_duration <- as.duration(visit_data$interval)

#Adding information about site: Y, Z
visit_data[substr(visit_data$feeder,1,1)=="Y","site"]<-"Y"
visit_data[substr(visit_data$feeder,1,1)=="Z","site"]<-"Z"

#Adding information about feeder quality
visit_data[substr(visit_data$feeder,2,3)=="HK","quality"]<-"high"
visit_data[substr(visit_data$feeder,2,3)=="LK","quality"]<-"low"

visit_data_knowledge <- visit_data
visit_data_knowledge2 <- visit_data_knowledge

#Add information about individuals
visit_data_knowledge$age <- all_individuals$age[match(visit_data_knowledge$JID, all_individuals$JID)]
visit_data_knowledge$sex <- all_individuals$sex[match(visit_data_knowledge$JID, all_individuals$JID)]

#Add day
visit_data_knowledge$day <- yday(visit_data_knowledge$start) #day of year

#Unique individuals per site
unique(visit_data_knowledge$JID)
unique(visit_data_knowledge$JID)
unique(visit_data_knowledge$JID)

visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="YLK1" | substr(visit_data_knowledge$feeder,1,4)=="YHK1","position"]<-"Y1"
visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="YLK2" | substr(visit_data_knowledge$feeder,1,4)=="YHK2","position"]<-"Y2"
visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="YLK3" | substr(visit_data_knowledge$feeder,1,4)=="YHK3","position"]<-"Y3"
visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="YLK4" | substr(visit_data_knowledge$feeder,1,4)=="YHK4","position"]<-"Y4"

visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="ZLK1" | substr(visit_data_knowledge$feeder,1,4)=="ZHK1","position"]<-"Z1"
visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="ZLK2" | substr(visit_data_knowledge$feeder,1,4)=="ZHK2","position"] <-"Z2"
visit_data_knowledge[substr(visit_data_knowledge$feeder,1,4)=="ZLK3" | substr(visit_data_knowledge$feeder,1,4)=="ZHK3","position"]<-"Z3"

visit_data_knowledge <- subset(visit_data_knowledge, !visit_data_knowledge$feeder == "ZLK3")
visit_data_knowledge <- subset(visit_data_knowledge, !visit_data_knowledge$feeder == "YLK2")
visit_data_knowledge <- subset(visit_data_knowledge, !visit_data_knowledge$feeder == "YLK4")

table(visit_data_knowledge$site)
table(visit_data_knowledge$position)
table(visit_data_knowledge$quality)

#Subset data per site
visit_data_knowledge_Y <- subset(visit_data_knowledge, site == "Y")
visit_data_knowledge_Z <- subset(visit_data_knowledge, site == "Z")

#Subset per feeder quality
visit_data_knowledge_LK <- subset(visit_data_knowledge, quality == "low")
visit_data_knowledge_LK_pre <- subset(visit_data_knowledge_LK, day < 180)
visit_data_knowledge_LK_post <- subset(visit_data_knowledge_LK, day > 179)

visit_data_knowledge_HK <- subset(visit_data_knowledge, quality == "high")

visit_data_knowledge_Y1 <- subset(visit_data_knowledge, position == "Y1")
visit_data_knowledge_Y2 <- subset(visit_data_knowledge, position == "Y2")
visit_data_knowledge_Y3 <- subset(visit_data_knowledge, position == "Y3")
visit_data_knowledge_Y4 <- subset(visit_data_knowledge, position == "Y4")

visit_data_knowledge_Z1 <- subset(visit_data_knowledge, position == "Z1")
visit_data_knowledge_Z2 <- subset(visit_data_knowledge, position == "Z2")
visit_data_knowledge_Z3 <- subset(visit_data_knowledge, position == "Z3")

unique(visit_data_knowledge_Z2$JID)

unique_ID = unique(visit_data_knowledge$JID)
all_individuals_lowq = data.frame(JID = unique_ID)
all_individuals_lowq$lowq <- "lowq"

#Visit number per individual data set
visit_number_knowledge <- as.data.frame(table(visit_data_knowledge$JID))
visit_number_knowledge$JID <- visit_number_knowledge$Var1
visit_number_knowledge$visit_number_knowledge <- visit_number_knowledge$Freq
visit_number_knowledge <- subset(visit_number_knowledge, select = -c(Var1, Freq))

visit_number_knowledge$age <- Ringed$min_age[match(visit_number_knowledge$JID, Ringed$ID)]

visit_number_knowledge$RFID <- LH_RFID$RFID[match(visit_number_knowledge$JID, LH_RFID$ID)]

visit_number_knowledge_post <- as.data.frame(table(visit_data_knowledge_post$JID))
visit_number_knowledge_post$JID <- visit_number_knowledge_post$Var1
visit_number_knowledge_post$visit_number_knowledge <- visit_number_knowledge_post$Freq
visit_number_knowledge_post <- subset(visit_number_knowledge_post, select = -c(Var1, Freq))

visit_number_knowledge_HK <- as.data.frame(table(visit_data_knowledge_HK$JID))
visit_number_knowledge_HK$JID <- visit_number_knowledge_HK$Var1
visit_number_knowledge_HK$visit_number <- visit_number_knowledge_HK$Freq
visit_number_knowledge_HK <- subset(visit_number_knowledge_HK, select = -c(Var1, Freq))
visit_number_knowledge_HK$knowledge_phase <- "post_knowledge"
visit_number_knowledge_HK$JID_phase <- paste(visit_number_knowledge_HK$JID, visit_number_knowledge_HK$knowledge_phase, sep = " ")

visit_number_knowledge_LK <- as.data.frame(table(visit_data_knowledge_LK$JID))
visit_number_knowledge_LK$JID <- visit_number_knowledge_LK$Var1
visit_number_knowledge_LK$visit_number <- visit_number_knowledge_LK$Freq
visit_number_knowledge_LK <- subset(visit_number_knowledge_LK, select = -c(Var1, Freq))
visit_number_knowledge_LK$knowledge_phase <- "post_knowledge"
visit_number_knowledge_LK$JID_phase <- paste(visit_number_knowledge_LK$JID, visit_number_knowledge_LK$knowledge_phase, sep = " ")

visit_number_knowledge_LK_pre <- as.data.frame(table(visit_data_knowledge_LK_pre$JID))
visit_number_knowledge_LK_pre$JID <- visit_number_knowledge_LK_pre$Var1
visit_number_knowledge_LK_pre$visit_number <- visit_number_knowledge_LK_pre$Freq
visit_number_knowledge_LK_pre <- subset(visit_number_knowledge_LK_pre, select = -c(Var1, Freq))
visit_number_knowledge_LK_pre$knowledge_phase <- "pre_knowledge"
visit_number_knowledge_LK_pre$JID_phase <- paste(visit_number_knowledge_LK_pre$JID, visit_number_knowledge_LK_pre$knowledge_phase, sep = " ")

visit_number_knowledge_LK_post <- as.data.frame(table(visit_data_knowledge_LK_post$JID))
visit_number_knowledge_LK_post$JID <- visit_number_knowledge_LK_post$Var1
visit_number_knowledge_LK_post$visit_number <- visit_number_knowledge_LK_post$Freq
visit_number_knowledge_LK_post <- subset(visit_number_knowledge_LK_post, select = -c(Var1, Freq))
visit_number_knowledge_LK_post$knowledge_phase <- "post_knowledge"
visit_number_knowledge_LK_post$JID_phase <- paste(visit_number_knowledge_LK_post$JID, visit_number_knowledge_LK_post$knowledge_phase, sep = " ")

#Visits per individual per site and per feeder, preferred sites
perindivpersite_knowledge <- visit_data_knowledge  %>%  count(JID, site)
perindivpersite_knowledge2 <- reshape(perindivpersite_knowledge, idvar = "JID", timevar = "site", direction = "wide")
perindivpersite_knowledge3 <- perindivpersite_knowledge2
perindivpersite_knowledge3[is.na(perindivpersite_knowledge3)] <- 0.99
perindivpersite_knowledge3$pref <- perindivpersite_knowledge3$n.Y / perindivpersite_knowledge3$n.Z
perindivpersite_knowledge3 <- perindivpersite_knowledge3 %>% 
  as_tibble() %>% 
  mutate(pref.site = if_else(pref > 0.5,"Y", "Z"))

visit_number_knowledge$pref_site <- perindivpersite_knowledge3$pref.site[match(visit_number_knowledge$JID, perindivpersite_knowledge3$JID)]

all_individuals_knowledge <- visit_number_knowledge

all_individuals_knowledge$knowledgeable <- all_individuals_knowledgeable$knowledge[match(all_individuals_knowledge$JID, all_individuals_knowledgeable$JID)]
all_individuals_knowledge <- all_individuals_knowledge %>%
  mutate(knowledge=replace_na(knowledgeable, "naive")) 

#Create a copy for subsetting
all_individuals_knowledge_sub <- all_individuals_knowledge

#Subset by number of visits at experimental setup
#all_individuals_knowledge_sub <- subset(all_individuals_knowledge, all_individuals_knowledge$visit_number > 5)
#all_individuals_knowledge_sub <- subset(all_individuals_knowledge, all_individuals_knowledge$visit_number > 10)

visit_data_knowledge$knowledge <- all_individuals_knowledge$knowledge[match(visit_data_knowledge$JID, all_individuals_knowledge$JID)]
visit_data_knowledge$quality_knowledge <- paste(visit_data_knowledge$quality, visit_data_knowledge$knowledge, sep = "_")

#Data for manuscript ----

write.csv(all_individuals_knowledge_sub4_long,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_knowledge_sub4_long.csv", row.names = FALSE)
write.csv(visit_data_knowledge_post_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/visit_data_knowledge_post_sub.csv", row.names = FALSE)
write.csv(all_individuals_knowledge_sub2,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_knowledge_sub2.csv", row.names = FALSE)
write.csv(dual_events_knowledge_dyads_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/dual_events_knowledge_dyads_sub.csv", row.names = FALSE)
write.csv(dyads_prepost_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/dyads_prepost_sub.csv", row.names = FALSE)

all_individuals_knowledge_sub4_long <- read.csv(all_individuals_knowledge_sub4_long,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_knowledge_sub4_long.csv", row.names = FALSE)
visit_data_knowledge_post_sub <- read.csv(visit_data_knowledge_post_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/visit_data_knowledge_post_sub.csv", row.names = FALSE)
all_individuals_knowledge_sub2 <- read.csv(all_individuals_knowledge_sub2,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_knowledge_sub2.csv", row.names = FALSE)
dual_events_knowledge_dyads_sub <- read.csv(dual_events_knowledge_dyads_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/dual_events_knowledge_dyads_sub.csv", row.names = FALSE)
dyads_prepost_sub <- read.csv(dyads_prepost_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/dyads_prepost_sub.csv", row.names = FALSE)


# (3) SOCIAL NETWORK ----

#Pre-knowledge network: 10 days until knowledge manipulation started
#visit_data_preknowledge <- subset(subset(visit_data, visit_data$day > 169 & visit_data$day < 180))

#Pre-knowledge network: 30 days until knowledge manipulation started
#visit_data_preknowledge <- subset(subset(visit_data, visit_data$day > 149 & visit_data$day < 180))

#Post-knowledge network: 10 days after knowledge manipulation ended
#visit_data_postknowledge <- subset(subset(visit_data, visit_data$day > 216 & visit_data$day < 227))

#Post-knowledge network: 30 days including 20 days of knowledge manipulation and 10 days after
#visit_data_postknowledge <- subset(subset(visit_data, visit_data$day > 196 & visit_data$day < 227))

#Last 63 days of study 
visit_data_prepostknowledge <- subset(subset(visit_data, visit_data$day > 164))
median(visit_data_prepostknowledge$day) #median 178

#Pre-knowledge network: days until knowledge manipulation started
visit_data_preknowledge <- subset(subset(visit_data, visit_data$day > 169 & visit_data$day < 180))
visit_data_preknowledge <- subset(subset(visit_data, visit_data$day > 165 & visit_data$day < 180))

#Post-knowledge network: knowledge manipulation and days after
visit_data_postknowledge <- subset(subset(visit_data, visit_data$day > 179 & visit_data$day < 230))

#Remove JID = NA
visit_data_preknowledge <- subset(visit_data_preknowledge, !is.na(visit_data_preknowledge$JID))
visit_data_postknowledge <- subset(visit_data_postknowledge, !is.na(visit_data_postknowledge$JID))

visit_number_preknowledge <- as.data.frame(table(visit_data_preknowledge$JID))
visit_number_preknowledge$JID <- visit_number_preknowledge$Var1
visit_number_preknowledge$visit_number <- visit_number_preknowledge$Freq
visit_number_preknowledge <- subset(visit_number_preknowledge, select = -c(Var1, Freq))

visit_number_postknowledge <- as.data.frame(table(visit_data_postknowledge$JID))
visit_number_postknowledge$JID <- visit_number_postknowledge$Var1
visit_number_postknowledge$visit_number <- visit_number_postknowledge$Freq
visit_number_postknowledge <- subset(visit_number_postknowledge, select = -c(Var1, Freq))

#Add knowledge status
visit_data_prepostknowledge$knowledgeable <- all_individuals_knowledgeable$knowledge[match(visit_data_prepostknowledge$JID, all_individuals_knowledgeable$JID)]
visit_data_prepostknowledge$knowledge <- ifelse(is.na(visit_data_prepostknowledge$knowledgeable), "naive", "knowledgeable")

#Add visit number at passive feeders
all_individuals_knowledge_sub$visit_number_passive_pre <- visit_number_preknowledge$visit_number[match(all_individuals_knowledge_sub$JID, visit_number_preknowledge$JID)]
all_individuals_knowledge_sub$visit_number_passive_post <- visit_number_postknowledge$visit_number[match(all_individuals_knowledge_sub$JID, visit_number_postknowledge$JID)]
all_individuals_knowledge_sub$visit_number_passive_diff_abs <- all_individuals_knowledge_sub$visit_number_passive_post - all_individuals_knowledge_sub$visit_number_passive_pre

#Add visit number at HK feeders 
all_individuals_knowledge_sub$visit_number_HK_pre <- NA
all_individuals_knowledge_sub$visit_number_HK_post <- visit_number_knowledge_HK$visit_number[match(all_individuals_knowledge_sub$JID, visit_number_knowledge_HK$JID)]

#Add visit number at LK feeders
all_individuals_knowledge_sub$visit_number_LK_pre <- visit_number_knowledge_LK_pre$visit_number[match(all_individuals_knowledge_sub$JID, visit_number_knowledge_LK_pre$JID)]
all_individuals_knowledge_sub$visit_number_LK_post <- visit_number_knowledge_LK_post$visit_number[match(all_individuals_knowledge_sub$JID, visit_number_knowledge_LK_post$JID)]

#Add visit number HK + LK feeders post 
all_individuals_knowledge_sub$visit_number_knowledge_post <- all_individuals_knowledge_sub$visit_number_LK_post + all_individuals_knowledge_sub$visit_number_HK_post

#Absolute change in degree and strength
all_individuals_knowledge_sub$degree_diff_abs <- all_individuals_knowledge_sub$degree_post - all_individuals_knowledge_sub$degree_pre
all_individuals_knowledge_sub$strength_diff_abs <- all_individuals_knowledge_sub$strength_post - all_individuals_knowledge_sub$strength_pre

#Residual change in degree and strength
all_individuals_knowledge_sub$degree_diff_resid <- glm(all_individuals_knowledge_sub$degree_post ~ all_individuals_knowledge_sub$degree_pre, family = negative.binomial(theta = 1))
all_individuals_knowledge_sub$strength_diff_resid <- all_individuals_knowledge_sub$strength_post - all_individuals_knowledge_sub$strength_pre

#Subset based on visit numbers at passive feeders
#all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, !is.na(all_individuals_knowledge_sub$visit_number_preknowledge))
#all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, !is.na(all_individuals_knowledge_sub$visit_number_postknowledge))

#all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, all_individuals_knowledge_sub$visit_number_preknowledge > 4)
#all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, all_individuals_knowledge_sub$visit_number_postknowledge > 4)

all_individuals_knowledge_sub2 <- all_individuals_knowledge_sub
all_individuals_knowledge_sub2 <- all_individuals_knowledge_sub2[, c("JID", "knowledge")]
all_individuals_knowledge_sub3 <- all_individuals_knowledge_sub2
all_individuals_knowledge_sub2$phase <- "pre_knowledge"
all_individuals_knowledge_sub3$phase <- "post_knowledge"

all_individuals_knowledge_sub4 <- rbind(all_individuals_knowledge_sub2, all_individuals_knowledge_sub3)
all_individuals_knowledge_sub4$JID_phase <- paste(all_individuals_knowledge_sub4$JID, all_individuals_knowledge_sub4$phase, sep = " ") 

all_individuals_knowledge_sub5 <- pivot_longer(all_individuals_knowledge_sub[, c("JID", "visit_number_passive_pre", "visit_number_passive_post")], cols = c("visit_number_passive_pre", "visit_number_passive_post"), names_to = "phase_knowledge", values_to = "visit_number_passive")
all_individuals_knowledge_sub5$phase_knowledge <- ifelse(all_individuals_knowledge_sub5$phase_knowledge == "visit_number_passive_pre", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub5$JID_phase <- paste(all_individuals_knowledge_sub5$JID, all_individuals_knowledge_sub5$phase_knowledge, sep = " ") 

all_individuals_knowledge_sub6 <- pivot_longer(all_individuals_knowledge_sub[, c("JID", "visit_number_HK_pre", "visit_number_HK_post")], cols = c("visit_number_HK_pre", "visit_number_HK_post"), names_to = "phase_knowledge", values_to = "visit_number_HK")
all_individuals_knowledge_sub6$phase_knowledge <- ifelse(all_individuals_knowledge_sub6$phase_knowledge == "visit_number_HK_pre", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub6$JID_phase <- paste(all_individuals_knowledge_sub6$JID, all_individuals_knowledge_sub6$phase_knowledge, sep = " ") 

all_individuals_knowledge_sub7 <- pivot_longer(all_individuals_knowledge_sub[, c("JID", "visit_number_LK_pre", "visit_number_LK_post")], cols = c("visit_number_LK_pre", "visit_number_LK_post"), names_to = "phase_knowledge", values_to = "visit_number_LK")
all_individuals_knowledge_sub7$phase_knowledge <- ifelse(all_individuals_knowledge_sub7$phase_knowledge == "visit_number_LK_pre", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub7$JID_phase <- paste(all_individuals_knowledge_sub7$JID, all_individuals_knowledge_sub7$phase_knowledge, sep = " ") 

all_individuals_knowledge_sub8 <- pivot_longer(all_individuals_knowledge_sub[, c("JID", "degree_pre", "degree_post")], cols = c("degree_pre", "degree_post"), names_to = "phase_knowledge", values_to = "degree")
all_individuals_knowledge_sub8$phase_knowledge <- ifelse(all_individuals_knowledge_sub8$phase_knowledge == "degree_pre", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub8$JID_phase <- paste(all_individuals_knowledge_sub8$JID, all_individuals_knowledge_sub8$phase_knowledge, sep = " ") 

all_individuals_knowledge_sub9 <- pivot_longer(all_individuals_knowledge_sub[, c("JID", "strength_pre", "strength_post")], cols = c("strength_pre", "strength_post"), names_to = "phase_knowledge", values_to = "strength")
all_individuals_knowledge_sub9$phase_knowledge <- ifelse(all_individuals_knowledge_sub9$phase_knowledge == "strength_pre", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub9$JID_phase <- paste(all_individuals_knowledge_sub9$JID, all_individuals_knowledge_sub9$phase_knowledge, sep = " ") 

all_individuals_knowledge_sub4$visit_number_passive <- all_individuals_knowledge_sub5$visit_number_passive[match(all_individuals_knowledge_sub4$JID_phase, all_individuals_knowledge_sub5$JID_phase)]
all_individuals_knowledge_sub4$visit_number_HK <- all_individuals_knowledge_sub6$visit_number_HK[match(all_individuals_knowledge_sub4$JID_phase, all_individuals_knowledge_sub6$JID_phase)]
all_individuals_knowledge_sub4$visit_number_LK <- all_individuals_knowledge_sub7$visit_number_LK[match(all_individuals_knowledge_sub4$JID_phase, all_individuals_knowledge_sub7$JID_phase)]
all_individuals_knowledge_sub4$degree <- all_individuals_knowledge_sub8$degree[match(all_individuals_knowledge_sub4$JID_phase, all_individuals_knowledge_sub8$JID_phase)]
all_individuals_knowledge_sub4$strength <- all_individuals_knowledge_sub9$strength[match(all_individuals_knowledge_sub4$JID_phase, all_individuals_knowledge_sub9$JID_phase)]
all_individuals_knowledge_sub4$mean_strength <- all_individuals_knowledge_sub4$strength/all_individuals_knowledge_sub4$degree

#all_individuals_knowledge_sub5 <- pivot_longer(all_individuals_knowledge_sub, cols = c("degree_preknowledge", "degree_postknowledge"), names_to = "phase_knowledge", values_to = "degree")
#all_individuals_knowledge_sub5$phase_knowledge <- ifelse(all_individuals_knowledge_sub5$phase_knowledge == "degree_preknowledge", "pre_knowledge", "post_knowledge")
#all_individuals_knowledge_sub5$JID_phase <- paste(all_individuals_knowledge_sub5$JID, all_individuals_knowledge_sub4$knowledge_phase, sep = " ") 
#all_individuals_knowledge_sub4$degree <- all_individuals_knowledge_sub5$degree[match(all_individuals_knowledge_sub5$JID_phase, all_individuals_knowledge_sub4$JID_phase)]
#all_individuals_knowledge_sub5 <- pivot_longer(all_individuals_knowledge_sub, cols = c("strength_preknowledge", "strength_postknowledge"), names_to = "phase_knowledge", values_to = "strength")
#all_individuals_knowledge_sub5$phase_knowledge <- ifelse(all_individuals_knowledge_sub5$phase_knowledge == "strength_preknowledge", "pre_knowledge", "post_knowledge")
#all_individuals_knowledge_sub5$JID_phase <- paste(all_individuals_knowledge_sub5$JID, all_individuals_knowledge_sub4$knowledge_phase, sep = " ") 
#all_individuals_knowledge_sub4$strength <- all_individuals_knowledge_sub5$strength[match(all_individuals_knowledge_sub5$JID_phase, all_individuals_knowledge_sub4$JID_phase)]

all_individuals_knowledge_sub2 <- all_individuals_knowledge_sub4
all_individuals_knowledge_sub2$knowledge_phase <- paste(all_individuals_knowledge_sub2$knowledge, all_individuals_knowledge_sub2$phase, sep = "_")
all_individuals_knowledge_sub2 <- subset(all_individuals_knowledge_sub2, !is.na(all_individuals_knowledge_sub2$degree))
all_individuals_knowledge_sub2 <- all_individuals_knowledge_sub2 %>% group_by(JID) %>% filter(n()>1)

all_individuals_knowledge_sub2_long <- pivot_longer(all_individuals_knowledge_sub2[, c("JID", "visit_number_LK", "visit_number_HK", "knowledge", "phase")], cols = c("visit_number_LK", "visit_number_HK"), names_to = "feeder", values_to = "visit_number")
all_individuals_knowledge_sub2_long <- subset(all_individuals_knowledge_sub2_long, all_individuals_knowledge_sub2_long$phase == "post_knowledge")

all_individuals_knowledge_sub4_long <- pivot_longer(all_individuals_knowledge_sub4[, c("JID", "visit_number_LK", "visit_number_HK", "knowledge", "phase")], cols = c("visit_number_LK", "visit_number_HK"), names_to = "feeder", values_to = "visit_number")
all_individuals_knowledge_sub4_long <- subset(all_individuals_knowledge_sub4_long, all_individuals_knowledge_sub4_long$phase == "post_knowledge")

all_individuals_knowledge_sub4_long[is.na(all_individuals_knowledge_sub4_long)] <- 0 

rm(all_individuals_knowledge_sub3)
rm(all_individuals_knowledge_sub4) 
rm(all_individuals_knowledge_sub5)
rm(all_individuals_knowledge_sub6)
rm(all_individuals_knowledge_sub7)

all_individuals_knowledge_sub$visit_number_LK <- all_individuals_knowledge_sub2$visit_number_LK[match(all_individuals_knowledge_sub$JID, all_individuals_knowledge_sub2$JID)]
all_individuals_knowledge_sub$visit_number_HK <- all_individuals_knowledge_sub2$visit_number_HK[match(all_individuals_knowledge_sub$JID, all_individuals_knowledge_sub2$JID)]

all_individuals_knowledge_sub_knowledgeable <- subset(all_individuals_knowledge_sub, all_individuals_knowledge_sub$knowledge == "knowledgeable")

#For network before/after knowledge data 
visit_data_RT_preknowledge <- subset(visit_data_RT, visit_data_RT$day > 169 & visit_data_RT$day < 180)
visit_data_RT_postknowledge <- subset(visit_data_RT, visit_data_RT$day > 179)

#Remove data with time slip
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, visit_data_RT_postknowledge$day < 300)

#Remove perches/feeders not considered
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, !feeder == "Y6.1")
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, !feeder == "Y1.2")
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, !feeder == "Y3.2")
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, !feeder == "Z2.2")
visit_data_RT_postknowledge <- subset(visit_data_RT_postknowledge, !feeder == "Z6.1")

visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Y1.2")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Y2.2")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Y3.2")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Z2.2")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Z5.1")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Z5.2")
visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, !feeder == "Z6.2")

visit_data_RT_preknowledge <- subset(visit_data_RT_preknowledge, day > 170)

#Gaussian Mixture Models 

#Prepare variables 
## (i) Time stamp: seconds since start of study period (1 s before first visit)
int_preknowledge <- interval(ymd_hms("2024-06-18 04:00:00 UTC"), ymd_hms(visit_data_RT_preknowledge$time))
visit_data_RT_preknowledge$time_gmm <- time_length(int_preknowledge, "second")

int_postknowledge <- interval(ymd_hms("2024-06-28 04:00:00 UTC"), ymd_hms(visit_data_RT_postknowledge$time))
visit_data_RT_postknowledge$time_gmm <- time_length(int_postknowledge, "second")

#int_knowledge <- interval(ymd_hms("2024-06-18 06:59:00 UTC"), ymd_hms(visit_data_knowledge$start))
#visit_data_knowledge$time <- time_length(int_knowledge, "second")

## (ii) Identity 
visit_data_RT_preknowledge$JID
visit_data_RT_postknowledge$JID
#visit_data_knowledge$JID

#Visit number per individual data set
global_ids_preknowledge <- sort(unique(visit_data_RT_preknowledge$JID))
global_ids_postknowledge <-  sort(unique(visit_data_RT_postknowledge$JID))
#global_ids_knowledge <- visit_number_knowledge$JID

## (iii) Location 
visit_data_RT_preknowledge$position  
visit_data_RT_postknowledge$position  
#visit_data_knowledge$position  

## (iv) Location in time
visit_data_RT_preknowledge$loc_date <-
  paste(visit_data_RT_preknowledge$position,
        visit_data_RT_preknowledge$day,sep="_")

visit_data_RT_postknowledge$loc_date <-
  paste(visit_data_RT_postknowledge$position,
        visit_data_RT_postknowledge$day,sep="_")

#visit_data_knowledge$loc_date <-
#  paste(visit_data_knowledge$position,
#        visit_data_knowledge$day,sep="_")

## Data set 
visit_data_RT_preknowledge_gmm <- visit_data_RT_preknowledge[, c("time_gmm", "JID", "position", "loc_date")]
visit_data_RT_preknowledge_gmm <- na.omit(visit_data_RT_preknowledge_gmm)

visit_data_RT_postknowledge_gmm <- visit_data_RT_postknowledge[, c("time_gmm", "JID", "position", "loc_date")]
visit_data_RT_postknowledge_gmm <- na.omit(visit_data_RT_postknowledge_gmm)

#visit_data_knowledge_gmm <- visit_data_knowledge[, c("time", "JID", "position", "loc_date")]
#visit_data_knowledge_gmm <- na.omit(visit_data_knowledge_gmm)

write.csv(visit_data_RT_preknowledge_gmm,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/visit_data_RT_preknowledge_gmm.csv", row.names = FALSE)
write.csv(visit_data_RT_postknowledge_gmm,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/visit_data_RT_postknowledge_gmm.csv", row.names = FALSE)

write.csv(global_ids_preknowledge,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/global_ids_preknowledge.csv", row.names = FALSE)
write.csv(global_ids_postknowledge,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/global_ids_postknowledge.csv", row.names = FALSE)

# Generate GMM data
gmm_data_preknowledge <- gmmevents(time= visit_data_RT_preknowledge_gmm$time_gmm,
                               identity=visit_data_RT_preknowledge_gmm$JID,
                               location=visit_data_RT_preknowledge_gmm$loc_date,
                               global_ids=global_ids_preknowledge)

gmm_data_postknowledge <- gmmevents(time= visit_data_RT_postknowledge_gmm$time_gmm,
                                   identity=visit_data_RT_postknowledge_gmm$JID,
                                   location=visit_data_RT_postknowledge_gmm$loc_date,
                                   global_ids=global_ids_postknowledge)

#gmm_data_knowledge <- gmmevents(time= visit_data_knowledge_gmm$time,
#                                    identity=visit_data_knowledge_gmm$JID,
#                                    location=visit_data_knowledge_gmm$loc_date,
#                                    global_ids=global_ids_knowledge)


# Extract output
gbi_preknowledge <- gmm_data_preknowledge$gbi
events_preknowledge <- gmm_data_preknowledge$metadata
observations_per_event_preknowledge <- gmm_data_preknowledge$B

gbi_postknowledge <- gmm_data_postknowledge$gbi
events_postknowledge <- gmm_data_postknowledge$metadata
observations_per_event_postknowledge <- gmm_data_postknowledge$B

#gbi_knowledge <- gmm_data_knowledge$gbi
#events_knowledge <- gmm_data_knowledge$metadata
#observations_per_event_knowledge <- gmm_data_knowledge$B

# Can also subset gbi to only individuals observed
# in the dataset to give same answer as if
# global_ids had not been provided
gbi <- gbi[,which(colSums(gbi)>0)]

gbi_preknowledge <- gbi_preknowledge[,which(colSums(gbi_preknowledge)>0)]

# Split up location and date data
tmp_preknowledge <- strsplit(events_preknowledge$Location,"_")
tmp_preknowledge <- do.call("rbind",tmp_preknowledge)
events_preknowledge$Location <- tmp_preknowledge[,1]
events_preknowledge$Date <- tmp_preknowledge[,2]

tmp_postknowledge <- strsplit(events_postknowledge$Location,"_")
tmp_postknowledge <- do.call("rbind",tmp_postknowledge)
events_postknowledge$Location <- tmp_postknowledge[,1]
events_postknowledge$Date <- tmp_postknowledge[,2]

#tmp_knowledge <- strsplit(events_knowledge$Location,"_")
#tmp_knowledge <- do.call("rbind",tmp_knowledge)
#events_knowledge$Location <- tmp_knowledge[,1]
#events_knowledge$Date <- tmp_knowledge[,2]

#Get the adjacency matrix from group by individual matrix 
am_preknowledge <- get_network(gbi_preknowledge, data_format = "GBI",
                  association_index = "SRI", identities = NULL,
                  which_identities = NULL, times = NULL, occurrences = NULL,
                  locations = NULL, which_locations = NULL, start_time = NULL,
                  end_time = NULL, classes = NULL, which_classes = NULL,
                  enter_time = NULL, exit_time = NULL)

am_postknowledge <- get_network(gbi_postknowledge, data_format = "GBI",
                               association_index = "SRI", identities = NULL,
                               which_identities = NULL, times = NULL, occurrences = NULL,
                               locations = NULL, which_locations = NULL, start_time = NULL,
                               end_time = NULL, classes = NULL, which_classes = NULL,
                               enter_time = NULL, exit_time = NULL)

#am_knowledge <- get_network(gbi_knowledge, data_format = "GBI",
#                                association_index = "SRI", identities = NULL,
#                                which_identities = NULL, times = NULL, occurrences = NULL,
#                                locations = NULL, which_locations = NULL, start_time = NULL,
#                                end_time = NULL, classes = NULL, which_classes = NULL,
#                                enter_time = NULL, exit_time = NULL)

write.csv(am,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/am.csv", row.names = FALSE)

am_preknowledge <- as.matrix(read.csv("Data/am_preknowledge.csv", header = T, stringsAsFactors = F))
am_postknowledge <- as.matrix(read.csv("Data/am_postknowledge.csv", header = T, stringsAsFactors = F))

#Get graph from the adjacency matrix
g_preknowledge <- graph_from_adjacency_matrix(am_preknowledge, mode= "undirected",weighted=TRUE,diag=FALSE)

g_postknowledge <- graph_from_adjacency_matrix(am_postknowledge, mode= "undirected",weighted=TRUE,diag=FALSE)

g_knowledge <- graph_from_adjacency_matrix(am_postknowledge, mode= "undirected",weighted=TRUE,diag=FALSE)

#Get sub-graph
edge_weight <- E(g)$weight
mean(E(g)$weight) * 2
summary(edge_weight)
hist(E(g)$weight, breaks = 100)
quantile(edge_weight, probs = seq(0, 1, 1/100))
gs <- subgraph.edges(g, E(g)[E(g)$weight > 0.028571429], del=F)
gs <- subgraph.edges(g, E(g)[E(g)$weight > 0.049019608], del=F)
gs2 <- delete_vertices(gs, degree(gs)==0)

#Vertex attributes
g <- set_vertex_attr(g, "pref_pos", value = perindivperpos2$pref_pos)
gs2 <- set_vertex_attr(gs2, "pref_pos", value = perindivperpos2$pref_pos)

g <- set_vertex_attr(g, "patch_binary", value = all_individuals_prepatch_sub$patch_binary)
gs <- set_vertex_attr(gs, "patch_binary", value = all_individuals_prepatch_sub$patch_binary)

g <- set_vertex_attr(g, "patch_number", value = all_individuals_prepatch_sub$patch_number)
gs <- set_vertex_attr(gs, "patch_number", value = all_individuals_prepatch_sub$patch_number)

g <- set_vertex_attr(g, "mean_OA", value = all_individuals_prepatch_sub$mean_OA_noNA)
gs <- set_vertex_attr(gs, "mean_OA", value = all_individuals_prepatch_sub$mean_OA_noNA)

g <- set_vertex_attr(g, "strength", value = all_individuals_prepatch_sub$strength)
gs <- set_vertex_attr(gs, "strength", value = all_individuals_prepatch_sub$strength)

g <- set_vertex_attr(g, "degree", value = all_individuals_prepatch_sub$degree)
gs <- set_vertex_attr(gs, "degree", value = all_individuals_prepatch_sub$degree)

V(g) <- sort(g)

V(g)$colour <- ifelse(V(g)$pref_pos == "Y1", "black", ifelse(V(g)$pref_pos == "Y2", "red", ifelse(V(g)$pref_pos == "Y3", "yellow", ifelse(V(g)$pref_pos == "Z2", "green", ifelse(V(g)$pref_pos == "Z5", "black", "white")))))
V(gs2)$colour <- ifelse(V(gs2)$pref_pos == "Y1", "blue", ifelse(V(gs2)$pref_pos == "Y2", "red", ifelse(V(gs2)$pref_pos == "Y3", "yellow", ifelse(V(gs2)$pref_pos == "Z2", "green", ifelse(V(gs2)$pref_pos == "Z5", "black", "white")))))

V(g)$colour <- ifelse(V(g)$patch_binary == 1, "black", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "black", "white")

V(g)$colour <- ifelse(V(g)$patch_binary == 1, "grey", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "grey", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "darkgoldenrod1", "white")

coords <- layout_(g, nicely())

plot(g, layout = coords, vertex.size =3, vertex.label = NA, vertex.color = "darkseagreen", edge.width = E(g)$weight * 15, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = "darkseagreen", edge.width = E(g)$weight * 15, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = "white", edge.width = E(g)$weight * 10, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

plot(g, vertex.size = 4, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

plot(gs, vertex.size = 4, vertex.label = NA, vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

V(gs)$label.cex = 0.5
plot(gs, vertex.size = 6, vertex.label = V(gs)$degree, vertex.label.color = "black", vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

V(gs)$label.cex = 0.75
plot(gs, vertex.size = V(gs)$patch_number + 5, vertex.label = V(gs)$patch_number, vertex.label.color = "black", vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

V(gs)$label.cex = 0.75
plot(gs, vertex.size = (5 + V(gs)$strength * 2), vertex.label = V(gs)$mean_OA, vertex.label.color = "black", vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

V(gs)$label.cex = 0.75
plot(gs, vertex.size = (5 + V(gs)$strength * 2), vertex.label = V(gs)$patch_number, vertex.label.color = "black", vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

plot(gs2, vertex.size = 4, vertex.label = NA, vertex.color = "darkseagreen", edge.width = E(gs2)$weight * 15, edge.curved = 0.35)
plot(gs2, vertex.size = 4, vertex.label = NA, vertex.color = V(gs2)$colour, edge.width = E(gs2)$weight * 15, edge.curved = 0.35)
plot(gs2, vertex.size = 4, vertex.label = NA, vertex.color = "white", edge.color = "black", edge.width = E(gs2)$weight * 10, edge.curved = 0.35)

degree_preknowledge <- as.data.frame(degree(g_preknowledge))
degree_preknowledge$JID <- rownames(degree_preknowledge)
degree_preknowledge$degree_preknowledge <- degree_preknowledge$`degree(g_preknowledge)`

strength_preknowledge <- as.data.frame(strength(g_preknowledge))
strength_preknowledge$JID <- rownames(strength_preknowledge)
strength_preknowledge$strength_preknowledge <- strength_preknowledge$`strength(g_preknowledge)`

degree_postknowledge <- as.data.frame(degree(g_postknowledge))
degree_postknowledge$JID <- rownames(degree_postknowledge)
degree_postknowledge$degree_postknowledge <- degree_postknowledge$`degree(g_postknowledge)`

strength_postknowledge <- as.data.frame(strength(g_postknowledge))
strength_postknowledge$JID <- rownames(strength_postknowledge)
strength_postknowledge$strength_postknowledge  <- strength_postknowledge$`strength(g_postknowledge)`

all_individuals_knowledge_sub$degree_pre <- degree_preknowledge$degree_preknowledge[match(all_individuals_knowledge_sub$JID, degree_preknowledge$JID)]
all_individuals_knowledge_sub$strength_pre <- strength_preknowledge$strength_preknowledge[match(all_individuals_knowledge_sub$JID, degree_preknowledge$JID)]

all_individuals_knowledge_sub$degree_post <- degree_postknowledge$degree_postknowledge[match(all_individuals_knowledge_sub$JID, degree_postknowledge$JID)]
all_individuals_knowledge_sub$strength_post <- strength_postknowledge$strength_postknowledge[match(all_individuals_knowledge_sub$JID, degree_postknowledge$JID)]

all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, !is.na(all_individuals_knowledge_sub$degree_preknowledge))
all_individuals_knowledge_sub <- subset(all_individuals_knowledge_sub, !is.na(all_individuals_knowledge_sub$degree_postknowledge))


# (4) DUAL EVENTS ----

#Remove data where no time information 
visit_data_knowledge2 <- visit_data_knowledge[!is.na(visit_data_knowledge$start),]

#Calculate visit interval
visit_data_knowledge2$interval <- interval(visit_data_knowledge2$start,visit_data_knowledge2$end)


#Identify dyadic events ----

## 0. Load packages
library(data.table)
library(lubridate)

## 1. Prepare data

#Calculate visit interval
visit_data_knowledge2$interval <- interval(visit_data_knowledge2$start,visit_data_knowledge2$end)

#If all feeders are used
#visit_data

setDT(visit_data_knowledge2)

#Ensure timestamps are POSIXct with milliseconds
visit_data_knowledge2[, start := ymd_hms(start)]
visit_data_knowledge2[, end   := ymd_hms(end)]
options(digits.secs = 3)

#Assign a unique event ID per visit
visit_data_knowledge2[, event_id := paste(feeder, format(start, "%Y-%m-%d %H:%M:%OS3"), JID, sep = " ")]

#Remove invalid intervals
visit_data_knowledge2 <- visit_data_knowledge2[end >= start]

#Sort for reproducibility
setorder(visit_data_knowledge2, position, start, end)

# 2. Split into primary vs secondary perches
visit_data_knowledge2_HK   <- visit_data_knowledge2[quality == "high"]
visit_data_knowledge2_LK <- visit_data_knowledge2[quality == "low"]  

# 3. Set keys for foverlaps
setkey(visit_data_knowledge2_HK, position, start, end)
setkey(visit_data_knowledge2_LK, position, start, end)

# 4. Overlap join: secondary onto primary
hits <- foverlaps(visit_data_knowledge2_LK, visit_data_knowledge2_HK, type="any", nomatch=0L)

# 5. Compute initiator / joiner and overlap
dual_events <- hits[tag != i.tag]

dual_events[, `:=`(
  initiator = fifelse(start <= i.start, tag, i.tag),
  joiner    = fifelse(start <= i.start, i.tag, tag),
  
  initiator_bout_start = fifelse(start <= i.start, start, i.start),
  initiator_bout_end   = fifelse(start <= i.start, end, i.end),
  
  joiner_bout_start = fifelse(start <= i.start, i.start, start),
  joiner_bout_end   = fifelse(start <= i.start, i.end, end),
  
  overlap_start = pmax(start, i.start),
  overlap_end   = pmin(end, i.end)
)]

# Compute duration in seconds
dual_events[, duration := as.numeric(overlap_end - overlap_start, units = "secs")]

#Add a small constant to all overlaps (optional)
#constant_add <- 0.5  # seconds
#dual_events[, duration := duration + constant_add]
#dual_events[, overlap_end := overlap_start + duration]

dual_events <- dual_events[duration >= 0]

# 6. Assign primary perch event ID
# Primary perch is always the second table in the join (i)
dual_events[, primary_event_id := event_id]
dual_events[, secondary_event_id := i.event_id]

# 7. Build final dual_events dataset
dual_events <- dual_events[, .(
  initiator,
  joiner,
  position,
  initiator_feeder = ifelse(initiator == tag, feeder, i.feeder),
  joiner_feeder    = ifelse(joiner == tag, feeder, i.feeder),
  start = overlap_start,
  end   = overlap_end,
  duration,
  initiator_bout_start,
  initiator_bout_end,
  joiner_bout_start,
  joiner_bout_end,
  year,
  event_id = primary_event_id,
  secondary_event_id = secondary_event_id
)]


#Create blank list
dual_events_knowledge_list <- list()

#For each visit
for (i in 1:(nrow(visit_data_knowledge2) - 1)) {
  progress(i, max.value = nrow(visit_data_knowledge2))
  
  # Find overlapping visits on the same array (position) but different feeders
  overlaps <- which(
    int_overlaps(
      visit_data_knowledge2$interval[i],
      visit_data_knowledge2$interval[i + 1:nrow(visit_data_knowledge2)]
    ) &
      visit_data_knowledge2$position[i] ==
      visit_data_knowledge2$position[i + 1:nrow(visit_data_knowledge2)] &
      visit_data_knowledge2$feeder[i] !=
      visit_data_knowledge2$feeder[i + 1:nrow(visit_data_knowledge2)]
  )
  
  # If there are any overlapping visits
  if (length(overlaps) > 0) {
    for (j in seq_along(overlaps)) {
      
      idx_other <- i + overlaps[j]
      
      # Prevent double-counting: only record each dyad once
      if (i < idx_other) {
        
        # Determine who arrived first
        if (visit_data_knowledge2$start[i] <= visit_data_knowledge2$start[idx_other]) {
          initiator_idx <- i
          joiner_idx <- idx_other
        } else {
          initiator_idx <- idx_other
          joiner_idx <- i
        }
        
        # Create data frame for the overlapping bout
        bouts <- data.frame(
          initiator = visit_data_knowledge2$tag[initiator_idx],
          joiner = visit_data_knowledge2$tag[joiner_idx],
          position = visit_data_knowledge2$position[i],
          initiator_feeder = visit_data_knowledge2$feeder[initiator_idx],
          joiner_feeder = visit_data_knowledge2$feeder[joiner_idx],
          initiator_bout = visit_data_knowledge2$interval[initiator_idx],
          joiner_bout = visit_data_knowledge2$interval[joiner_idx],
          start = max(
            c(visit_data_knowledge2$start[i],
              visit_data_knowledge2$start[idx_other])
          ),
          end = min(
            c(visit_data_knowledge2$end[i],
              visit_data_knowledge2$end[idx_other])
          ),
          duration = min(
            c(visit_data_knowledge2$end[i],
              visit_data_knowledge2$end[idx_other])
          ) - max(
            c(visit_data_knowledge2$start[i],
              visit_data_knowledge2$start[idx_other])
          )
        )
        
        # Add the bout to the growing list
        dual_events_knowledge_list[[length(dual_events_knowledge_list) + 1]] <- bouts
      }
    }
  }
}


for (i in 1:dim(visit_data_knowledge2)[1]) {
  progress(i, max.value = dim(visit_data_knowledge2)[1])
  
  # Find other same array/different feeder overlapping visits
  overlaps <- which(int_overlaps(visit_data_knowledge2$interval[i], visit_data_knowledge2$interval[i+1:dim(visit_data_knowledge2)[1]]) & visit_data_knowledge2$position[i] == visit_data_knowledge2$position[i+1:dim(visit_data_knowledge2)[1]] & visit_data_knowledge2$feeder[i] != visit_data_knowledge2$feeder[i+1:dim(visit_data_knowledge2)[1]])
  
  if(length(overlaps) > 0) {  for (j in 1:length(overlaps)) {  # If there are any
    
    ## Create data-frame for those overlapping bouts, including the two birds and the duration of overlap
    
    bouts <- data.frame(initiator = visit_data_knowledge2$tag[i], joiner = visit_data_knowledge2$tag[i+overlaps[j]], position = visit_data_knowledge2$position[i], initiator_feeder = visit_data_knowledge2$feeder[i], joiner_feeder = visit_data_knowledge2$feeder[i+overlaps[j]], initiator_bout = visit_data_knowledge2$interval[i], joiner_bout = visit_data_knowledge2$interval[i+overlaps[j]], start = max(c(visit_data_knowledge2$start[i],visit_data_knowledge2$start[i+overlaps[j]])), end = min(c(visit_data_knowledge2$end[i],visit_data_knowledge2$end[i+overlaps[j]])), duration = min(c(visit_data_knowledge2$end[i],visit_data_knowledge2$end[i+overlaps[j]])) - max(c(visit_data_knowledge2$start[i],visit_data_knowledge2$start[i+overlaps[j]])))
  }
    # add those into the list
    dual_events_knowledge_list[[i]] <- bouts
  }
}

# Create data frame of dual events
dual_events_knowledge_df <- do.call(rbind, dual_events_knowledge_list)

#Remove cases in which individual is both initiator and joiner
dual_events_knowledge <- dual_events_knowledge_df %>% 
  as_tibble() %>% 
  mutate(duplicates = if_else(initiator == joiner,TRUE,FALSE)) %>% 
  filter(duplicates == FALSE)

#Final dual events data set
dual_events_knowledge <- dual_events_knowledge[!duplicated(dual_events_knowledge),]

#Add initiator and joiner duration
dual_events_knowledge$initiator_duration <- as.duration(dual_events_knowledge$initiator_bout)
dual_events_knowledge$joiner_duration <- as.duration(dual_events_knowledge$joiner_bout)

#Add difference between initiator and joiner arrival
dual_events_knowledge$initiator_arriv <- substr(dual_events_knowledge$initiator_bout, 1, 23)
dual_events_knowledge$joiner_arriv <- substr(dual_events_knowledge$joiner_bout, 1, 23)
dual_events_knowledge$initiator_arriv <- ymd_hms(dual_events_knowledge$initiator_arriv)
dual_events_knowledge$joiner_arriv <- ymd_hms(dual_events_knowledge$joiner_arriv)
dual_events_knowledge$arriv_diff <- as.duration(dual_events_knowledge$joiner_arriv - dual_events_knowledge$initiator_arriv)

#Add difference between initiator and joiner departure
dual_events_knowledge$initiator_dep <- substr(dual_events_knowledge$initiator_bout, 26, 48)
dual_events_knowledge$joiner_dep <- substr(dual_events_knowledge$joiner_bout, 26, 48)
dual_events_knowledge$initiator_dep <- ymd_hms(dual_events_knowledge$initiator_dep)
dual_events_knowledge$joiner_dep<- ymd_hms(dual_events_knowledge$joiner_dep)
dual_events_knowledge$dep_diff <- as.duration(dual_events_knowledge$joiner_dep - dual_events_knowledge$initiator_dep)

#Add initiator and joiner JID 
dual_events_knowledge$initiator_JID <- LH$ID[match(dual_events_knowledge$initiator,LH$RFID)]
dual_events_knowledge$joiner_JID <- LH$ID[match(dual_events_knowledge$joiner,LH$RFID)]

#Add initiator-joiner ID (considers both variants per dyad ID)
dual_events_knowledge$initiator_joiner_ID <- paste(dual_events_knowledge$initiator_JID, dual_events_knowledge$joiner_JID)

#Add unique dyad ID
dual_events_knowledge$dyad_ID <- apply(dual_events_knowledge[c("initiator_JID", "joiner_JID")], 1, function(x) paste(sort(x), collapse=" "))

#Add feeder quality 
dual_events_knowledge$initiator_feeder_quality <- substr(dual_events_knowledge$initiator_feeder,2,2)
dual_events_knowledge$joiner_feeder_quality <- substr(dual_events_knowledge$joiner_feeder,2,2)

#Add high-low quality ID for dyads 
dual_events_knowledge$HKLK_dyad_ID <- ifelse(dual_events_knowledge$initiator_feeder_quality == "H", paste(paste(dual_events_knowledge$initiator_JID, dual_events_knowledge$joiner_JID)), paste(dual_events_knowledge$joiner_JID, dual_events_knowledge$initiator_JID))

dual_events_knowledge$HK_JID <- ifelse(dual_events_knowledge$initiator_feeder_quality == "H", paste(dual_events_knowledge$initiator_JID), paste(dual_events_knowledge$joiner_JID))
dual_events_knowledge$LK_JID <- ifelse(dual_events_knowledge$initiator_feeder_quality == "L", paste(dual_events_knowledge$initiator_JID), paste(dual_events_knowledge$joiner_JID))

dual_events_knowledge$HK_knowledge <- all_individuals_knowledgeable$knowledge[match(dual_events_knowledge$HK_JID, all_individuals_knowledgeable$JID)]
dual_events_knowledge$HK_knowledge <- ifelse(!is.na(dual_events_knowledge$HK_knowledge), "knowledge", "naive")
dual_events_knowledge$LK_knowledge <- all_individuals_knowledgeable$knowledge[match(dual_events_knowledge$LK_JID, all_individuals_knowledgeable$JID)]
dual_events_knowledge$LK_knowledge <- ifelse(!is.na(dual_events_knowledge$LK_knowledge), "knowledge", "naive")

dual_events_knowledge$HKLK_knowledge_comb <- paste(dual_events_knowledge$HK_knowledge, dual_events_knowledge$LK_knowledge, sep = " ")

#Dataframe with number of dyadic events per dyad 
#Each dyad occurs twice for reversed roles (high vs low quality feeder)

dual_events_knowledge_dyads <- as.data.frame(table(dual_events_knowledge$HKLK_dyad_ID))
dual_events_knowledge_dyads$HKLK_dyad_ID <- dual_events_knowledge_dyads$Var1
dual_events_knowledge_dyads$HKLK_no_dual_events <- dual_events_knowledge_dyads$Freq
dual_events_knowledge_dyads$dyad_ID <- dual_events_knowledge$dyad_ID[match(dual_events_knowledge_dyads$HKLK_dyad_ID, dual_events_knowledge$HKLK_dyad_ID)]

dual_events_knowledge_dyads$HK_knowledge <- dual_events_knowledge$HK_knowledge[match(dual_events_knowledge_dyads$HKLK_dyad_ID, dual_events_knowledge$HKLK_dyad_ID)]
dual_events_knowledge_dyads$LK_knowledge <- dual_events_knowledge$LK_knowledge[match(dual_events_knowledge_dyads$HKLK_dyad_ID, dual_events_knowledge$HKLK_dyad_ID)]

dual_events_knowledge_dyads$HKLK_knowledge_comb <- paste(dual_events_knowledge_dyads$HK_knowledge, dual_events_knowledge_dyads$LK_knowledge, sep = " ")

#Remove instances where knowledgeable at LK and naive at HK
dual_events_knowledge_dyads <- subset(dual_events_knowledge_dyads, !HKLK_knowledge_comb == "naive knowledge")

dual_events_knowledge_dyads <- dual_events_knowledge_dyads %>%
  group_by(dyad_ID) %>%           
  mutate(HKLK_no_dual_events_mean = mean(HKLK_no_dual_events))

dual_events_knowledge_dyads <- dual_events_knowledge_dyads%>%
  group_by(dyad_ID) %>% 
  slice(1:1)

dual_events_knowledge_dyads$edge_weight_post <- dyad_edge_post$g_edge_weights_post[match(dual_events_knowledge_dyads$dyad_ID, dyad_edge_post$g_dyad_ID_post)]
dual_events_knowledge_dyads$edge_weight_pre <- dyad_edge_pre$g_edge_weights_pre[match(dual_events_knowledge_dyads$dyad_ID, dyad_edge_pre$g_dyad_ID_pre)]

dual_events_knowledge_dyads$edge_new <- ifelse(is.na(dual_events_knowledge_dyads$edge_weight_pre) & !is.na(dual_events_knowledge_dyads$edge_weight_post), 1, 0)

dual_events_knowledge_dyads$HKLK_knowledge_comb4 <- ifelse(dual_events_knowledge_dyads$HKLK_knowledge_comb == "naive naive", "naive", "knowledge")


#Full dataset on dyadic associations and interactions

dyads_prepost <- dyad_edge_post
dyads_prepost$dyad_ID <- dyads_prepost$g_dyad_ID_post
dyads_prepost$g_edge_weights_pre <- dyad_edge_pre$g_edge_weights_pre[match(dyads_prepost$g_dyad_ID_post, dyad_edge_pre$g_dyad_ID_pre)]
dyads_prepost$g_edge_weights_pre <- ifelse(!is.na(dyads_prepost$g_edge_weights_pre), dyads_prepost$g_edge_weights_pre, 0)
dyads_prepost$HKLK_no_dual_events_mean <- dual_events_knowledge_dyads$HKLK_no_dual_events_mean[match(dyads_prepost$dyad_ID, dual_events_knowledge_dyads$dyad_ID)]
dyads_prepost$HKLK_no_dual_events_mean <- ifelse(!is.na(dyads_prepost$HKLK_no_dual_events_mean), dyads_prepost$HKLK_no_dual_events_mean, 0)
dyads_prepost$HKLK_no_dual_events_binary <- ifelse(dyads_prepost$HKLK_no_dual_events_mean > 0, 1, 0)
dyads_prepost$edge_diff_abs <- dyads_prepost$g_edge_weights_post - dyads_prepost$g_edge_weights_pre
dyads_prepost$edge_diff_resid <- resid(lm(g_edge_weights_post ~ g_edge_weights_pre, data = dyads_prepost))
dyads_prepost$knowledge_comb3 <- ifelse(dyads_prepost$knowledge_comb2 == "knowledge naive", "different", "same")
dyads_prepost$knowledge_comb4 <- ifelse(dyads_prepost$knowledge_comb2 == "naive naive", "naive", "knowledge")
dyads_prepost$ind1_visit_no_pre <- visit_number_preknowledge$visit_number[match(dyads_prepost$ind1_ID, visit_number_preknowledge$JID)]
dyads_prepost$ind2_visit_no_pre <- visit_number_preknowledge$visit_number[match(dyads_prepost$ind2_ID, visit_number_preknowledge$JID)]
dyads_prepost$ind1_visit_no_post <- visit_number_postknowledge$visit_number[match(dyads_prepost$ind1_ID, visit_number_postknowledge$JID)]
dyads_prepost$ind2_visit_no_post <- visit_number_postknowledge$visit_number[match(dyads_prepost$ind2_ID, visit_number_postknowledge$JID)]
dyads_prepost$new_edge <- ifelse(dyads_prepost$g_edge_weights_pre == 0 & dyads_prepost$g_edge_weights_post >0, "new", "old")

mean(dyads_prepost$edge_diff_abs)
3 * sd(dyads_prepost$edge_diff_abs)

dyads_prepost_sub <- subset(dyads_prepost, !is.na(ind1_visit_no_pre))
dyads_prepost_sub <- subset(dyads_prepost_sub, !is.na(ind2_visit_no_pre))
dyads_prepost_sub <- subset(dyads_prepost_sub, !is.na(ind1_visit_no_post))
dyads_prepost_sub <- subset(dyads_prepost_sub, !is.na(ind2_visit_no_post))

dyads_prepost_sub <- subset(dyads_prepost_sub, edge_diff_abs <  0.1507949)
dyads_prepost_sub <- subset(dyads_prepost_sub, edge_diff_abs >  -0.09390224)
dyads_prepost_sub <- subset(dyads_prepost_sub, !g_edge_weights_pre == 0)

dyads_prepost_sub <- subset(dyads_prepost, new_edge == "old")
dyads_prepost_sub$edge_diff_resid <- resid(lm(g_edge_weights_post ~ g_edge_weights_pre, data = dyads_prepost_sub))

dyads_prepost_sub$knowledge_comb2_HKLK_dual_events_binary <- paste(dyads_prepost_sub$knowledge_comb2, dyads_prepost_sub$HKLK_no_dual_events_binary, sep = "_")

#Dataset focusing on how often knowledgeable individuals were observed performing their "skill"
dual_events_knowledge_knowledgable <- subset(dual_events_knowledge, HK_knowledge == "knowledge")

dual_events_knowledge_knowledgable <- as.data.frame(table(dual_events_knowledge_knowledgable$HK_JID))
dual_events_knowledge_knowledgable$JID <- dual_events_knowledge_knowledgable$Var1
dual_events_knowledge_knowledgable$HKLK_no_dual_events <- dual_events_knowledge_knowledgable$Freq

all_individuals_knowledgeable$HKLK_no_dual_events <- dual_events_knowledge_knowledgable$HKLK_no_dual_events[match(all_individuals_knowledgeable$JID, dual_events_knowledge_knowledgable$JID)]
all_individuals_knowledgeable$visit_number_HK <- all_individuals_knowledge_sub$visit_number_HK_post[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$visit_number_passive_pre <- all_individuals_knowledge_sub$visit_number_passive_pre[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$visit_number_passive_post <- all_individuals_knowledge_sub$visit_number_passive_post[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$degree_pre <- all_individuals_knowledge_sub$degree_pre[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$degree_post <- all_individuals_knowledge_sub$degree_post[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$strength_pre <- all_individuals_knowledge_sub$strength_pre[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]
all_individuals_knowledgeable$strength_post <- all_individuals_knowledge_sub$strength_post[match(all_individuals_knowledgeable$JID, all_individuals_knowledge_sub$JID)]

all_individuals_knowledgeable$degree_diff_abs <- all_individuals_knowledgeable$degree_post - all_individuals_knowledgeable$degree_pre
all_individuals_knowledgeable$strength_diff_abs <- all_individuals_knowledgeable$strength_post - all_individuals_knowledgeable$strength_pre

all_individuals_knowledge_sub$degree_diff_abs <- all_individuals_knowledge_sub$degree_post - all_individuals_knowledge_sub$degree_pre
all_individuals_knowledge_sub$strength_diff_abs <- all_individuals_knowledge_sub$strength_post - all_individuals_knowledge_sub$strength_pre


# (5) DISPLACEMENTS ----

# Sort by feeder and start_time
visit_data_knowledge_displace <- visit_data_knowledge %>% arrange(feeder, start)
visit_data_knowledge_displace_post <- subset(visit_data_knowledge_displace, day > 179)

# Detect potential displacements
displacements_knowledge_post <- visit_data_knowledge_displace_post %>%
  group_by(feeder) %>%
  mutate(next_JID = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_JID, start, time_diff)

displacements_knowledge_post$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_knowledge_post$displaced_JID, all_individuals_knowledgeable$JID)]
displacements_knowledge_post$displaced_knowledge <- ifelse(is.na(displacements_knowledge_post$displaced_knowledge), "naive", "knowledge")

displacements_knowledge_post$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_knowledge_post$displacer_JID, all_individuals_knowledgeable$JID)]
displacements_knowledge_post$displacer_knowledge <- ifelse(is.na(displacements_knowledge_post$displacer_knowledge), "naive", "knowledge")

displacements_knowledge_post$knowledge_comb <- paste(displacements_knowledge_post$displaced_knowledge, displacements_knowledge_post$displacer_knowledge, sep = "_")
table(displacements_knowledge_post$knowledge_comb)

displacements_knowledge_post$feeder_JID_time <- paste(displacements_knowledge_post$feeder, displacements_knowledge_post$displaced_JID, displacements_knowledge_post$start, sep = "_")

displacements_knowledge_post$displaced <- 1

visit_data_knowledge_post$feeder_JID_time <- paste(visit_data_knowledge_post$feeder, visit_data_knowledge_post$JID, visit_data_knowledge_post$start, sep = "_")

visit_data_knowledge_post$displaced <- displacements_knowledge_post$displaced[match(visit_data_knowledge_post$feeder_JID_time, displacements_knowledge_post$feeder_JID_time)]
visit_data_knowledge_post$displaced <- ifelse(is.na(visit_data_knowledge_post$displaced), 0, 1)

table(visit_data_knowledge_post$displaced, visit_data_knowledge_post$knowledge, visit_data_knowledge_post$quality)

visit_data_knowledge_post <- visit_data_knowledge_post %>%
  mutate(next_JID = lead(JID))

visit_data_knowledge_post$next_knowledge <- all_individuals_knowledgeable$knowledge[match(visit_data_knowledge_post$next_JID, all_individuals_knowledgeable$JID)]
visit_data_knowledge_post$next_knowledge <- ifelse(is.na(visit_data_knowledge_post$next_knowledge), "naive", "knowledge")

visit_data_knowledge_post$knowledge_comb <- paste(visit_data_knowledge_post$knowledge, visit_data_knowledge_post$next_knowledge, sep = "_")

visit_data_knowledge_post$visit_number_knowledge_post <- visit_number_knowledge_post$visit_number_knowledge[match(visit_data_knowledge_post$JID, visit_number_knowledge_post$JID)]

visit_data_knowledge_post$dyad_ID <- ifelse(visit_data_knowledge_post$JID < visit_data_knowledge_post$next_JID, paste(visit_data_knowledge_post$JID, visit_data_knowledge_post$next_JID, sep = " "), paste(visit_data_knowledge_post$next_JID, visit_data_knowledge_post$JID, sep =  " "))

visit_data_knowledge_post$edge_weight_pre <- dyad_edge_pre$g_edge_weights_pre[match(visit_data_knowledge_post$dyad_ID, dyad_edge_pre$g_dyad_ID_pre)]

visit_data_knowledge_post$same_diff_indiv <- ifelse(visit_data_knowledge_post$JID == visit_data_knowledge_post$next_JID, "same", "different")

visit_data_knowledge_post_sub <- subset(visit_data_knowledge_post, !same_diff_indiv == "same")
visit_data_knowledge_post_sub <- subset(visit_data_knowledge_post, !is.na(edge_weight_pre))
#visit_data_knowledge_post_sub <- sample_n(visit_data_knowledge_post, 15000) 

visit_data_knowledge_post_sub$age <- Ringed$min_age[match(visit_data_knowledge_post_sub$JID, Ringed$ID)]
visit_data_knowledge_post_sub$next_age <- Ringed$min_age[match(visit_data_knowledge_post_sub$next_JID, Ringed$ID)]
visit_data_knowledge_post_sub$JID[visit_data_knowledge_post_sub$JID == "J55774"] <- "J5774"
visit_data_knowledge_post_sub$next_JID[visit_data_knowledge_post_sub$next_JID == "J55774"] <- "J5774"

visit_data_knowledge_post_sub$age[visit_data_knowledge_post_sub$JID == "J5084"] <- 3
visit_data_knowledge_post_sub$age[visit_data_knowledge_post_sub$JID == "J5657B"] <- 1
visit_data_knowledge_post_sub$age[visit_data_knowledge_post_sub$JID == "J5774"] <- 1

visit_data_knowledge_post_sub$next_age[visit_data_knowledge_post_sub$next_JID == "J5084"] <- 3
visit_data_knowledge_post_sub$next_age[visit_data_knowledge_post_sub$next_JID == "J5657B"] <- 1
visit_data_knowledge_post_sub$next_age[visit_data_knowledge_post_sub$next_JID == "J5774"] <- 1

visit_data_knowledge_post_sub$age_class <- ifelse(visit_data_knowledge_post_sub$age > 1, "adult", "juvenile")
visit_data_knowledge_post_sub$next_age_class <- ifelse(visit_data_knowledge_post_sub$next_age > 1, "adult", "juvenile")
visit_data_knowledge_post_sub$age_class_comb <- paste(visit_data_knowledge_post_sub$age_class, visit_data_knowledge_post_sub$next_age_class, sep = "_")

visit_data_knowledge_post_sub$dyad_ID <- ifelse(visit_data_knowledge_post_sub$JID < visit_data_knowledge_post_sub$next_JID, paste(visit_data_knowledge_post_sub$JID, visit_data_knowledge_post_sub$next_JID, sep = " "), paste(visit_data_knowledge_post_sub$next_JID, visit_data_knowledge_post_sub$JID, sep = " "))

displacements_knowledge_post_summary <- visit_data_knowledge_post %>%
  group_by(JID, knowledge) %>%
  summarise(
    n_visits = n(),
    n_displaced = sum(displaced),
    .groups = "drop"
  )

# Sort by feeder and start_time
visit_data_knowledge_displace <- visit_data_knowledge %>% arrange(feeder, start)

# Detect potential displacements
displacements <- visit_data_knowledge_displace %>%
  group_by(feeder) %>%
  mutate(next_JID = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_JID, start, time_diff)

displacements$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements$displaced_JID, all_individuals_knowledgeable$JID)]
displacements$displaced_knowledge <- ifelse(is.na(displacements$displaced_knowledge), "naive", "knowledge")

displacements$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements$displacer_JID, all_individuals_knowledgeable$JID)]
displacements$displacer_knowledge <- ifelse(is.na(displacements$displacer_knowledge), "naive", "knowledge")

displacements$knowledge_comb <- paste(displacements$displaced_knowledge, displacements$displacer_knowledge, sep = "_")
table(displacements$knowledge_comb)

displacements$feeder_JID_time <- paste(displacements$feeder, displacements$displaced_JID, displacements$start, sep = "_")

displacements$displaced <- 1

visit_data_knowledge$feeder_JID_time <- paste(visit_data_knowledge$feeder, visit_data_knowledge$JID, visit_data_knowledge$start, sep = "_")

visit_data_knowledge$displaced <- displacements$displaced[match(visit_data_knowledge$feeder_JID_time, displacements$feeder_JID_time)]
visit_data_knowledge$displaced <- ifelse(is.na(visit_data_knowledge$displaced), 0, 1)

table(visit_data_knowledge$displaced, visit_data_knowledge$knowledge, visit_data_knowledge$quality)

# Sort by feeder and start_time
visit_data_knowledge_HK_displace <- visit_data_knowledge_HK %>% arrange(feeder, start)

# Detect potential displacements
displacements_HK <- visit_data_knowledge_HK_displace %>%
  group_by(feeder) %>%
  mutate(next_JID = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_JID, start, time_diff)

displacements_HK$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_HK$displaced_JID, all_individuals_knowledgeable$JID)]
displacements_HK$displaced_knowledge <- ifelse(is.na(displacements_HK$displaced_knowledge), "naive", "knowledge")

displacements_HK$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_HK$displacer_JID, all_individuals_knowledgeable$JID)]
displacements_HK$displacer_knowledge <- ifelse(is.na(displacements_HK$displacer_knowledge), "naive", "knowledge")

displacements_HK$knowledge_comb <- paste(displacements_HK$displaced_knowledge, displacements_HK$displacer_knowledge, sep = "_")
table(displacements_HK$knowledge_comb)

displacements_HK$feeder_JID_time <- paste(displacements_HK$feeder, displacements_HK$displaced_JID, displacements_HK$start, sep = "_")

displacements_HK$displaced <- 1

visit_data_knowledge_HK$feeder_JID_time <- paste(visit_data_knowledge_HK$feeder, visit_data_knowledge_HK$JID, visit_data_knowledge_HK$start, sep = "_")

visit_data_knowledge_HK$displaced <- displacements$displaced[match(visit_data_knowledge_HK$feeder_JID_time, displacements_HK$feeder_JID_time)]
visit_data_knowledge_HK$displaced <- ifelse(is.na(visit_data_knowledge_HK$displaced), 0, 1)

table(visit_data_knowledge_HK$displaced, visit_data_knowledge_HK$knowledge)


# Sort by feeder and start_time
visit_data_knowledge_LK_displace <- visit_data_knowledge_LK %>% arrange(feeder, start)

# Detect potential displacements
displacements_LK <- visit_data_knowledge_LK_displace %>%
  group_by(feeder) %>%
  mutate(next_JID = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_JID, start, time_diff)

displacements_LK$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_LK$displaced_JID, all_individuals_knowledgeable$JID)]
displacements_LK$displaced_knowledge <- ifelse(is.na(displacements_LK$displaced_knowledge), "naive", "knowledge")

displacements_LK$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_LK$displacer_JID, all_individuals_knowledgeable$JID)]
displacements_LK$displacer_knowledge <- ifelse(is.na(displacements_LK$displacer_knowledge), "naive", "knowledge")

displacements_LK$knowledge_comb <- paste(displacements_LK$displaced_knowledge, displacements_LK$displacer_knowledge, sep = "_")

table(displacements_LK$knowledge_comb)

displacements_LK$feeder_JID_time <- paste(displacements_LK$feeder, displacements_LK$displaced_JID, displacements_LK$start, sep = "_")

displacements_LK$displaced <- 1

visit_data_knowledge_LK$feeder_JID_time <- paste(visit_data_knowledge_LK$feeder, visit_data_knowledge_LK$JID, visit_data_knowledge_LK$start, sep = "_")

visit_data_knowledge_LK$displaced <- displacements_LK$displaced[match(visit_data_knowledge_LK$feeder_JID_time, displacements_LK$feeder_JID_time)]
visit_data_knowledge_LK$displaced <- ifelse(is.na(visit_data_knowledge_LK$displaced), 0, 1)

visit_data_knowledge_LK$period <- ifelse(visit_data_knowledge_LK$day < 180, "preknowledge", "postknowledge")

table(visit_data_knowledge_LK$displaced, visit_data_knowledge_LK$knowledge, visit_data_knowledge_LK$period)

visit_data_knowledge_LK <- visit_data_knowledge_LK %>%
  mutate(next_JID = lead(JID))

visit_data_knowledge_LK$next_knowledge <- all_individuals_knowledgeable$knowledge[match(visit_data_knowledge_LK$next_JID, all_individuals_knowledgeable$JID)]
visit_data_knowledge_LK$next_knowledge <- ifelse(is.na(visit_data_knowledge_LK$next_knowledge), "naive", "knowledge")

visit_data_knowledge_LK$knowledge_comb <- paste(visit_data_knowledge_LK$knowledge, visit_data_knowledge_LK$next_knowledge, sep = "_")

visit_data_knowledge_LK$dyad_ID <- ifelse(visit_data_knowledge_LK$JID < visit_data_knowledge_LK$next_JID, paste(visit_data_knowledge_LK$JID, visit_data_knowledge_LK$next_JID, sep = " "), paste(visit_data_knowledge_LK$next_JID, visit_data_knowledge_LK$JID, sep = " "))

visit_data_knowledge_LK_sub <- subset(visit_data_knowledge_LK, day < 190)

visit_data_knowledge_LK_sub <- visit_data_knowledge_post_sub$age_class_comb[match()]

visit_data_knowledge_LK_sub$age <- Ringed$min_age[match(visit_data_knowledge_LK_sub$JID, Ringed$ID)]
visit_data_knowledge_LK_sub$age[visit_data_knowledge_LK_sub$JID == "J5084"] <- 3
visit_data_knowledge_LK_sub$age_class <- ifelse(visit_data_knowledge_LK_sub$age > 1, "adult", "juvenile")

visit_data_knowledge_LK_sub$next_age <- Ringed$min_age[match(visit_data_knowledge_LK_sub$next_JID, Ringed$ID)]
visit_data_knowledge_LK_sub$next_age[visit_data_knowledge_LK_sub$next_JID == "J5084"] <- 3
visit_data_knowledge_LK_sub$next_age_class <- ifelse(visit_data_knowledge_LK_sub$next_age > 1, "adult", "juvenile")

visit_data_knowledge_LK_sub$age_class_comb <- paste(visit_data_knowledge_LK_sub$age_class, visit_data_knowledge_LK_sub$next_age_class, sep = "_")

visit_data_knowledge_LK_sub$dyad_ID <- ifelse(visit_data_knowledge_LK_sub$JID < visit_data_knowledge_LK_sub$next_JID, paste(visit_data_knowledge_LK_sub$JID, visit_data_knowledge_LK_sub$next_JID, sep = " "), paste(visit_data_knowledge_LK_sub$next_JID, visit_data_knowledge_LK_sub$JID, sep = " "))

displacement_knowledge_LK_summary <- visit_data_knowledge_LK_pre %>%
  group_by(JID, knowledge) %>%
  summarise(
    n_visits = n(),
    n_displaced = sum(displaced),
    .groups = "drop"
  )

# Sort by feeder and start_time
visit_data_prepostknowledge_displace <- visit_data_prepostknowledge %>% arrange(feeder, start)

# Detect potential displacements
displacements_prepostknowledge <- visit_data_prepostknowledge_displace %>%
  group_by(feeder) %>%
  mutate(next_JID = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_JID, start, time_diff)

displacements_prepostknowledge$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_prepostknowledge$displaced_JID, all_individuals_knowledgeable$JID)]
displacements_prepostknowledge$displaced_knowledge <- ifelse(is.na(displacements_prepostknowledge$displaced_knowledge), "naive", "knowledge")

displacements_prepostknowledge$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_prepostknowledge$displacer_JID, all_individuals_knowledgeable$JID)]
displacements_prepostknowledge$displacer_knowledge <- ifelse(is.na(displacements_prepostknowledge$displacer_knowledge), "naive", "knowledge")

displacements_prepostknowledge$knowledge_comb <- paste(displacements_prepostknowledge$displaced_knowledge, displacements_prepostknowledge$displacer_knowledge, sep = "_")

table(displacements_prepostknowledge$knowledge_comb)

displacements_prepostknowledge$feeder_JID_time <- paste(displacements_prepostknowledge$feeder, displacements_prepostknowledge$displaced_JID, displacements_prepostknowledge$start, sep = "_")

displacements_prepostknowledge$displaced <- 1

visit_data_prepostknowledge$feeder_JID_time <- paste(visit_data_prepostknowledge$feeder, visit_data_prepostknowledge$JID, visit_data_prepostknowledge$start, sep = "_")

visit_data_prepostknowledge$displaced <- displacements_prepostknowledge$displaced[match(visit_data_prepostknowledge$feeder_JID_time, displacements_prepostknowledge$feeder_JID_time)]
visit_data_prepostknowledge$displaced <- ifelse(is.na(visit_data_prepostknowledge$displaced), 0, 1)

visit_data_prepostknowledge$prepostknowledge <- ifelse(visit_data_prepostknowledge$day < 180, "preknowledge", "postknowledge")

table(visit_data_prepostknowledge$displaced, visit_data_prepostknowledge$knowledge, visit_data_prepostknowledge$prepostknowledge)

visit_data_prepostknowledge <- visit_data_prepostknowledge %>%
  mutate(next_JID = lead(JID))

visit_data_prepostknowledge$next_knowledge <- all_individuals_knowledgeable$knowledge[match(visit_data_prepostknowledge$next_JID, all_individuals_knowledgeable$JID)]
visit_data_prepostknowledge$next_knowledge <- ifelse(!is.na(visit_data_prepostknowledge$next_knowledge), "knowledge", "naive")

visit_data_prepostknowledge <- subset(visit_data_prepostknowledge, visit_data_prepostknowledge$next_knowledge == "naive")

displacements_prepostknowledge_summary <- visit_data_prepostknowledge[visit_data_prepostknowledge$prepostknowledge =="postknowledge", ] %>%
  group_by(JID, knowledge) %>%
  summarise(
    n_visits = n(),
    n_displaced = sum(displaced),
    .groups = "drop"
  )


# Sort by feeder and start_time
visit_data_preknowledge_displace <- visit_data_preknowledge %>% arrange(feeder, start)

# Detect potential displacements
displacements_preknowledge <- visit_data_preknowledge_displace %>%
  group_by(feeder) %>%
  mutate(next_bird = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced_JID = JID, displacer_JID = next_bird, start, time_diff)

displacements_preknowledge$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_preknowledge$displaced, all_individuals_knowledgeable$JID)]
displacements_preknowledge$displaced_knowledge <- ifelse(is.na(displacements_preknowledge$displaced_knowledge), "naive", "knowledge")

displacements_preknowledge$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_preknowledge$displacer, all_individuals_knowledgeable$JID)]
displacements_preknowledge$displacer_knowledge <- ifelse(is.na(displacements_preknowledge$displacer_knowledge), "naive", "knowledge")

displacements_preknowledge$knowledge_comb <- paste(displacements_preknowledge$displaced_knowledge, displacements_preknowledge$displacer_knowledge, sep = "_")

table(displacements_preknowledge$knowledge_comb)


# Sort by feeder and start_time
visit_data_postknowledge_displace <- visit_data_postknowledge %>% arrange(feeder, start)

# Detect potential displacements
displacements_postknowledge <- visit_data_postknowledge_displace %>%
  group_by(feeder) %>%
  mutate(next_bird = lead(JID),
         next_start = lead(start),
         time_diff = next_start - end) %>%
  filter(!is.na(time_diff) & time_diff >= 0 & time_diff <= 2) %>%
  dplyr::select(feeder, displaced = JID, displacer = next_bird, start, time_diff)

displacements_postknowledge$displaced_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_postknowledge$displaced, all_individuals_knowledgeable$JID)]
displacements_postknowledge$displaced_knowledge <- ifelse(is.na(displacements_postknowledge$displaced_knowledge), "naive", "knowledge")

displacements_postknowledge$displacer_knowledge <- all_individuals_knowledgeable$knowledge[match(displacements_postknowledge$displacer, all_individuals_knowledgeable$JID)]
displacements_postknowledge$displacer_knowledge <- ifelse(is.na(displacements_postknowledge$displacer_knowledge), "naive", "knowledge")

displacements_postknowledge$knowledge_comb <- paste(displacements_postknowledge$displaced_knowledge, displacements_postknowledge$displacer_knowledge, sep = "_")

table(displacements_postknowledge$knowledge_comb)

table(displacements$displaced, displacements$displaced_knowledge)
unique(sort(c(displacements$displaced, displacements$displacer)))



# (6) FIRST EXPLORATION ----

table(visit_data_knowledge$knowledge, visit_data_knowledge$quality)
table(all_individuals_knowledge_sub$knowledge)

ggplot(visit_data_knowledge, aes(x = day, colour = quality_knowledge)) + geom_point(stat = "count") 
ggplot(visit_data, aes(x = day)) + geom_point(stat = "count") 

mean(all_individuals_knowledge_sub$degree_preknowledge, na.rm = TRUE)
mean(all_individuals_knowledge_sub$degree_postknowledge, na.rm = TRUE)
mean(all_individuals_knowledge_sub$strength_preknowledge, na.rm = TRUE)
mean(all_individuals_knowledge_sub$strength_postknowledge, na.rm = TRUE)

t.test(all_individuals_knowledge_sub$strength_preknowledge, all_individuals_knowledge_sub$strength_postknowledge)

aggregate(all_individuals_knowledge_sub$degree_preknowledge, by = list(all_individuals_knowledge_sub$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub$degree_postknowledge, by = list(all_individuals_knowledge_sub$knowledge), FUN = "mean", na.rm = TRUE)

aggregate(all_individuals_knowledge_sub$strength_preknowledge, by = list(all_individuals_knowledge_sub$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub$strength_postknowledge, by = list(all_individuals_knowledge_sub$knowledge), FUN = "mean", na.rm = TRUE)

mean(visit_data_knowledge_HK$visit_duration)
sd(visit_data_knowledge_HK$visit_duration)

mean(visit_data_knowledge_LK$visit_duration)
sd(visit_data_knowledge_LK$visit_duration)

aggregate(visit_data_knowledge_HK$visit_duration, by = list(visit_data_knowledge_HK$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(visit_data_knowledge_LK$visit_duration, by = list(visit_data_knowledge_LK$knowledge), FUN = "mean", na.rm = TRUE)


#(7) STATISTICAL MODELLING ----

#(7.1) Feeder visits ---- 

#Number of visits at both HK and LK feeders (filtered to include only data after manipulation)

all_individuals_knowledge_sub4_long$obs <- row.names(all_individuals_knowledge_sub4_long)

#Model just with prior
default_prior()

visit_no_brm1.1_prior <- brm(visit_number ~ knowledge * feeder + (1|JID), 
                              data = all_individuals_knowledge_sub4_long, 
                              family = poisson(link = "log"),
                              prior = c(
                                set_prior("normal(log(90), 1)", class = "Intercept"),  
                                set_prior("normal(0, 0.5)", class = "b")),              
                              sample_prior = "only")

visit_no_brm1.2_prior <- brm(visit_number ~ knowledge * feeder + (1|JID) + (1|obs), 
                           data = all_individuals_knowledge_sub4_long, 
                           family = poisson(link = "log"),
                           prior = c(
                             set_prior("normal(log(90), 1)", class = "Intercept"),  
                             set_prior("normal(0, 0.5)", class = "b")),              
                           sample_prior = "only")

visit_no_brm1.3_prior <- brm(visit_number ~ knowledge * feeder + (1|JID), 
                           data = all_individuals_knowledge_sub4_long, 
                           family = negbinomial(link = "log"),
                           prior = c(
                             set_prior("normal(log(90), 1)", class = "Intercept"),  
                             set_prior("normal(0, 0.5)", class = "b")),              
                           sample_prior = "only")

#Prior predictive checks 
pp_check(visit_no_brm1_prior, ndraws = 1000)

#Model

#Poisson
visit_no_brm1.1 <- brm(visit_number ~ knowledge * feeder + (1|JID), 
                     data = all_individuals_knowledge_sub4_long, 
                     family = poisson(link = "log"),
                     prior = c(
                       set_prior("normal(log(90), 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),              
                     save_pars = save_pars(all = TRUE)
)

#Poisson with observation-level random effect
visit_no_brm1.2 <- brm(visit_number ~ knowledge * feeder + (1|JID) + (1|obs), 
                       data = all_individuals_knowledge_sub4_long, 
                       family = poisson(link = "log"),
                       prior = c(
                         set_prior("normal(log(90), 1)", class = "Intercept"),  
                         set_prior("normal(0, 0.5)", class = "b")),              
                       save_pars = save_pars(all = TRUE)
)

#Negative binomial 
visit_no_brm1.3 <- brm(visit_number ~ knowledge * feeder + (1|JID), 
                        data = all_individuals_knowledge_sub4_long, 
                        family = negbinomial(link = "log"),
                        prior = c(
                          set_prior("normal(log(90), 1)", class = "Intercept"),  
                          set_prior("normal(0, 0.5)", class = "b")),
                       save_pars = save_pars(all = TRUE)
                        
)

visit_no_brm1.1_loo <- loo(visit_no_brm1.1, moment_match = TRUE)
visit_no_brm1.2_loo <- loo(visit_no_brm1.2, moment_match = TRUE)
visit_no_brm1.3_loo <- loo(visit_no_brm1.3, moment_match = TRUE)

loo_compare(visit_no_brm1.2_loo, visit_no_brm1.3_loo)

#Poisson OLRE-model performs best

visit_no_brm1 <- visit_no_brm1.2

summary(visit_no_brm1)

check_collinearity(visit_no_brm1)

#Posterior predictive checks
pp_check(visit_no_brm1, type = "dens_overlay", ndraws = 100) 
pp_check(visit_no_brm1, type = "dens_overlay", ndraws = 100) + scale_x_log10()

pp_check(visit_no_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(visit_no_brm1, type = "bars", ndraws = 1000)
pp_check(visit_no_brm1, type = "hist", ndraws = 1000)

y <- model.frame(visit_no_brm1)[[1]]
yrep <- posterior_predict(visit_no_brm1)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(visit_no_brm1)
launch_shinystan(visit_no_brm1)

#Posterior distribution
as_draws_df(visit_no_brm1)
mcmc_areas(visit_no_brm1)
mcmc_intervals(visit_no_brm1)

#Evaluation and interpretation
fitted(visit_no_brm1, scale = "response")
conditional_effects(visit_no_brm1)
conditional_effects(visit_no_brm1,effects = "knowledge:feeder")
bayes_R2(visit_no_brm1)
emmeans(visit_no_brm1, ~ knowledge * feeder, type = "response") |> pairs()
hypothesis(visit_no_brm1, "knowledgenaive > 0") 
hypothesis(visit_no_brm1, "feedervisit_number_LK > 0") 

#Sensitivity analysis and prior checks 
prior_summary(visit_no_brm1)

#Extract predictions
fitted(visit_no_brm1)

newdata <- expand.grid(knowledge = c("knowledgeable","naive"), feeder = c("visit_number_LK","visit_number_HK"))
preds <- posterior_epred(visit_no_brm1, newdata = newdata)
apply(preds, 2, mean)   # posterior mean per combination
apply(preds, 2, quantile, c(0.025, 0.975))

LK_visit_boxplot <- ggplot(all_individuals_knowledge_sub2, aes(x = knowledge, y = visit_number_LK, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  labs(x="Knowledge", y = "Visit number low quality feeder") +
  ylim(0, 2000) +
  theme_classic(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  
LK_visit_boxplot

HK_visit_boxplot <- ggplot(all_individuals_knowledge_sub2, aes(x = knowledge, y = visit_number_HK, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = colours) +
  labs(x="Knowledge", y = "Visit number high quality feeder") +
  theme_classic(base_size = 14)+
  ylim(0, 2000) +
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  
HK_visit_boxplot

LK_HK_visit_plot <- ggarrange(LK_visit_boxplot, HK_visit_boxplot, ncol = 2, nrow = 1,  widths = c(1, 1))
LK_HK_visit_plot

posterior_marginal <- all_individuals_knowledge_sub4_long %>%
  add_epred_draws(
    object = visit_no_brm1,
    re_formula = NA)

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge, feeder) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_draw_knowledge_comb$knowledge_feeder <- paste(posterior_draw_knowledge_comb$knowledge, posterior_draw_knowledge_comb$feeder, sep = "_")

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge, feeder) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  compare_levels(epred, by = knowledge_feeder)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = knowledge_feeder,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

write.csv(pairwise_contrasts_summary,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/pairwise_contrasts_summary.csv", row.names = FALSE)

#Halfeye plot of model predictions
visit_no_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = knowledge_feeder, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "darkgrey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Visit number", y = "Performance combination") 
visit_no_halfeye_plot

#Visit duration at both HK and LK feeders (filtered to include only data after manipulation)

visit_data_knowledge_post <- subset(visit_data_knowledge, day > 179)

# Make a new column for modelling
visit_data_knowledge_post$visit_duration2 <- visit_data_knowledge_post$visit_duration
visit_data_knowledge_post$visit_duration2 <- as.numeric(visit_data_knowledge_post$visit_duration2)

# Set censored values to an approximate detection limit
visit_data_knowledge_post$visit_duration2[visit_data_knowledge_post$visit_duration == 0] <- 1
summary(visit_data_knowledge_post$visit_duration)

# Create a censoring indicator
visit_data_knowledge_post$visit_duration_cens <- ifelse(visit_data_knowledge_post$visit_duration == 0, "left", "none")
table(visit_data_knowledge_post$visit_duration_cens)
visit_data_knowledge_post$visit_duration_cens <- factor(visit_data_knowledge_post$visit_duration_cens)

#Model just with prior
default_prior()

visit_dur_brm1_prior <- brm(visit_duration2 | cens(visit_duration_cens) ~ knowledge * quality + (1|JID), 
                           data = visit_data_knowledge_post, 
                           family = Gamma(link = "log"),
                           prior = c(
                             set_prior("normal(log(25), 1)", class = "Intercept"),  
                             set_prior("normal(0, 0.5)", class = "b"),
                             set_prior("gamma(2, 0.5)", class = "shape")),              
                           sample_prior = "only")

#Prior predictive checks 
pp_check(visit_dur_brm1_prior, ndraws = 1000)

#Model
visit_dur_brm1 <- brm(visit_duration2 | cens(visit_duration_cens) ~ knowledge * quality + (1|JID) + (1|feeder), 
                     data = visit_data_knowledge_post, 
                     family = Gamma(link = "log"),
                     prior = c(
                       set_prior("normal(log(25), 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b"),
                       set_prior("gamma(2, 0.5)", class = "shape")),              
                     
)

summary(visit_dur_brm1)

check_collinearity(visit_dur_brm1)

#Posterior predictive checks
pp_check(visit_dur_brm1, type = "dens_overlay", ndraws = 100) + scale_x_log10()

pp_check(visit_dur_brm1, type = "bars", ndraws = 1000)
pp_check(visit_dur_brm1, type = "hist", ndraws = 1000)

#Plot model
plot(visit_dur_brm1)
launch_shinystan(visit_dur_brm1)

#Posterior distribution
as_draws_df(visit_dur_brm1)
mcmc_areas(visit_dur_brm1)
mcmc_intervals(visit_dur_brm1)

#Evaluation and interpretation
loo(visit_dur_brm1, moment_match = TRUE)
fitted(visit_dur_brm1, scale = "response")
conditional_effects(visit_dur_brm1)
conditional_effects(visit_dur_brm1,effects = "knowledge:quality")
bayes_R2(visit_dur_brm1)
hypothesis(visit_dur_brm1, "knowledgenaive:qualitylow = 0") 
hypothesis(visit_dur_brm1, "knowledgenaive + knowledgenaive:qualitylow = 0")
emmeans(visit_dur_brm1, ~ knowledge * quality, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(visit_dur_brm1)

#Extract predictions
fitted(visit_dur_brm1)

posterior_marginal <- visit_data_knowledge_post %>%
  add_epred_draws(
    object = visit_dur_brm1,
    re_formula = NA)

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge, quality) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_draw_knowledge_comb$knowledge_quality <- paste(posterior_draw_knowledge_comb$knowledge, posterior_draw_knowledge_comb$quality, sep = "_")

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge, quality) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  compare_levels(epred, by = knowledge_quality)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = knowledge_quality,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary


#Halfeye plot of model predictions
visit_no_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = knowledge_feeder, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "darkgrey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Visit number", y = "Performance combination") 
visit_no_halfeye_plot


#(7.2) Displacements ---- 

#Displacements at both HK and LK

#Model knowledge_comb * quality

#visit_data_knowledge_post_sub should exclude instances in which two subsequent visits by the same individual (i.e., no displacement possible)
visit_data_knowledge_post_sub <- subset(visit_data_knowledge_post, !same_diff_indiv == "same")

#Model just with prior
default_prior()

displace_brm1_prior <- brm(displaced ~ knowledge * quality + (1 | mm(JID, next_JID)) + (1|feeder), 
                           data = visit_data_knowledge_post_sub, 
                           family = bernoulli(link = "logit"),
                           prior = c(
                             set_prior("normal(-2, 1)", class = "Intercept"),  
                             set_prior("normal(0, 0.5)", class = "b")),
                           sample_prior = "only"
)

#Prior predictive checks 
pp_check(displace_brm1_prior, ndraws = 1000)

displace_brm1 <- brm(displaced ~ knowledge * quality + (1|mm(JID, next_JID)) + (1|feeder), 
                     data = visit_data_knowledge_post_sub, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(-2, 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),
                     control = list(adapt_delta = 0.90))

summary(displace_brm1)

pairs(displace_brm1)

check_collinearity(displace_brm1)

#Posterior predictive checks
pp_check(displace_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(displace_brm1)
launch_shinystan(displace_brm1)

#Posterior distribution
as_draws_df(displace_brm1)
mcmc_areas(displace_brm1)
mcmc_intervals(displace_brm1)

#Evaluation and interpretation
loo(displace_brm1, moment_match = TRUE)
fitted(displace_brm1, scale = "response")
conditional_effects(displace_brm1)
conditional_effects(displace_brm1,effects = "knowledge_comb:quality")

bayes_R2(displace_brm1)
emmeans(displace_brm1, ~ knowledge_comb * quality, type = "response") |> pairs()
emm <- emmeans(displace_brm1, ~ knowledge_comb * quality, type = "response")
as.data.frame(emm)

#Sensitivity analysis and prior checks 
prior_summary(displace_brm1)

#Extract predictions
fitted(displace_brm1)


#Model knowledge_comb * quality

#visit_data_knowledge_post_sub should exclude instances in which two subsequent visits by the same individual (i.e., no displacement possible)
visit_data_knowledge_post_sub <- subset(visit_data_knowledge_post, !same_diff_indiv == "same")

#Model just with prior
default_prior()

displace_brm2_prior <- brm(displaced ~ knowledge_comb * quality + age_class_comb +
                           (1 | mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                           data = visit_data_knowledge_post_sub, 
                           family = bernoulli(link = "logit"),
                           sample_prior = "only"
)

#Prior predictive checks 
pp_check(displace_brm2_prior, ndraws = 1000)

displace_brm2 <- brm(displaced ~ knowledge_comb * quality + age_class_comb +
                     (1|mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                     data = visit_data_knowledge_post_sub, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(-2, 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),
                     control = list(adapt_delta = 0.90),
                     iter = 8000, warmup = 2000
)

visit_data_knowledge_post_sub$day_z <- scale(visit_data_knowledge_post_sub$day)

displace_glmm2 <- glmmTMB(displaced ~ knowledge_comb * quality + age_class_comb +
                     day_z + (1|JID) + (1|next_JID) + (1|dyad_ID) + (1|feeder), 
                     data = visit_data_knowledge_post_sub, 
                     family = binomial(link = "logit"))

displace_brm2.1 <- brm(displaced ~ knowledge_comb * quality + day_z * quality +
                     age_class_comb + (1|mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                     data = visit_data_knowledge_post_sub, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(-2, 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),
                     control = list(adapt_delta = 0.90),
                     iter = 8000, warmup = 2000
)

displace_brm2.2 <- brm(displaced ~ knowledge_comb * quality * day_z +
                         age_class_comb + (1|mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                       data = visit_data_knowledge_post_sub, 
                       family = bernoulli(link = "logit"),
                       prior = c(
                         set_prior("normal(-2, 1)", class = "Intercept"),  
                         set_prior("normal(0, 0.5)", class = "b")),
                       control = list(adapt_delta = 0.90),
                       iter = 8000, warmup = 2000
)

summary(displace_brm2)
summary(displace_brm2.1)
summary(displace_glmm2)
Anova(displace_glmm2)

check_collinearity(displace_brm2)

#Posterior predictive checks
pp_check(displace_brm2, type = "dens_overlay", ndraws = 100)

#Plot model
plot(displace_brm2)
launch_shinystan(displace_brm2)

#Posterior distribution
as_draws_df(displace_brm2)
mcmc_areas(displace_brm2)
mcmc_intervals(displace_brm2)

#Evaluation and interpretation
loo(displace_brm2, moment_match = TRUE)
fitted(displace_brm2, scale = "response")
conditional_effects(displace_brm2)
conditional_effects(displace_brm2.1)
conditional_effects(displace_brm2,effects = "knowledge_comb:quality")

bayes_R2(displace_brm2)
emmeans(displace_brm2, ~ knowledge_comb * quality, type = "response") |> pairs()
emmeans(displace_brm2.1, ~ day_z * quality, type = "response") |> pairs()

emm <- emmeans(displace_brm2, ~ knowledge_comb * quality, type = "response")
as.data.frame(emm)

hypothesis(displace_brm2.1, "qualitylow * day_z - day_z <  0") 
hypothesis(displace_brm2.1, "qualitylow * day_z > 0") 
hypothesis(displace_brm2.1, "day_z >  0") 

#Sensitivity analysis and prior checks 
prior_summary(displace_brm2)

#Extract predictions
fitted(displace_brm2)

posterior_marginal <- visit_data_knowledge_post_sub %>%
  add_epred_draws(
    object = displace_brm2,
    re_formula = NA,
    ndraws = 1000
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge_comb, quality) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_cells <- posterior_draw_knowledge_comb %>%
  mutate(cell = interaction(knowledge_comb, quality, sep = " | "))

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge_comb, quality) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_cells %>%
  group_by(.draw) %>%
  compare_levels(epred, by = cell) %>%
  rename(
    contrast = cell,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

write.csv(pairwise_contrasts_summary,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/pairwise_contrasts_summary.csv", row.names = FALSE)
write.csv(posterior_summary,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/posterior_summary.csv", row.names = FALSE)

posterior_draw_knowledge_comb_HK <- subset(posterior_draw_knowledge_comb, posterior_draw_knowledge_comb$quality == "high")
posterior_draw_knowledge_comb_LK <- subset(posterior_draw_knowledge_comb, posterior_draw_knowledge_comb$quality == "low")

#Halfeye plot of model predictions
displacement_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = knowledge_comb, x = epred, fill = quality)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  scale_fill_manual(values = c("low" = "lightgrey", "high" = "darkgrey")) +              
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of displacement", y = "Performance combination") +
  xlim(c(0, 0.5)) +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displacement_halfeye_plot

displacement_halfeye_plot_HK <- ggplot(posterior_draw_knowledge_comb_HK, aes(y = knowledge_comb, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "darkgrey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of displacement", y = "Performance combination") +
  xlim(c(0, 0.5)) +
  geom_vline(xintercept=c(	
    0.14535221), linetype="dotted") +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displacement_halfeye_plot_HK

displacement_halfeye_plot_LK <- ggplot(posterior_draw_knowledge_comb_LK, aes(y = knowledge_comb, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95), fill = "lightgrey") +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of displacement", y = "Performance combination") +
  xlim(c(0, 0.5)) +
  geom_vline(xintercept=c(	
    0.13489253), linetype="dotted") +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displacement_halfeye_plot_LK

ggarrange(displacement_halfeye_plot_HK, displacement_halfeye_plot_LK, widths = 1)


#Model knowledge * phase (LK)

#Model just with prior
default_prior()

displace_brm3_prior <- brm(displaced ~ knowledge_comb * period + age_class_comb +
                           (1|mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                           data = visit_data_knowledge_LK_sub, 
                           family = bernoulli(link = "logit"),
                           prior = c(
                             set_prior("normal(-2, 1)", class = "Intercept"),  
                             set_prior("normal(0, 0.5)", class = "b")),
                           sample_prior = "only"
)

#Prior predictive checks 
pp_check(displace_brm3_prior, ndraws = 1000)

displace_brm3 <- brm(displaced ~ knowledge_comb * period + age_class_comb +
                     (1|mm(JID, next_JID)) + (1|dyad_ID) + (1|feeder), 
                     data = visit_data_knowledge_LK_sub, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(-2, 1)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),
                     control = list(adapt_delta = 0.99))

summary(displace_brm3)

check_collinearity(displace_brm3)

#Posterior predictive checks
pp_check(displace_brm3, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(displace_brm3)
launch_shinystan(displace_brm3)

#Posterior distribution
as_draws_df(displace_brm3)
mcmc_areas(displace_brm3)
mcmc_intervals(displace_brm3)

#Evaluation and interpretation
loo(displace_brm3, moment_match = TRUE)
fitted(displace_brm3, scale = "response")
conditional_effects(displace_brm3)
conditional_effects(displace_brm3,effects = "knowledge:period")

bayes_R2(displace_brm3)
emmeans(displace_brm3, ~ knowledge * period, type = "response") |> pairs()
emm <- emmeans(displace_brm3, ~ knowledge * period, type = "response")
as.data.frame(emm)

#Sensitivity analysis and prior checks 
prior_summary(displace_brm3)


ce <- conditional_effects(displace_brm2, effects = "knowledge:period", re_formula = NA)
newdat <- as.data.frame(ce[["knowledge:period"]]) %>%
  dplyr::select(knowledge, period)
posterior <- posterior_epred(displace_brm2, newdata = newdat, re_formula = NA)

# Convert to long format
posterior_long <- as.data.frame(posterior) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "draw", values_to = "pred") %>%
  mutate(row = as.numeric(gsub("V", "", draw))) %>%
  left_join(newdat %>% mutate(row = row_number()), by = "row")

posterior_long$period <- factor(posterior_long$period, levels = c("preknowledge", "postknowledge"))
posterior_long$period_knowledge <- paste(posterior_long$knowledge, posterior_long$period, sep = "_")
posterior_long$period_knowledge <- factor(posterior_long$period_knowledge, levels = c("knowledgeable_preknowledge", "naive_preknowledge", "knowledgeable_postknowledge", "naive_postknowledge"))


displace_plot3 <- ggplot(posterior_long, aes(x = period_knowledge, y = pred, fill = knowledge)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values =c("black", "grey")) +
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8)) +
  labs(
    x = "Knowledge/period",
    y = "Predicted probability of displacement",
    fill = "Knowledge") +
  scale_x_discrete(labels=c('informed pre', 'naive pre', 'informed post', 'naive post'))

displace_plot3

#Extract predictions
fitted(displace_brm2)

posterior_marginal <- visit_data_knowledge_LK_sub %>%
  add_epred_draws(
    object = displace_brm2,
    re_formula = NA,
    ndraws = 1000
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge_comb, period) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_cells <- posterior_draw_knowledge_comb %>%
  mutate(cell = interaction(knowledge_comb, period, sep = " | "))

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge_comb, period) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_cells %>%
  group_by(.draw) %>%
  compare_levels(epred, by = cell) %>%
  rename(
    contrast = cell,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary



#(7.3) Network centrality (passive feeders) ----

#Model degree ~ knowledge * phase

mean_count <- mean(all_individuals_knowledge_sub2$degree)
var_count  <- var(all_individuals_knowledge_sub2$degree)
overdispersion_ratio <- var_count / mean_count
overdispersion_ratio

all_individuals_knowledge_sub2$obs <- row.names(all_individuals_knowledge_sub2)

all_individuals_knowledge_sub2$visit_number_passive_z <- scale(all_individuals_knowledge_sub2$visit_number_passive)

all_individuals_knowledge_sub2$visit_number_HK <- ifelse(!is.na(all_individuals_knowledge_sub2$visit_number_HK), all_individuals_knowledge_sub2$visit_number_HK, 0)
all_individuals_knowledge_sub2$visit_number_HK_z <- scale(all_individuals_knowledge_sub2$visit_number_HK)

all_individuals_knowledge_sub3 <- subset(all_individuals_knowledge_sub2, all_individuals_knowledge_sub2$visit_number_passive > 19)
all_individuals_knowledge_sub3 <- all_individuals_knowledge_sub3 %>%
  group_by(JID) %>%
  filter(n() > 1) %>%
  ungroup()

#Model just with prior
default_prior()

degree_brm1_prior <- brm(degree ~ knowledge * phase + s(visit_number_passive) + (1|JID), 
                           data = all_individuals_knowledge_sub2, 
                           family = negbinomial(link = "log"),
                           prior = c(
                            set_prior("normal(3, 0.5)", class = "Intercept"),  
                            set_prior("normal(0, 0.25)", class = "b"), 
                            set_prior("exponential(1)", class = "shape")),
                           sample_prior = "only"
)

#Prior predictive checks 
pp_check(degree_brm1_prior, ndraws = 1000)

degree_brm1.1 <- brm(degree ~ knowledge * phase + s(visit_number_passive) + (1|JID), 
                     data = all_individuals_knowledge_sub2, 
                     family = negbinomial(link = "log"),
                     prior = c(
                       set_prior("normal(3, 0.5)", class = "Intercept"),  
                       set_prior("normal(0, 0.25)", class = "b"), 
                       set_prior("exponential(1)", class = "shape")), 
                     iter = 6000, warmup = 2000,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE))

degree_brm1.2 <- brm(degree ~ knowledge * phase + s(visit_number_passive) + (1|JID), 
                     data = all_individuals_knowledge_sub2, 
                     family = poisson(link = "log"),
                     prior = c(
                       set_prior("normal(3, 0.5)", class = "Intercept"),  
                       set_prior("normal(0, 0.25)", class = "b")), 
                     iter = 6000, warmup = 2000,
                     control = list(adapt_delta = 0.99),
                     save_pars = save_pars(all = TRUE))

degree_brm1.3 <- brm(degree ~ knowledge * phase +
                     s(visit_number_passive_z) + (1|JID) + (1|obs), 
                     data = all_individuals_knowledge_sub2, 
                     family = poisson(link = "log"),
                     prior = c(
                       set_prior("normal(3, 0.5)", class = "Intercept"),  
                       set_prior("normal(0, 0.25)", class = "b")), 
                     iter = 6000, warmup = 2000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     save_pars = save_pars(all = TRUE))

#Current model in manuscript
degree_brm1.3.1 <- brm(degree ~ knowledge * phase + s(visit_number_passive) + (1|JID) 
                       + (1|obs), 
                       data = all_individuals_knowledge_sub2, 
                       family = poisson(link = "log"),
                       prior = c(
                         set_prior("normal(3, 0.5)", class = "Intercept"),  
                         set_prior("normal(0, 0.25)", class = "b")), 
                       iter = 6000, warmup = 2000,
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       save_pars = save_pars(all = TRUE))

#Exclude individuals with less than 10 visits
degree_brm1.3.1 <- brm(degree ~ knowledge * phase + s(visit_number_passive) + (1|JID) 
                     + (1|obs), 
                     data = all_individuals_knowledge_sub3, 
                     family = poisson(link = "log"),
                     prior = c(
                       set_prior("normal(3, 0.5)", class = "Intercept"),  
                       set_prior("normal(0, 0.25)", class = "b")), 
                     iter = 6000, warmup = 2000,
                     control = list(adapt_delta = 0.99, max_treedepth = 15),
                     save_pars = save_pars(all = TRUE))

degree_brm1 <- degree_brm1.1
degree_brm1 <- degree_brm1.2
degree_brm1 <- degree_brm1.3
degree_brm1 <- degree_brm1.3.1

summary(degree_brm1)

check_collinearity(degree_brm1)

#Posterior predictive checks
pp_check(degree_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(degree_brm1)
launch_shinystan(degree_brm1)

#Posterior distribution
as_draws_df(degree_brm1)
mcmc_areas(degree_brm1)
mcmc_intervals(degree_brm1)

#Evaluation and interpretation
degree_brm1.1_loo <- loo(degree_brm1.1, moment_match = TRUE)
degree_brm1.2_loo <- loo(degree_brm1.2, moment_match = TRUE)
degree_brm1.3_loo <- loo(degree_brm1.3, moment_match = TRUE)

loo(degree_brm1.1, degree_brm1.2, degree_brm1.3)
fitted(degree_brm1, scale = "response")
conditional_effects(degree_brm1)
conditional_effects(degree_brm1,effects = "knowledge:phase")
conditional_effects(degree_brm1,effects = "visit_number_HK_post:knowledge")

bayes_R2(degree_brm1)
emmeans(degree_brm1, ~ knowledge * phase, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(degree_brm1)

#Extract predictions
fitted(degree_brm1)

posterior_marginal <- all_individuals_knowledge_sub2 %>%
  add_epred_draws(
    object = degree_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge, phase) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_cells <- posterior_draw_knowledge_comb %>%
  mutate(cell = interaction(knowledge, phase, sep = " | "))

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge, phase) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_cells %>%
  group_by(.draw) %>%
  compare_levels(epred, by = cell) %>%
  rename(
    contrast = cell,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

degree_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = phase, x = epred, fill = knowledge)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  scale_fill_manual(values = c("naive" = "lightgrey", "knowledgeable" = "darkgrey")) +              
  theme_few(base_size = 14) +
  theme(
    legend.position = "top",
    text = element_text(family = "Garamond")) +
  labs(x = "Degree", y = "Performance combination") 
degree_halfeye_plot


#Model strength ~ knowledge * phase

#Model just with prior
default_prior()

strength_brm1_prior <- brm(strength ~ knowledge * phase + s(visit_number_passive) +(1|JID), 
                         data = all_individuals_knowledge_sub2, 
                         family = gaussian,
                         prior = c(
                           set_prior("normal(0.5, 0.1)", class = "Intercept"),  
                           set_prior("normal(0, 0.1)", class = "b"),
                           set_prior("exponential(1)", class = "sigma")), 
                         sample_prior = "only"
)

#Prior predictive checks 
pp_check(strength_brm1_prior, ndraws = 100)

strength_brm1 <- brm(strength ~ knowledge * phase + s(visit_number_passive) + (1|JID), 
                           data = all_individuals_knowledge_sub2, 
                           family = gaussian,
                           prior = c(
                             set_prior("normal(0.5, 0.1)", class = "Intercept"),  
                             set_prior("normal(0, 0.1)", class = "b"),
                             set_prior("exponential(1)", class = "sigma")), 
                           iter = 8000, warmup = 2000, control = list(adapt_delta = 0.999, max_treedepth = 15)
                           )
                           
summary(strength_brm1)

check_collinearity(strength_brm1)

#Posterior predictive checks
pp_check(strength_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(strength_brm1)
launch_shinystan(strength_brm1)

#Posterior distribution
as_draws_df(strength_brm1)
mcmc_areas(strength_brm1)
mcmc_intervals(strength_brm1)

#Evaluation and interpretation
loo(strength_brm1, moment_match = TRUE)
fitted(strength_brm1, scale = "response")
conditional_effects(strength_brm1)
conditional_effects(strength_brm1,effects = "knowledge:phase")

bayes_R2(strength_brm1)
emmeans(strength_brm1, ~ knowledge * phase, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(strength_brm1)

#Extract predictions
fitted(strength_brm1)

hypothesis(strength_brm1, "phasepre_knowledge >  0") 

posterior_marginal <- all_individuals_knowledge_sub2 %>%
  add_epred_draws(
    object = strength_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge, phase) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_cells <- posterior_draw_knowledge_comb %>%
  mutate(cell = interaction(knowledge, phase, sep = " | "))

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge, phase) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_cells %>%
  group_by(.draw) %>%
  compare_levels(epred, by = cell) %>%
  rename(
    contrast = cell,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

strength_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = phase, x = epred, fill = knowledge)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  scale_fill_manual(values = c("naive" = "lightgrey", "knowledgeable" = "darkgrey")) +              
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "strength", y = "Performance combination") 
strength_halfeye_plot


#Model degree change ~ knowledge * visit number HK post

#Model just with prior
default_prior()

degree_brm2_prior <- brm(degree_diff_abs ~ knowledge * visit_number_HK_post + s(visit_number_passive_diff_abs), 
                         data = all_individuals_knowledge_sub, 
                         family = gaussian(),
                         prior = c(
                             set_prior("normal(0, 0.25)", class = "Intercept"),  
                             set_prior("normal(0, 0.25)", class = "b")),
                         sample_prior = "only")

#Prior predictive checks 
pp_check(degree_brm2_prior, ndraws = 1000)

degree_brm2 <- brm(degree_diff_abs ~ knowledge * visit_number_HK_post + s(visit_number_passive_diff_abs), 
                   data = all_individuals_knowledge_sub, 
                   family = gaussian(),
                   prior = c(
                     set_prior("normal(0, 0.25)", class = "Intercept"),  
                     set_prior("normal(0, 0.25)", class = "b")))


#degree_brm2 <- brm(degree_diff_abs ~ HKLK_no_dual_events + s(visit_number_passive_post), 
#                   data = all_individuals_knowledgeable, 
#                   family = gaussian(),
#                   prior = c(
#                     set_prior("normal(0, 0.5)", class = "Intercept"),  
#                     set_prior("normal(0, 0.25)", class = "b")))


summary(degree_brm2)

check_collinearity(degree_brm2)

#Posterior predictive checks
pp_check(degree_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(degree_brm2)
launch_shinystan(degree_brm2)

#Posterior distribution
as_draws_df(degree_brm2)
mcmc_areas(degree_brm2)
mcmc_intervals(degree_brm2)

#Evaluation and interpretation
loo(degree_brm2, moment_match = TRUE)
fitted(degree_brm2, scale = "response")
conditional_effects(degree_brm2)
conditional_effects(degree_brm2,effects = "visit_number_HK_post:knowledge")

bayes_R2(degree_brm2)
emmeans(degree_brm2, ~ knowledge * phase, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(degree_brm2)

#Extract predictions
fitted(degree_brm2)


#Model strength change ~ knowledge * visit number HK post

#Model just with prior
default_prior()

strength_brm2_prior <- brm(strength_diff_abs ~ knowledge * visit_number_knowledge_post + s(visit_number_passive_diff_abs), 
                         data = all_individuals_knowledge_sub, 
                         family = gaussian(),
                         prior = c(
                           set_prior("normal(0, 0.15)", class = "Intercept"),  
                           set_prior("normal(0, 0.15)", class = "b")),
                         sample_prior = "only")

#Prior predictive checks 
pp_check(strength_brm2_prior, ndraws = 1000)

strength_brm2 <- brm(strength_diff_abs ~ knowledge * visit_number_knowledge_post + s(visit_number_passive_diff_abs), 
                   data = all_individuals_knowledge_sub, 
                   family = gaussian(),
                   prior = c(
                     set_prior("normal(0, 0.15)", class = "Intercept"),  
                     set_prior("normal(0, 0.15)", class = "b")),
                   sample_prior = "only")


strength_brm2 <- brm(strength_diff_abs ~ HKLK_no_dual_events + s(visit_number_passive_post), 
                   data = all_individuals_knowledgeable, 
                   family = gaussian(),
                   prior = c(
                     set_prior("normal(0, 0.5)", class = "Intercept"),  
                     set_prior("normal(0, 0.25)", class = "b")))


summary(strength_brm2)

check_collinearity(strength_brm2)

#Posterior predictive checks
pp_check(strength_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(strength_brm2)
launch_shinystan(strength_brm2)

#Posterior distribution
as_draws_df(strength_brm2)
mcmc_areas(strength_brm2)
mcmc_intervals(strength_brm2)

#Evaluation and interpretation
loo(strength_brm2, moment_match = TRUE)
fitted(strength_brm2, scale = "response")
conditional_effects(strength_brm2)
conditional_effects(strength_brm2,effects = "knowledge:phase")
conditional_effects(strength_brm2,effects = "visit_number_HK_post:knowledge")

bayes_R2(strength_brm2)
emmeans(strength_brm2, ~ knowledge * phase, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(strength_brm2)

#Extract predictions
fitted(strength_brm2)


#(7.4) Dyadic events and relationships ----

#Are new dyadic associations at network feeders formed after dyadic 
  #events at experimental feeders?  

#Only dyads with dyadic events at HKLK feeders that had no previous/recent association at network feeders
dual_events_knowledge_dyads_sub <- subset(dual_events_knowledge_dyads, is.na(dual_events_knowledge_dyads$edge_weight_pre))
dual_events_knowledge_dyads_sub$HKLK_no_dual_events_mean_z <- scale(dual_events_knowledge_dyads_sub$HKLK_no_dual_events_mean)
dual_events_knowledge_dyads_sub <- subset(dual_events_knowledge_dyads_sub, !HKLK_knowledge_comb == "knowledge knowledge")

dual_events_knowledge_dyads_sub <- dual_events_knowledge_dyads_sub %>%
  group_by(HKLK_knowledge_comb4) %>%
  mutate(HKLK_no_dual_events_mean_gc = HKLK_no_dual_events_mean - mean(HKLK_no_dual_events_mean),
         HKLK_no_dual_events_mean_gc_z = scale(HKLK_no_dual_events_mean_gc))

dual_events_knowledge_dyads_sub <- dual_events_knowledge_dyads_sub %>%
  group_by(HKLK_knowledge_comb) %>%
  mutate(HKLK_no_dual_events_mean_gc = HKLK_no_dual_events_mean - mean(HKLK_no_dual_events_mean),
         HKLK_no_dual_events_mean_gc_z = scale(HKLK_no_dual_events_mean_gc))

dual_events_knowledge_dyads_sub$HKLK_no_dual_events_mean_z

dual_events_knowledge_dyads_sub <- dual_events_knowledge_dyads_sub %>%
                              separate_wider_delim(
                              col = dyad_ID,
                              delim = " ",
                              names = c("JID1", "JID2")
)

dual_events_knowledge_dyads_sub$phase <- "postknowledge"
dual_events_knowledge_dyads_sub$dyad_ID_phase <- paste(dual_events_knowledge_dyads_sub$dyad_ID, dual_events_knowledge_dyads_sub$phase, sep = "_")
dual_events_knowledge_dyads_sub$edge_cooccurrence <- edge_list_knowledge$edge[match(dual_events_knowledge_dyads_sub$dyad_ID_phase, edge_list_knowledge$dyad_ID_phase)]
dual_events_knowledge_dyads_sub$edge_binary <- ifelse(dual_events_knowledge_dyads_sub$edge_cooccurrence > 0, 1, 0)

dual_events_knowledge_dyads_sub$JID1_phase <- paste(dual_events_knowledge_dyads_sub$JID1, dual_events_knowledge_dyads_sub$phase, sep = "_")
dual_events_knowledge_dyads_sub$JID2_phase <- paste(dual_events_knowledge_dyads_sub$JID2, dual_events_knowledge_dyads_sub$phase, sep = "_")

dual_events_knowledge_dyads_sub$ID1_events <- edge_list_knowledge$ID1_events[match(dual_events_knowledge_dyads_sub$JID1_phase, edge_list_knowledge$ID1_phase)]
dual_events_knowledge_dyads_sub$ID2_events <- edge_list_knowledge$ID2_events[match(dual_events_knowledge_dyads_sub$JID2_phase, edge_list_knowledge$ID2_phase)]

dual_events_knowledge_dyads_sub$ID1_events_z <- scale(dual_events_knowledge_dyads_sub$ID1_events)
dual_events_knowledge_dyads_sub$ID2_events_z <- scale(dual_events_knowledge_dyads_sub$ID2_events)

plot(dual_events_knowledge_dyads_sub$edge_new, dual_events_knowledge_dyads_sub$edge_binary)

#Model just with prior
default_prior()

edge_brm1_prior <- brm(edge_new ~ HKLK_knowledge_comb * HKLK_no_dual_events_mean_gc_z 
                           + (1 | mm(JID1, JID2)), 
                           data = dual_events_knowledge_dyads_sub, 
                           family = bernoulli(link = "logit"),
                           prior = c(
                             set_prior("normal(0, 0.5)", class = "b"),        
                             set_prior("normal(-2, 0.75)", class = "Intercept"),
                             set_prior("student_t(3, 0, 2.5)", class = "sd")),
                       sample_prior = "only"
)
                             
#Prior predictive checks 
pp_check(edge_brm1_prior, ndraws = 100)

#Model currently used in manuscript 
edge_brm1 <- brm(edge_new ~ HKLK_knowledge_comb * HKLK_no_dual_events_mean_gc_z
                     + (1 | mm(JID1, JID2)), data = dual_events_knowledge_dyads_sub, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(0, 0.5)", class = "b"),        
                       set_prior("normal(-2, 0.75)", class = "Intercept")),
                     iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

edge_brm1 <- brm(edge_new ~ HKLK_knowledge_comb * HKLK_no_dual_events_mean_gc_z +
                 ID1_events_z + ID2_events_z + (1 | mm(JID1, JID2)), 
                 data = dual_events_knowledge_dyads_sub, 
                 family = bernoulli(link = "logit"),
                 prior = c(
                   set_prior("normal(0, 0.5)", class = "b"),        
                   set_prior("normal(-2, 0.75)", class = "Intercept")),
                 iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

edge_brm1 <- brm(edge_binary ~ HKLK_knowledge_comb * HKLK_no_dual_events_mean_gc_z
                 + ID1_events_z + ID2_events_z + (1 | mm(JID1, JID2)), data = dual_events_knowledge_dyads_sub, 
                 family = bernoulli(link = "logit"),
                 prior = c(
                   set_prior("normal(0, 0.5)", class = "b"),        
                   set_prior("normal(-2, 0.75)", class = "Intercept")),
                 iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

edge_brm1 <- brm(edge_binary ~ knowledge_comb * HKLK_dyadic_events
                 + ID1_events + ID2_events + (1 | mm(ID1, ID2)), data = edge_list_know_sub, 
                 family = bernoulli(link = "logit"),
                 prior = c(
                   set_prior("normal(0, 0.5)", class = "b"),        
                   set_prior("normal(-2, 0.75)", class = "Intercept")),
                 iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

edge_brm1 <- brm(edge_new ~ HKLK_knowledge_comb 
                 + (1 | mm(JID1, JID2)) + (1|HKLK_dyad_ID), data = dual_events_knowledge_dyads_sub, 
                 family = bernoulli(link = "logit"),
                 prior = c(
                   set_prior("normal(0, 0.5)", class = "b"),        
                   set_prior("normal(-2, 0.75)", class = "Intercept")),
                 iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

edge_brm1 <- brm(edge_new ~ HKLK_knowledge_comb4 * HKLK_no_dual_events_mean_gc_z 
                 + (1|mm(JID1, JID2)), data = dual_events_knowledge_dyads_sub, 
                 family = bernoulli(link = "logit"),
                 prior = c(
                   set_prior("normal(0, 0.5)", class = "b"),        
                   set_prior("normal(-2, 0.75)", class = "Intercept")),
                 iter = 8000, warmup = 2000, control = list(adapt_delta = 0.95)
                 
)

summary(edge_brm1)

check_collinearity(edge_brm1)

#Posterior predictive checks
pp_check(edge_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(edge_brm1)
launch_shinystan(edge_brm1)

#Posterior distribution
as_draws_df(edge_brm1)
mcmc_areas(edge_brm1)
mcmc_intervals(edge_brm1)

#Evaluation and interpretation
loo(edge_brm1, moment_match = TRUE)
fitted(edge_brm1, scale = "response")
conditional_effects(edge_brm1)
ce <- conditional_effects(edge_brm1,effects = "HKLK_knowledge_comb4")

ce_df <- ce$HKLK_knowledge_comb4

ggplot(ce_df, aes(x = HKLK_knowledge_comb4, y = estimate__, ymin = lower__, ymax = upper__)) +
  geom_point() +
  geom_errorbar(width = 0.05) +
  ylim(0, 0.25) +
  labs(y = "Predicted probability", x = "Knowledge level") +
  theme_minimal()

conditional_effects(edge_brm1,effects = "HKLK_no_dual_events_mean_gc_z:HKLK_knowledge_comb4")
conditional_effects(edge_brm1,effects = "HKLK_no_dual_events_mean_gc_z")

bayes_R2(edge_brm1)
emmeans(edge_brm1, ~ HKLK_knowledge_comb, type = "response") |> pairs()
emmeans(edge_brm1, ~ HKLK_knowledge_comb) |> pairs()
pairs(emmeans(edge_brm1, ~ HKLK_knowledge_comb4, type = "response"))
hypothesis(edge_brm1, "HKLK_knowledge_combknowledgenaive - HKLK_knowledge_combnaivenaive > 0")
hypothesis(edge_brm1, "HKLK_knowledge_combknowledgenaive > 0")
hypothesis(edge_brm1, "HKLK_knowledge_combnaivenaive > 0")

hypothesis(edge_brm1, "HKLK_no_dual_events_mean_gc_z > 0")

hypothesis(edge_brm1, "HKLK_knowledge_comb4naive < 0")
hypothesis(edge_brm1, "HKLK_knowledge_comb4naive < 0")

emmeans(edge_brm1, ~ HKLK_knowledge_comb:HKLK_no_dual_events_mean_gc_z, type = "response") |> pairs()
emmeans(edge_brm1, ~ knowledge_comb, type = "response") |> pairs()

contrast(emm, interaction = "pairwise")

#Sensitivity analysis and prior checks 
prior_summary(edge_brm1)

#Extract predictions
fitted(edge_brm1)

#Posterior marginal 
posterior_marginal <- dual_events_knowledge_dyads_sub %>%
  add_epred_draws(
    object = edge_brm1,
    re_formula = NA
  )

posterior_marginal <- edge_list_know_sub %>%
  add_epred_draws(
    object = edge_brm1,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, HKLK_knowledge_comb) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge_comb) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(HKLK_knowledge_comb) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge_comb) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  compare_levels(epred, by = HKLK_knowledge_comb)

pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  compare_levels(epred, by = knowledge_comb)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = HKLK_knowledge_comb,
    diff = epred
  )

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = knowledge_comb,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Halfeye plot of model predictions
edge_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = HKLK_knowledge_comb, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of new association", y = "Performance combination") +
  xlim(c(0, 0.2)) +
  geom_vline(xintercept=c(	
    0.04572497), linetype="dotted") +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive naive'))
edge_halfeye_plot

edge_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = knowledge_comb, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of new edge", y = "Performance combination") +
  xlim(c(0, 1)) +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive naive'))
edge_halfeye_plot


ce <- conditional_effects(edge_brm1, effects = "HKLK_knowledge_comb4", re_formula = NA)
newdat <- as.data.frame(ce[["HKLK_knowledge_comb4"]]) %>%
  dplyr::select(HKLK_knowledge_comb4, HKLK_no_dual_events_mean_gc_z)
posterior <- posterior_epred(edge_brm1, newdata = newdat, re_formula = NA)

# Convert to long format
posterior_long <- as.data.frame(posterior) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "draw", values_to = "pred") %>%
  mutate(row = as.numeric(gsub("V", "", draw))) %>%
  left_join(newdat %>% mutate(row = row_number()), by = "row")

edge_plot1 <- ggplot(posterior_long, aes(x = HKLK_knowledge_comb4, y = pred, fill = HKLK_knowledge_comb4)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values =c("black", "grey")) +
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8)) +
  labs(
    x = "Knowledge combination",
    y = "Predicted probability of 'new' association")+
  scale_x_discrete(labels=c('informed', 'naive')) +
  ylim(c(0, 0.25))
edge_plot1


#How do pre-existing associations change 

dyads_prepost_sub$HKLK_no_dual_events_binary_c <- dyads_prepost_sub$HKLK_no_dual_events_binary - 0.5

dyads_prepost_sub <- dyads_prepost_sub %>%
  group_by(knowledge_comb4) %>%
  mutate(HKLK_no_dual_events_mean_gc = HKLK_no_dual_events_mean - mean(HKLK_no_dual_events_mean),
         HKLK_no_dual_events_mean_gc_z = scale(HKLK_no_dual_events_mean_gc))

dyads_prepost_sub$HKLK_no_dual_events_mean_z <- scale(dyads_prepost_sub$HKLK_no_dual_events_mean)

dyads_prepost_sub$g_edge_weights_pre_z <- scale(dyads_prepost_sub$g_edge_weights_pre)

dyads_prepost_sub$edge_diff_resid_10 <- dyads_prepost_sub$edge_diff_resid * 10

dyads_prepost_sub2 <- subset(dyads_prepost_sub, dyads_prepost_sub$HKLK_no_dual_events_mean >0)

dyads_prepost_sub2 <- subset(dyads_prepost_sub2, dyads_prepost_sub2$HKLK_no_dual_events_mean_gc_z < 5)

#Model just with prior
default_prior()

edge_brm2_prior <- brm(edge_diff_resid ~ knowledge_comb2 * HKLK_no_dual_events_mean_gc_z,
                       data = dyads_prepost_sub2, 
                       family = gaussian(),
                       prior = c(
                         set_prior("normal(0, 0.025)", class = "b"),        
                         set_prior("normal(0, 0.025)", class = "Intercept"),
                         set_prior("exponential(50)", class = "sigma")),
                         sample_prior = "only"
)

#Prior predictive checks 
pp_check(edge_brm2_prior, ndraws = 100)

#model currently in manuscript
edge_brm2 <- brm(g_edge_weights_post ~ g_edge_weights_pre + knowledge_comb2 * HKLK_no_dual_events_mean_gc_z, 
                 data = dyads_prepost_sub, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.02)", class = "b"),        
                   set_prior("normal(0, 0.02)", class = "Intercept"),
                   set_prior("exponential(200)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

edge_brm2 <- brm(edge_diff_resid ~ knowledge_comb2 * HKLK_no_dual_events_mean_gc_z, 
                 data = dyads_prepost_sub, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.02)", class = "b"),        
                   set_prior("normal(0, 0.02)", class = "Intercept"),
                   set_prior("exponential(200)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

edge_brm2 <- brm(edge_diff_resid ~ knowledge_comb4 * HKLK_no_dual_events_mean_gc_z, 
                       data = dyads_prepost_sub, 
                       family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.02)", class = "b"),        
                   set_prior("normal(0, 0.02)", class = "Intercept"),
                   set_prior("exponential(200)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

edge_brm2 <- brm(edge_diff_resid_10 ~ knowledge_comb4 * HKLK_no_dual_events_mean_gc_z, 
                 data = dyads_prepost_sub, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.2)", class = "b"),        
                   set_prior("normal(0, 0.2)", class = "Intercept"),
                   set_prior("exponential(2000)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

#Remove influential datapoints 
dyads_prepost_sub2 <- subset(dyads_prepost_sub, dyads_prepost_sub$HKLK_no_dual_events_mean_gc_z < 10)

edge_brm2 <- brm(edge_diff_resid_10 ~ knowledge_comb4 * HKLK_no_dual_events_mean_gc_z, 
                 data = dyads_prepost_sub2, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.2)", class = "b"),        
                   set_prior("normal(0, 0.2)", class = "Intercept"),
                   set_prior("exponential(2000)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

edge_brm2 <- brm(g_edge_weights_post ~ g_edge_weights_pre + knowledge_comb2 * HKLK_no_dual_events_mean_gc_z, 
                 data = dyads_prepost_sub2, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.02)", class = "b"),        
                   set_prior("normal(0, 0.02)", class = "Intercept"),
                   set_prior("exponential(200)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha"))
)

#Only dyads that engaged in dyadic events at experimental locations
dyads_prepost_sub2 <- subset(dyads_prepost_sub, dyads_prepost_sub$HKLK_no_dual_events_mean > 0)

#Remove dyads that engaged in many dyadic events
dyads_prepost_sub3 <- subset(dyads_prepost_sub2, dyads_prepost_sub2$HKLK_no_dual_events_mean < 40)
dyads_prepost_sub3 <- subset(dyads_prepost_sub2, dyads_prepost_sub2$HKLK_no_dual_events_mean_gc_z < 5)

#Remove dyads with high edge weight 
dyads_prepost_sub3 <- subset(dyads_prepost_sub2, dyads_prepost_sub2$g_edge_weights_post < 0.15)
dyads_prepost_sub3 <- subset(dyads_prepost_sub3, dyads_prepost_sub3$g_edge_weights_post < 0.15)

#Current model in manuscript 
edge_brm2 <- brm(g_edge_weights_post ~ g_edge_weights_pre_z + knowledge_comb2 * 
                 HKLK_no_dual_events_mean_gc_z + (1 | mm(ind1_ID, ind2_ID)), 
                 data = dyads_prepost_sub3, 
                 family = skew_normal(),
                 prior = c(
                   set_prior("normal(0, 0.02)", class = "b"),        
                   set_prior("normal(0, 0.02)", class = "Intercept"),
                   set_prior("exponential(200)", class = "sigma"),
                   set_prior("normal(-2, 1)", class = "alpha")),
                 iter = 8000, warmup = 2000
)

summary(edge_brm2)

check_collinearity(edge_brm2)

#Posterior predictive checks
pp_check(edge_brm2, type = "dens_overlay", ndraws = 100)
pp_check(edge_brm2, type = "scatter_avg")

#Plot model
plot(edge_brm2)
launch_shinystan(edge_brm2)

#Posterior distribution
as_draws_df(edge_brm2)
mcmc_areas(edge_brm2)
mcmc_intervals(edge_brm2)

#Evaluation and interpretation
loo(edge_brm2, moment_match = TRUE)
fitted(edge_brm2, scale = "response")
conditional_effects(edge_brm2)
conditional_effects(edge_brm2,effects = "knowledge_comb2")
conditional_effects(edge_brm2,effects = "knowledge_comb4")
conditional_effects(edge_brm2,effects = "HKLK_no_dual_events_mean_gc_z:knowledge_comb2")
conditional_effects(edge_brm2,effects = "HKLK_no_dual_events_binary_c:knowledge_comb4")
conditional_effects(edge_brm2,effects = "HKLK_no_dual_events_mean_gc_z:knowledge_comb4")

bayes_R2(edge_brm2)
emmeans(edge_brm2, ~ knowledge_comb2, type = "response") |> pairs()
emmeans(edge_brm2, ~ knowledge_comb2:HKLK_no_dual_events_mean_gc_z, type = "response") |> pairs()
hypothesis(edge_brm2, "knowledge_comb2knowledgenaive - knowledge_comb2naivenaive > 0")
hypothesis(edge_brm2, "knowledge_comb2naivenaive > 0")
hypothesis(edge_brm2, "knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z < 0")
hypothesis(edge_brm2, "knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z - knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z < 0")
hypothesis(edge_brm2, "knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z < 0")
hypothesis(edge_brm2, "knowledge_comb4naive:HKLK_no_dual_events_mean_gc_z < 0")

#Do slopes differ from 0
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z + knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z + knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z > 0")

#Do slopes differ from each other
hypothesis(edge_brm2, "knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z - knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z - knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z - knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z > 0")

hypothesis(edge_brm2, "knowledge_comb2knowledgenaive:HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z > 0")

hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z + knowledge_comb2naivenaive:HKLK_no_dual_events_mean_gc_z > 0")

hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z > 0")
hypothesis(edge_brm2, "HKLK_no_dual_events_mean_gc_z + knowledge_comb4naive:HKLK_no_dual_events_mean_gc_z > 0")

#Sensitivity analysis and prior checks 
prior_summary(edge_brm2)

#Extract predictions
fitted(edge_brm2)

#Posterior marginal 
posterior_marginal <- dyads_prepost_sub2 %>%
  add_epred_draws(
    object = edge_brm2,
    re_formula = NA
  )

#Posterior summary for text reporting 
posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge_comb2) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_summary <- posterior_draw_knowledge_comb %>%
  group_by(knowledge_comb2) %>%
  median_qi(epred, .width = c(0.5, 0.8, 0.95))
posterior_summary

#Contrasts for text reporting
pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  compare_levels(epred, by = knowledge_comb2)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = knowledge_comb2,
    diff = epred
  )

pairwise_contrasts_summary <- pairwise_contrasts %>%
  group_by(contrast) %>%
  summarise(
    median = median(diff),
    mean = mean(diff),
    lower_80 = quantile(diff, 0.1),
    upper_80 = quantile(diff, 0.9),
    lower_95 = quantile(diff, 0.025),
    upper_95 = quantile(diff, 0.975),
    P_gt_0 = mean(diff > 0),
    .groups = "drop"
  ) %>%
  arrange(contrast)
pairwise_contrasts_summary

#Halfeye plot of model predictions
edge_halfeye_plot <- ggplot(posterior_draw_knowledge_comb, aes(y = HKLK_knowledge_comb, x = epred)) +
  stat_halfeye(.width = c(0.5, 0.8, 0.95)) +
  theme_few(base_size = 14) +
  theme(
    legend.position = "none",
    text = element_text(family = "Garamond")) +
  labs(x = "Probability of new association", y = "Performance combination") +
  xlim(c(0, 0.2)) +
  geom_vline(xintercept=c(	
    0.04572497), linetype="dotted") +
  scale_y_discrete(labels=c('informed informed', 'informed naive', 'naive naive'))
edge_halfeye_plot


# Conditional effects for interaction
ce <- conditional_effects(edge_brm2, effects = "HKLK_no_dual_events_mean_gc_z:knowledge_comb4", spaghetti = FALSE)
# Extract the first effect (categorical x continuous)
ce_df <- ce[[1]]

# Set colors for the two groups
cols <- c("naive" = "grey", "knowledge" = "black")  # change labels if your factor levels differ

ggplot() +
  # Raw data points
  geom_point(data = dyads_prepost_sub,
             aes(x = HKLK_no_dual_events_mean_gc_z, y = edge_diff_resid * 10, color = knowledge_comb4),
             alpha = 0.6) +
  # Regression lines
  geom_line(data = ce_df,
            aes(x = HKLK_no_dual_events_mean_gc_z, y = estimate__, color = knowledge_comb4),
            size = 1.2) +
  # Shaded 95% CI
  geom_ribbon(data = ce_df,
              aes(x = HKLK_no_dual_events_mean_gc_z, ymin = lower__, ymax = upper__, fill = knowledge_comb4),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_minimal() +
  labs(x = "Continuous predictor", y = "Residualized change", color = "Group", fill = "Group") +
  theme(legend.position = "top")




edge_m1 <- lmer(edge_diff_abs ~ knowledge_comb2 * HKLK_no_dual_events_mean + (1|ind1_ID) + (1|ind2_ID), data = dyads_prepost_sub)
edge_m1 <- lmer(edge_diff_resid ~ knowledge_comb2 * HKLK_no_dual_events_mean + (1|ind1_ID) + (1|ind2_ID), data = dyads_prepost_sub)

edge_m1 <- glmmTMB(edge_new ~ HKLK_knowledge_comb * HKLK_no_dual_events_mean, data = dual_events_knowledge_dyads[is.na(dual_events_knowledge_dyads$edge_weight_pre),], family = binomial)
edge_m1 <- glmmTMB(edge_new ~ HKLK_knowledge_comb + (1|dyad_ID), data = dual_events_knowledge_dyads[is.na(dual_events_knowledge_dyads$edge_weight_pre),], family = binomial)

summary(edge_m1)
Anova(edge_m1, type = "II")

library(interactions)
interact_plot(model = edge_m1, pred = HKLK_no_dual_events_mean,
              modx = knowledge_comb2)

emmeans(edge_m1, ~ knowledge_comb2 * HKLK_no_dual_events_mean, type = "response") |> pairs()


#Dual/following events (passive feeders)

#Dyadic relationships

#(4) Public information use/dual events (experimental feeders)





#Create dataframe with dyads and edge weights

#Preknowledge
g_dyads_pre <- igraph::ends(g_preknowledge, igraph::E(g_preknowledge), names=TRUE)
g_edge_weights_pre <- igraph::E(g_preknowledge)$weight

dyad_edge_pre <- data.frame(matrix(NA, nrow = 1968, ncol = 2))
dyad_edge_pre$X1 <- g_dyads_pre
dyad_edge_pre$X2 <- g_edge_weights_pre

dyad_edge_pre$ind1_ID <- dyad_edge_pre$X1[,1]
dyad_edge_pre$ind2_ID <- dyad_edge_pre$X1[,2]

dyad_edge_pre$g_dyad_ID_pre <- paste(dyad_edge_pre$ind1_ID, dyad_edge_pre$ind2_ID)
dyad_edge_pre$g_edge_weights_pre <- dyad_edge_pre$X2

dyad_edge_pre <- dyad_edge_pre[, c("g_dyad_ID_pre", "g_edge_weights_pre", "ind1_ID", "ind2_ID")]

dyad_edge_pre$ind1_knowledge <- all_individuals_total$knowledge[match(dyad_edge_pre$ind1_ID, all_individuals_total$JID)]
dyad_edge_pre$ind2_knowledge <- all_individuals_total$knowledge[match(dyad_edge_pre$ind2_ID, all_individuals_total$JID)]

dyad_edge_pre$knowledge_comb <- paste(dyad_edge_pre$ind1_knowledge, dyad_edge_pre$ind2_knowledge, sep = " ")

dyad_edge_pre$knowledge_comb2 <- ifelse(dyad_edge_pre$knowledge_comb == "knowledge knowledge", "knowledge knowledge", ifelse(dyad_edge_pre$knowledge_comb == "naive naive", "naive naive", "knowledge naive"))


#Postknowledge
g_dyads_post <- igraph::ends(g_postknowledge, igraph::E(g_postknowledge), names=TRUE)
g_edge_weights_post <- igraph::E(g_postknowledge)$weight

dyad_edge_post <- data.frame(matrix(NA, nrow = 2965, ncol = 2))
dyad_edge_post$X1 <- g_dyads_post
dyad_edge_post$X2 <- g_edge_weights_post

dyad_edge_post$ind1_ID <- dyad_edge_post$X1[,1]
dyad_edge_post$ind2_ID <- dyad_edge_post$X1[,2]

dyad_edge_post$g_dyad_ID_post <- paste(dyad_edge_post$ind1_ID, dyad_edge_post$ind2_ID)
dyad_edge_post$g_edge_weights_post <- dyad_edge_post$X2

dyad_edge_post <- dyad_edge_post[, c("g_dyad_ID_post", "g_edge_weights_post", "ind1_ID", "ind2_ID")]

dyad_edge_post$ind1_knowledge <- all_individuals_total$knowledge[match(dyad_edge_post$ind1_ID, all_individuals_total$JID)]
dyad_edge_post$ind2_knowledge <- all_individuals_total$knowledge[match(dyad_edge_post$ind2_ID, all_individuals_total$JID)]

dyad_edge_post$knowledge_comb <- paste(dyad_edge_post$ind1_knowledge, dyad_edge_post$ind2_knowledge, sep = " ")

dyad_edge_post$knowledge_comb2 <- ifelse(dyad_edge_post$knowledge_comb == "knowledge knowledge", "knowledge knowledge", ifelse(dyad_edge_post$knowledge_comb == "naive naive", "naive naive", "knowledge naive"))


#Create dataframe with dyads and edge weights

#Preknowledge
g_dyads_pre <- igraph::ends(g_preknowledge, igraph::E(g_preknowledge), names=TRUE)
g_edge_weights_pre <- igraph::E(g_preknowledge)$weight

dyad_edge_pre <- data.frame(matrix(NA, nrow = 1968, ncol = 2))
dyad_edge_pre$X1 <- g_dyads_pre
dyad_edge_pre$X2 <- g_edge_weights_pre

dyad_edge_pre$ind1_ID <- dyad_edge_pre$X1[,1]
dyad_edge_pre$ind2_ID <- dyad_edge_pre$X1[,2]

dyad_edge_pre$g_dyad_ID_pre <- paste(dyad_edge_pre$ind1_ID, dyad_edge_pre$ind2_ID)
dyad_edge_pre$g_edge_weights_pre <- dyad_edge_pre$X2

dyad_edge_pre <- dyad_edge_pre[, c("g_dyad_ID_pre", "g_edge_weights_pre", "ind1_ID", "ind2_ID")]

dyad_edge_pre$ind1_knowledge <- all_individuals_total$knowledge[match(dyad_edge_pre$ind1_ID, all_individuals_total$JID)]
dyad_edge_pre$ind2_knowledge <- all_individuals_total$knowledge[match(dyad_edge_pre$ind2_ID, all_individuals_total$JID)]

dyad_edge_pre$knowledge_comb <- paste(dyad_edge_pre$ind1_knowledge, dyad_edge_pre$ind2_knowledge, sep = " ")

dyad_edge_pre$knowledge_comb2 <- ifelse(dyad_edge_pre$knowledge_comb == "knowledge knowledge", "knowledge knowledge", ifelse(dyad_edge_pre$knowledge_comb == "naive naive", "naive naive", "knowledge naive"))


#Postknowledge
g_dyads_post <- igraph::ends(g_postknowledge, igraph::E(g_postknowledge), names=TRUE)
g_edge_weights_post <- igraph::E(g_postknowledge)$weight

dyad_edge_post <- data.frame(matrix(NA, nrow = 2965, ncol = 2))
dyad_edge_post$X1 <- g_dyads_post
dyad_edge_post$X2 <- g_edge_weights_post

dyad_edge_post$ind1_ID <- dyad_edge_post$X1[,1]
dyad_edge_post$ind2_ID <- dyad_edge_post$X1[,2]

dyad_edge_post$g_dyad_ID_post <- paste(dyad_edge_post$ind1_ID, dyad_edge_post$ind2_ID)
dyad_edge_post$g_edge_weights_post <- dyad_edge_post$X2

dyad_edge_post <- dyad_edge_post[, c("g_dyad_ID_post", "g_edge_weights_post", "ind1_ID", "ind2_ID")]

dyad_edge_post$ind1_knowledge <- all_individuals_total$knowledge[match(dyad_edge_post$ind1_ID, all_individuals_total$JID)]
dyad_edge_post$ind2_knowledge <- all_individuals_total$knowledge[match(dyad_edge_post$ind2_ID, all_individuals_total$JID)]

dyad_edge_post$knowledge_comb <- paste(dyad_edge_post$ind1_knowledge, dyad_edge_post$ind2_knowledge, sep = " ")

dyad_edge_post$knowledge_comb2 <- ifelse(dyad_edge_post$knowledge_comb == "knowledge knowledge", "knowledge knowledge", ifelse(dyad_edge_post$knowledge_comb == "naive naive", "naive naive", "knowledge naive"))


cooccur_preknow <- t(gbi_preknowledge) %*% gbi_preknowledge

diag(cooccur_preknow) <- 0

library(reshape2)

edges_preknow <- melt(cooccur_preknow,
              varnames = c("ID1", "ID2"),
              value.name = "cooccurrence")
edges_preknow$ID1 <- as.character(edges_preknow$ID1)
edges_preknow$ID2 <- as.character(edges_preknow$ID2)
edges_preknow$dyad_ID <- ifelse(edges_preknow$ID1 < edges_preknow$ID2, paste(edges_preknow$ID1, edges_preknow$ID2, sep = "_"), paste(edges_preknow$ID2, edges_preknow$ID1, sep = "_"))
edges_preknow <- subset(edges_preknow, !edges_preknow$ID1 == edges_preknow$ID2)
edges_preknow <- edges_preknow[!duplicated(edges_preknow$dyad_ID), ]
edges_preknow$phase <- "preknowledge"

cooccur_postknow <- t(gbi_postknowledge) %*% gbi_postknowledge

diag(cooccur_postknow) <- 0

library(reshape2)

edges_postknow <- melt(cooccur_postknow,
                      varnames = c("ID1", "ID2"),
                      value.name = "cooccurrence")
edges_postknow$ID1 <- as.character(edges_postknow$ID1)
edges_postknow$ID2 <- as.character(edges_postknow$ID2)
edges_postknow$dyad_ID <- ifelse(edges_postknow$ID1 < edges_postknow$ID2, paste(edges_postknow$ID1, edges_postknow$ID2, sep = "_"), paste(edges_postknow$ID2, edges_postknow$ID1, sep = "_"))
edges_postknow <- subset(edges_postknow, !edges_postknow$ID1 == edges_postknow$ID2)
edges_postknow <- edges_postknow[!duplicated(edges_postknow$dyad_ID), ]
edges_postknow$phase <- "postknowledge"

edges_know <- bind_rows(edges_preknow, edges_postknow)

edges_know_long <- edges_know

edges_know_long <-
  edges_know %>%
  complete(ID1, ID2, phase,
           fill = list(cooccurrence = 0))

edges_know_long$dyad_ID <- ifelse(edges_know_long$ID1 < edges_know_long$ID2, paste(edges_know_long$ID1, edges_know_long$ID2, sep = "_"), paste(edges_know_long$ID2, edges_know_long$ID1, sep = "_"))

edges_know_long <- subset(edges_know_long, !edges_know_long$ID1 == edges_know_long$ID2)

edges_know_long$ID1_phase <- paste(edges_know_long$ID1, edges_know_long$phase, sep = "_")
edges_know_long$ID2_phase <- paste(edges_know_long$ID2, edges_know_long$phase, sep = "_")
edges_know_long$dyad_ID_phase <- paste(edges_know_long$dyad_ID, edges_know_long$phase, sep = "_")
edges_know_long <- edges_know_long[!duplicated(edges_know_long$dyad_ID_phase), ]

events_per_indiv_preknow <- as.data.frame(colSums(gbi_preknowledge))
events_per_indiv_preknow$JID <- row.names(events_per_indiv_preknow)
events_per_indiv_preknow <- events_per_indiv_preknow %>%
  rename(events = `colSums(gbi_preknowledge)` )
events_per_indiv_preknow$phase <- "preknowledge"
events_per_indiv_preknow$JID_phase <- paste(events_per_indiv_preknow$JID, events_per_indiv_preknow$phase, sep = "_")

events_per_indiv_postknow <- as.data.frame(colSums(gbi_postknowledge))
events_per_indiv_postknow$JID <- row.names(events_per_indiv_postknow)
events_per_indiv_postknow <- events_per_indiv_postknow %>%
  rename(events = `colSums(gbi_postknowledge)` )
events_per_indiv_postknow$phase <- "postknowledge"
events_per_indiv_postknow$JID_phase <- paste(events_per_indiv_postknow$JID, events_per_indiv_postknow$phase, sep = "_")

events_per_indiv_know <- bind_rows(events_per_indiv_preknow, events_per_indiv_postknow)

edges_know_long$ID1_events <- events_per_indiv_know$events[match(edges_know_long$ID1_phase, events_per_indiv_know$JID_phase)]
edges_know_long$ID2_events <- events_per_indiv_know$events[match(edges_know_long$ID2_phase, events_per_indiv_know$JID_phase)]

edges_know_long$ID1_events[is.na(edges_know_long$ID1_events)] <- 0
edges_know_long$ID2_events[is.na(edges_know_long$ID2_events)] <- 0

edges_know_long$ID1_knowledge <- all_individuals_knowledge$knowledgeable[match(edges_know_long$ID1, all_individuals_knowledge$JID)]
edges_know_long$ID2_knowledge <- all_individuals_knowledge$knowledgeable[match(edges_know_long$ID2, all_individuals_knowledge$JID)]

edges_know_long$ID1_knowledge[is.na(edges_know_long$ID1_knowledge )] <- "naive"
edges_know_long$ID2_knowledge[is.na(edges_know_long$ID2_knowledge )] <- "naive"

edges_know_long$knowledge_comb <- paste(edges_know_long$ID1_knowledge, edges_know_long$ID2_knowledge, sep = "_")

edges_know_long$knowledge_comb <- ifelse(!edges_know_long$knowledge_comb == "naive_knowledgeable", edges_know_long$knowledge_comb, "knowledgeable_naive")

edges_know_long <- edges_know_long %>%
  group_by(dyad_ID) %>%
  mutate(cooccurrence_sum = sum(cooccurrence, na.rm = TRUE))

edges_know_long <- edges_know_long %>%
  group_by(dyad_ID) %>%
  add_count()

edges_know_long_sub <- subset(edges_know_long, edges_know_long$cooccurrence_sum > 0)
edges_know_long_sub <- subset(edges_know_long_sub, edges_know_long_sub$n > 1)
edges_know_long_sub <- subset(edges_know_long_sub, edges_know_long_sub$ID1_events > 4)
edges_know_long_sub <- subset(edges_know_long_sub, edges_know_long_sub$ID2_events > 4)

edges_know_long_sub <- edges_know_long_sub %>%
  group_by(dyad_ID) %>%
  add_count()

edges_know_long_sub <- subset(edges_know_long_sub, edges_know_long_sub$nn > 1)

boxplot(edges_know_long_sub$cooccurrence ~ edges_know_long_sub$knowledge_comb * edges_know_long_sub$phase)

tapply(edges_know_long_sub$cooccurrence[edges_know_long_sub$knowledge_comb =="naive_naive"], edges_know_long_sub$phase[edges_know_long_sub$knowledge_comb =="naive_naive"], mean)
tapply(edges_know_long_sub$cooccurrence[edges_know_long_sub$knowledge_comb =="knowledgeable_naive"], edges_know_long_sub$phase[edges_know_long_sub$knowledge_comb =="knowledgeable_naive"], mean)
tapply(edges_know_long_sub$cooccurrence[edges_know_long_sub$knowledge_comb =="knowledgeable_knowledgeable"], edges_know_long_sub$phase[edges_know_long_sub$knowledge_comb =="knowledgeable_knowledgeable"], mean)
tapply(edges_know_long_sub$cooccurrence, edges_know_long_sub$knowledge_comb, mean)
tapply(edges_know_long_sub$ID1_events, edges_know_long_sub$knowledge_comb, mean)
tapply(edges_know_long_sub$ID2_events, edges_know_long_sub$knowledge_comb, mean)

cooccur_glmm1 <- glmmTMB(cooccurrence ~ phase * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edges_know_long_sub, family = poisson(link = "log"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

library(interactions)
cat_plot(cooccur_glmm1, phase, knowledge_comb)

emmeans(cooccur_glmm1, ~ phase * knowledge_comb, type = "response") |> pairs()



ids <- all_individuals_knowledge_sub2$JID[1:154]
ids <- as.character(ids)

edge_list_knowledge <- t(combn(ids, 2))
edge_list_knowledge <- as.data.frame(edge_list)

colnames(edge_list) <- c("ID1", "ID2")

edge_list_preknow <- edge_list_knowledge
edge_list_preknow$phase <- "preknowledge"

edge_list_postknow <- edge_list_knowledge
edge_list_postknow$phase <- "postknowledge"

edge_list_know <- rbind(edge_list_preknow, edge_list_postknow)

edge_list_know$dyad_ID <- ifelse(edge_list_know$ID1 < edge_list_know$ID2, paste(edge_list_know$ID1, edge_list_know$ID2, sep = "_"), paste (edge_list_know$ID2, edge_list_know$ID1, sep = "_"))
edge_list_know$dyad_ID_phase <- paste(edge_list_know$dyad_ID, edge_list_know$phase, sep = "_")
edge_list_know$ID1_phase <- paste(edge_list_know$ID1, edge_list_know$phase, sep = "_")
edge_list_know$ID2_phase <- paste(edge_list_know$ID2, edge_list_know$phase, sep = "_")

edge_list_know$edge <- edges_know_long_sub$cooccurrence[match(edge_list_know$dyad_ID_phase, edges_know_long_sub$dyad_ID_phase)]
edge_list_know$edge <- ifelse(!is.na(edge_list_know$edge), edge_list_know$edge, 0)

edge_list_know$ID1_events <- edges_know_long_sub$ID1_events[match(edge_list_know$ID1_phase, edges_know_long_sub$ID1_phase)]
edge_list_know$ID2_events <- edges_know_long_sub$ID2_events[match(edge_list_know$ID2_phase, edges_know_long_sub$ID2_phase)]
edge_list_know$ID1_events <- ifelse(!is.na(edge_list_know$ID1_events), edge_list_know$ID1_events, 0)
edge_list_know$ID2_events <- ifelse(!is.na(edge_list_know$ID2_events), edge_list_know$ID2_events, 0)

edge_list_know$ID1_knowledge <- all_individuals_knowledge_sub$knowledge[match(edge_list_know$ID1, all_individuals_knowledge_sub$JID)]
edge_list_know$ID2_knowledge <- all_individuals_knowledge_sub$knowledge[match(edge_list_know$ID2, all_individuals_knowledge_sub$JID)]
edge_list_know$knowledge_comb <- paste(edge_list_know$ID1_knowledge, edge_list_know$ID2_knowledge, sep = "_")
edge_list_know$knowledge_comb <- ifelse(!edge_list_know$knowledge_comb == "naive_knowledgeable", edge_list_know$knowledge_comb, "knowledgeable_naive")

table(edge_list_know$knowledge_comb)
boxplot(edge_list_know$edge ~ edge_list_know$knowledge_comb * edge_list_know$phase)

dual_events_knowledge_dyads_sub$dyad_ID <- ifelse(dual_events_knowledge_dyads_sub$JID1 < dual_events_knowledge_dyads_sub$JID2, paste(dual_events_knowledge_dyads_sub$JID1, dual_events_knowledge_dyads_sub$JID2, sep = "_"), paste(dual_events_knowledge_dyads_sub$JID2, dual_events_knowledge_dyads_sub$JID1, sep = "_"))
edge_list_know$HKLK_dyadic_events <- dual_events_knowledge_dyads_sub$HKLK_no_dual_events_mean[match(edge_list_know$dyad_ID, dual_events_knowledge_dyads_sub$dyad_ID)]
edge_list_know$HKLK_dyadic_events <- ifelse(!is.na(edge_list_know$HKLK_dyadic_events), edge_list_know$HKLK_dyadic_events, 0)  
  
cooccur_glmm1 <- glmmTMB(edge ~ phase * knowledge_comb + HKLK_dyadic_events * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edge_list_know, family = poisson(link = "log"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

library(interactions)
cat_plot(cooccur_glmm1, phase, knowledge_comb)

emmeans(cooccur_glmm1, ~ phase * knowledge_comb, type = "response") |> pairs()

edge_list_preknow_sub <- subset(edge_list_know, edge_list_know$phase == "preknowledge")
edge_list_preknow_sub <- subset(edge_list_preknow_sub, edge_list_preknow_sub$edge > 0)

edge_list_postknow_sub <- subset(edge_list_know, edge_list_know$phase == "postknowledge")

edge_list_know_sub <- rbind(edge_list_preknow_sub, edge_list_postknow_sub)

edge_list_know_sub <- edge_list_know_sub %>%
  group_by(dyad_ID) %>%
  add_count()

edge_list_know_sub <- subset(edge_list_know_sub, edge_list_know_sub$n > 1)

cooccur_glmm1 <- glmmTMB(edge ~ phase * knowledge_comb + HKLK_dyadic_events * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edge_list_know_sub, family = poisson(link = "log"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

cat_plot(cooccur_glmm1, phase, knowledge_comb)

emmeans(cooccur_glmm1, ~ phase * knowledge_comb, type = "response") |> pairs()

edge_list_know_sub <- subset(edge_list_know_sub, edge_list_know_sub$phase == "postknowledge")

edge_list_know_sub$edge_binary <- ifelse(edge_list_know_sub$edge > 0, 1, 0)

cooccur_glmm1 <- glmmTMB(edge_binary ~ HKLK_dyadic_events * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edge_list_know_sub, family = binomial(link = "logit"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

emmeans(cooccur_glmm1, ~ knowledge_comb, type = "response") |> pairs()


#Edge list from individuals that visited experimental locations 

all_individuals_knowledge_sub$JID <- as.character(all_individuals_knowledge_sub$JID)
ids <- all_individuals_knowledge_sub$JID
ids <- as.character(ids)

edge_list_knowledge <- t(combn(ids, 2))
edge_list_knowledge <- as.data.frame(edge_list_knowledge)

colnames(edge_list_knowledge) <- c("ID1", "ID2")

edge_list_knowledge_pre <- edge_list_knowledge
edge_list_knowledge_pre$phase <- "preknowledge"

edge_list_knowledge_post <- edge_list_knowledge
edge_list_knowledge_post$phase <- "postknowledge"

edge_list_knowledge <- rbind(edge_list_knowledge_pre, edge_list_knowledge_post)

edge_list_knowledge$dyad_ID <- ifelse(edge_list_knowledge$ID1 < edge_list_knowledge$ID2, paste(edge_list_knowledge$ID1, edge_list_knowledge$ID2, sep = "_"), paste (edge_list_knowledge$ID2, edge_list_knowledge$ID1, sep = "_"))
edge_list_knowledge$dyad_ID_phase <- paste(edge_list_knowledge$dyad_ID, edge_list_knowledge$phase, sep = "_")
edge_list_knowledge$ID1_phase <- paste(edge_list_knowledge$ID1, edge_list_knowledge$phase, sep = "_")
edge_list_knowledge$ID2_phase <- paste(edge_list_knowledge$ID2, edge_list_knowledge$phase, sep = "_")

edge_list_knowledge$edge <- edges_know_long_sub$cooccurrence[match(edge_list_knowledge$dyad_ID_phase, edges_know_long_sub$dyad_ID_phase)]
edge_list_knowledge$edge <- ifelse(!is.na(edge_list_knowledge$edge), edge_list_knowledge$edge, 0)

edge_list_knowledge$ID1_events <- edges_know_long_sub$ID1_events[match(edge_list_knowledge$ID1_phase, edges_know_long_sub$ID1_phase)]
edge_list_knowledge$ID2_events <- edges_know_long_sub$ID2_events[match(edge_list_knowledge$ID2_phase, edges_know_long_sub$ID2_phase)]
edge_list_knowledge$ID1_events <- ifelse(!is.na(edge_list_knowledge$ID1_events), edge_list_knowledge$ID1_events, 0)
edge_list_knowledge$ID2_events <- ifelse(!is.na(edge_list_knowledge$ID2_events), edge_list_knowledge$ID2_events, 0)

edge_list_knowledge$ID1_knowledge <- all_individuals_knowledge_sub$knowledge[match(edge_list_knowledge$ID1, all_individuals_knowledge_sub$JID)]
edge_list_knowledge$ID2_knowledge <- all_individuals_knowledge_sub$knowledge[match(edge_list_knowledge$ID2, all_individuals_knowledge_sub$JID)]
edge_list_knowledge$knowledge_comb <- paste(edge_list_knowledge$ID1_knowledge, edge_list_knowledge$ID2_knowledge, sep = "_")
edge_list_knowledge$knowledge_comb <- ifelse(!edge_list_knowledge$knowledge_comb == "naive_knowledgeable", edge_list_knowledge$knowledge_comb, "knowledgeable_naive")

table(edge_list_knowledge$knowledge_comb)
boxplot(edge_list_knowledge$edge ~ edge_list_knowledge$knowledge_comb * edge_list_knowledge$phase)

dual_events_knowledge_dyads_sub$dyad_ID <- ifelse(dual_events_knowledge_dyads_sub$JID1 < dual_events_knowledge_dyads_sub$JID2, paste(dual_events_knowledge_dyads_sub$JID1, dual_events_knowledge_dyads_sub$JID2, sep = "_"), paste(dual_events_knowledge_dyads_sub$JID2, dual_events_knowledge_dyads_sub$JID1, sep = "_"))
edge_list_knowledge$HKLK_dyadic_events <- dual_events_knowledge_dyads_sub$HKLK_no_dual_events_mean[match(edge_list_knowledge$dyad_ID, dual_events_knowledge_dyads_sub$dyad_ID)]
edge_list_knowledge$HKLK_dyadic_events <- ifelse(!is.na(edge_list_knowledge$HKLK_dyadic_events), edge_list_knowledge$HKLK_dyadic_events, 0)  

edge_list_knowledge_pre_sub <- subset(edge_list_knowledge, edge_list_knowledge$phase == "preknowledge")
edge_list_knowledge_pre_sub <- subset(edge_list_knowledge_pre_sub, edge_list_knowledge_pre_sub$edge > 0)

edge_list_knowledge_post_sub <- subset(edge_list_knowledge, edge_list_knowledge$phase == "postknowledge")

edge_list_knowledge_sub <- rbind(edge_list_knowledge_pre_sub, edge_list_knowledge_post_sub)

edge_list_knowledge_sub <- edge_list_knowledge_sub %>%
  group_by(dyad_ID) %>%
  add_count()

edge_list_knowledge_sub <- subset(edge_list_knowledge_sub, edge_list_knowledge_sub$n > 1)

edge_list_knowledge_sub <- subset(edge_list_knowledge, edge_list_knowledge$HKLK_dyadic_events > 0)

cooccur_glmm1 <- glmmTMB(edge ~ phase * knowledge_comb + HKLK_dyadic_events * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edge_list_knowledge, family = poisson(link = "log"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

library(interactions)
cat_plot(cooccur_glmm1, phase, knowledge_comb)
interact_plot(cooccur_glmm1, HKLK_dyadic_events, knowledge_comb)

emmeans(cooccur_glmm1, ~ phase * knowledge_comb, type = "response") |> pairs()

edge_list_knowledge_sub <- subset(edge_list_knowledge_sub, edge_list_knowledge_sub$phase == "postknowledge")
edge_list_knowledge_sub$edge_binary <- ifelse(edge_list_knowledge_sub$edge > 0, 1, 0)

cooccur_glmm1 <- glmmTMB(edge_binary ~ HKLK_dyadic_events * knowledge_comb + ID1_events + ID2_events + (1|ID1) + (1|ID2) + (1|dyad_ID), data = edge_list_know_sub, family = binomial(link = "logit"))
summary(cooccur_glmm1)
Anova(cooccur_glmm1)

emmeans(cooccur_glmm1, ~ knowledge_comb, type = "response") |> pairs()

emmeans(cooccur_glmm1, ~ knowledge_comb, type = "response") |> pairs()


# (8) FIGURES ----


# (9) ARCHIVE ----

all_individuals_knowledge_sub6 <- pivot_longer(all_individuals_knowledge_sub, cols = c("visit_number_preknowledge", "visit_number_postknowledge"), names_to = "phase_knowledge", values_to = "visit_number_passive")
all_individuals_knowledge_sub6$phase_knowledge <- ifelse(all_individuals_knowledge_sub6$phase_knowledge == "visit_number_preknowledge", "pre_knowledge", "post_knowledge")
all_individuals_knowledge_sub6$JID_phase <- paste(all_individuals_knowledge_sub6$JID, all_individuals_knowledge_sub4$knowledge_phase, sep = " ") 

all_individuals_knowledge_sub4$visit_number_passive <- all_individuals_knowledge_sub6$visit_number_passive[match(all_individuals_knowledge_sub6$JID_phase, all_individuals_knowledge_sub4$JID_phase)]

#Number of visits at HK feeders 

#Model just with prior
default_prior()

visit_no_HK_brm1_prior <- brm(visit_number_HK ~ knowledge, 
                              data = all_individuals_knowledge_sub2, 
                              family = negbinomial(link = "log"),
                              prior = c(
                                set_prior("normal(log(90), 1)", class = "Intercept"),  
                                set_prior("normal(0, 0.5)", class = "b")),              
                              sample_prior = "only")

#Prior predictive checks 
pp_check(visit_no_HK_brm1_prior, ndraws = 1000)

#Model
visit_no_HK_brm1 <- brm(visit_number_HK ~ knowledge, 
                        data = all_individuals_knowledge_sub2[!is.na(all_individuals_knowledge_sub2$visit_number_HK),], 
                        family = negbinomial(link = "log"),
                        prior = c(
                          set_prior("normal(log(90), 1)", class = "Intercept"),  
                          set_prior("normal(0, 0.5)", class = "b")),              
                        
)

summary(visit_no_HK_brm1)

check_collinearity(visit_no_HK_brm1)

#Posterior predictive checks
pp_check(visit_no_HK_brm1, type = "dens_overlay", ndraws = 1000)

pp_check(visit_no_HK_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(visit_no_HK_brm1, type = "bars", ndraws = 1000)
pp_check(visit_no_HK_brm1, type = "hist", ndraws = 1000)

y <- model.frame(visit_no_HK_brm1)[[1]]
yrep <- posterior_predict(visit_no_HK_brm1)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(visit_no_HK_brm1)
launch_shinystan(visit_no_HK_brm1)

#Posterior distribution
as_draws_df(visit_no_HK_brm1)
mcmc_areas(visit_no_HK_brm1)
mcmc_intervals(visit_no_HK_brm1)

#Evaluation and interpretation
loo(visit_no_HK_brm1, moment_match = TRUE)
fitted(visit_no_HK_brm1, scale = "response")
conditional_effects(visit_no_HK_brm1)
bayes_R2(visit_no_HK_brm1)

#Sensitivity analysis and prior checks 
prior_summary(visit_no_HK_brm1)

#Extract predictions
fitted(visit_no_HK_brm1)


#Number of visits at LK feeders 

mean(all_individuals_knowledge_sub2$visit_number_LK, na.rm = TRUE)

#Model just with prior
default_prior()

visit_no_LK_brm1_prior <- brm(visit_number_LK ~ knowledge, 
                              data = all_individuals_knowledge_sub2, 
                              family = negbinomial(link = "log"),
                              prior = c(
                                set_prior("normal(log(70), 1)", class = "Intercept"),  
                                set_prior("normal(0, 0.5)", class = "b")),              
                              sample_prior = "only")

#Prior predictive checks 
pp_check(visit_no_LK_brm1_prior, ndraws = 1000)

#Model
visit_no_LK_brm1 <- brm(visit_number_LK ~ knowledge, 
                        data = all_individuals_knowledge_sub2[!is.na(all_individuals_knowledge_sub2$visit_number_LK) & all_individuals_knowledge_sub2$phase == "pre_knowledge",], 
                        family = negbinomial(link = "log"),
                        prior = c(
                          set_prior("normal(log(70), 1)", class = "Intercept"),  
                          set_prior("normal(0, 0.5)", class = "b")), 
                        control = list(max_treedepth = 15)
                        
)

summary(visit_no_LK_brm1)

check_collinearity(visit_no_LK_brm1)

#Posterior predictive checks
pp_check(visit_no_LK_brm1, type = "dens_overlay", ndraws = 1000)

pp_check(visit_no_LK_brm1, type = "stat", stat = "var", ndraws = 1000)
pp_check(visit_no_LK_brm1, type = "bars", ndraws = 1000)
pp_check(visit_no_LK_brm1, type = "hist", ndraws = 1000)

y <- model.frame(visit_no_LK_brm1)[[1]]
yrep <- posterior_predict(visit_no_LK_brm1)
# Mean variance ratio across posterior draws:
dispersion <- apply(yrep, 1, function(x) var(x - y))
mean(dispersion)

#Plot model
plot(visit_no_LK_brm1)
launch_shinystan(visit_no_LK_brm1)

#Posterior distribution
as_draws_df(visit_no_LK_brm1)
mcmc_areas(visit_no_LK_brm1)
mcmc_intervals(visit_no_HK_brm1)

#Evaluation and interpretation
loo(visit_no_HK_brm1, moment_match = TRUE)
fitted(visit_no_HK_brm1, scale = "response")
conditional_effects(visit_no_HK_brm1)
conditional_effects(visit_no_HK_brm1, effects = "knowledge:phase")

bayes_R2(visit_no_HK_brm1)

#Sensitivity analysis and prior checks 
prior_summary(visit_no_HK_brm1)

#Extract predictions
fitted(visit_no_HK_brm1)
#Changes in feeder visits - learning

#Visits at passive feeders
aggregate(all_individuals_knowledge_sub2$visit_number_passive, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub2$visit_number_passive, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "sd", na.rm = TRUE)

knowledge_passive_brm0 <- brm(visit_number_passive ~ knowledge + (1|JID), data = all_individuals_knowledge_sub2, family = poisson())
knowledge_passive_brm0 <- brm(visit_number_passive ~ knowledge_phase * knowledge + (1|JID), data = all_individuals_knowledge_sub2, family = poisson())
summary(knowledge_passive_brm0)
plot(knowledge_passive_brm0)

interaction.plot(x.factor = all_individuals_knowledge_sub2$phase, 
                 trace.factor = all_individuals_knowledge_sub2$knowledge,  
                 response = all_individuals_knowledge_sub2$visit_number_passive, fun = mean)

colours <- c("#FFFFFF", "#CCCCCC")

passive_visit_boxplot <- ggplot(all_individuals_knowledge_sub2, aes(x = phase_knowledge, y = visit_number_passive, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Knowledge", y = "Visit number") +
  theme_base(base_size = 14)+
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  
passive_visit_boxplot


#Visits at low quality feeders

aggregate(all_individuals_knowledge_sub2$visit_number_LK, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub2$visit_number_LK, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "sd", na.rm = TRUE)

knowledge_LK_brm0 <- brm(visit_number_LK ~ knowledge, data = all_individuals_knowledge_sub2, family = poisson())
summary(knowledge_LK_brm0)
plot(knowledge_LK_brm0)

#Visits at high quality feeders

aggregate(all_individuals_knowledge_sub2$visit_number_HK, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub2$visit_number_HK, by = list(all_individuals_knowledge_sub2$knowledge), FUN = "sd", na.rm = TRUE)

knowledge_HK_brm0 <- brm(visit_number_HK ~ knowledge, data = all_individuals_knowledge_sub2, family = poisson())
summary(knowledge_HK_brm0)
plot(knowledge_HK_brm0)

aggregate(visit_data_knowledge_HK$visit_duration, by = list(visit_data_knowledge_HK$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(visit_data_knowledge_HK$visit_duration, by = list(visit_data_knowledge_HK$knowledge), FUN = "sd", na.rm = TRUE)

aggregate(visit_data_knowledge_LK$visit_duration, by = list(visit_data_knowledge_LK$knowledge), FUN = "mean", na.rm = TRUE)
aggregate(visit_data_knowledge_LK$visit_duration, by = list(visit_data_knowledge_LK$knowledge), FUN = "sd", na.rm = TRUE)

HK_visit_boxplot <- ggplot(visit_data_knowledge_HK, aes(x = knowledge, y = visit_duration, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Knowledge", y = "Visit duration high quality feeder") +
  theme_classic(base_size = 14)+
  ylim(0, 650) +
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  
HK_visit_boxplot

LK_visit_boxplot <- ggplot(visit_data_knowledge_LK, aes(x = knowledge, y = visit_duration, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  scale_fill_manual(values = colours) +
  labs(x="Knowledge", y = "Visit duration low quality feeder") +
  theme_classic(base_size = 14)+
  ylim(0, 650) +
  geom_jitter(alpha = 0.3, shape=16, position=position_jitter(0.1), size = 1.5) + theme(legend.position="none") +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  
LK_visit_boxplot

#Visit duration at HK feeders
knowledge_HK_brm0 <- brm(visit_duration ~ knowledge + (1|JID), data = visit_data_knowledge_HK, family = gaussian())
knowledge_HK_glmm0 <- glmmTMB(visit_duration + 1 ~ knowledge + (1|JID), data = visit_data_knowledge_HK, family = Gamma(link = "log"))

summary(knowledge_HK_brm0)
summary(knowledge_HK_glmm0)

plot(knowledge_HK_brm0)

all_individuals_knowledge_sub2

all_individuals_knowledge_sub3 <- subset(all_individuals_knowledge_sub2, all_individuals_knowledge_sub2$knowledge == "knowledgeable")

all_individuals_knowledge_sub2$phase_knowledge <- paste(all_individuals_knowledge_sub2$phase, all_individuals_knowledge_sub2$knowledge, sep = "_")

#reorder phase_knowledge factor to show pre first and post second
all_individuals_knowledge_sub2$phase_knowledge <- factor(all_individuals_knowledge_sub2$phase_knowledge, levels=c("pre_knowledge_knowledgeable", "pre_knowledge_naive", "post_knowledge_knowledgeable", "post_knowledge_naive"))

aggregate(all_individuals_knowledge_sub2$degree, by = list(all_individuals_knowledge_sub2$phase_knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub2$degree, by = list(all_individuals_knowledge_sub2$phase_knowledge), FUN = "sd", na.rm = TRUE)

aggregate(all_individuals_knowledge_sub2$strength, by = list(all_individuals_knowledge_sub2$phase_knowledge), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_knowledge_sub2$strength, by = list(all_individuals_knowledge_sub2$phase_knowledge), FUN = "sd", na.rm = TRUE)

knowledge_degree_brm0 <- brm(degree ~ phase * knowledge + visit_number_passive + (1|JID), data = all_individuals_knowledge_sub2, family = poisson())

summary(knowledge_degree_brm0)

knowledge_strength_brm0 <- brm(strength ~ phase * knowledge + visit_number_passive + (1|JID), data = all_individuals_knowledge_sub2, family = gaussian())

summary(knowledge_strength_brm0)

interaction.plot(x.factor = all_individuals_knowledge_sub2$phase, 
                 trace.factor = all_individuals_knowledge_sub2$knowledge,  
                 response = all_individuals_knowledge_sub2$strength, fun = mean)

interaction.plot(x.factor = all_individuals_knowledge_sub2$phase, 
                 trace.factor = all_individuals_knowledge_sub2$knowledge,  
                 response = all_individuals_knowledge_sub2$degree, fun = mean)



strength_boxplot <- ggplot(all_individuals_knowledge_sub2, aes(x = phase_knowledge, y = strength, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  geom_point(alpha = 0.3) +
  scale_fill_manual(values = colours) +
  labs(x="Phase/knowledge", y = "Strength") +
  scale_x_discrete(labels = c("knowledgeable (pre)","naive (pre)", "knowledgeable (post)", "naive (post)")) +
  theme_classic(base_size = 14) +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  +
  theme(legend.position= "none") 
strength_boxplot

degree_boxplot <- ggplot(all_individuals_knowledge_sub2, aes(x = phase_knowledge, y = degree, fill = knowledge)) + 
  geom_violin(width = 0.75) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +
  geom_point(alpha = 0.3) +
  scale_fill_manual(values = colours) +
  labs(x="Phase/knowledge", y = "Degree") +
  scale_x_discrete(labels = c("knowledgeable (pre)","naive (pre)", "knowledgeable (post)", "naive (post)")) +
  theme_classic(base_size = 14) +
  stat_summary(fun=mean, geom="point", shape=8, size=3, col = "coral")  +
  theme(legend.position= "none") 
degree_boxplot


posterior <- as.matrix(displace_brm2)

mcmc_areas(posterior,
           pars = c("b_knowledge_combknowledgeable_naive:qualitylow", "b_knowledge_combnaive_knowledge:qualitylow"),
           prob = 0.8) 

cond <- conditional_effects(displace_brm1, effects = "knowledge:quality", re_formula = NA)
plot(cond)


#Model 1: knowledge * quality 

#Model just with prior
default_prior()

displace_brm1_prior <- brm(displaced ~ knowledge * quality + (1|JID), 
                           data = visit_data_knowledge_post, 
                           family = bernoulli(link = "logit"),
                           sample_prior = "only"
)

#Prior predictive checks 
pp_check(displace_brm1_prior, ndraws = 1000)

displace_brm1 <- brm(displaced ~ knowledge * quality + (1 | mm(JID, next_JID)) + (1|feeder), 
                     data = visit_data_knowledge_post, 
                     family = bernoulli(link = "logit"),
                     prior = c(
                       set_prior("normal(-1.75, 0.5)", class = "Intercept"),  
                       set_prior("normal(0, 0.5)", class = "b")),
                     control = list(adapt_delta = 0.99))

summary(displace_brm1)

check_collinearity(displace_brm1)

#Posterior predictive checks
pp_check(displace_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(displace_brm1)
launch_shinystan(displace_brm1)

#Posterior distribution
as_draws_df(displace_brm1)
mcmc_areas(displace_brm1)
mcmc_intervals(displace_brm1)

#Evaluation and interpretation
loo(displace_brm1, moment_match = TRUE)
fitted(displace_brm1, scale = "response")
conditional_effects(displace_brm1)
conditional_effects(displace_brm1,effects = "knowledge:quality")

bayes_R2(displace_brm1)
emmeans(displace_brm1, ~ knowledge * quality, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(displace_brm1)

#Extract predictions
fitted(displace_brm1)

cond <- conditional_effects(displace_brm1, effects = "knowledge:quality")
p <- plot(cond, points=F)
displace_plot1 <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = knowledge, y = displaced), 
    data = visit_data_knowledge_post, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Knowledge", y = "Probability of displacement")
displace_plot1

#Does pattern hold with binomial data and not bernoulli data
displace_brm1.2 <- brm(n_displaced | trials(n_visits) ~ knowledge, 
                       data = displacements_knowledge_post_summary, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "Intercept"),  
                         set_prior("normal(0, 0.5)", class = "b")),
                       control = list(adapt_delta = 0.99))


summary(displace_brm1.2)

displace_brm1.3 <- brm(n_displaced | trials(n_visits) ~ knowledge, 
                       data = displacement_knowledge_LK_summary, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "Intercept"),  
                         set_prior("normal(0, 0.5)", class = "b")),
                       control = list(adapt_delta = 0.99))

summary(displace_brm1.3)

#Model edge change ~ knowledge_comb * HKLK_no_dual_events_mean

#Only pre-existing edges

#Model just with prior
default_prior()

edge_brm1_prior <- brm(edge_diff_resid ~ knowledge_comb2 * HKLK_no_dual_events_mean, 
                       data = dyads_prepost_sub, 
                       family = gaussian(),
                       prior = c(
                         set_prior("normal(0, 0.15)", class = "Intercept"),  
                         set_prior("normal(0, 0.15)", class = "b")),
                       sample_prior = "only")

#Prior predictive checks 
pp_check(edge_brm1_prior, ndraws = 1000)

edge_brm1 <- brm(edge_diff_resid ~ knowledge_comb2 * HKLK_no_dual_events_mean, 
                 data = dyads_prepost_sub, 
                 family = gaussian(),
                 prior = c(
                   set_prior("normal(0, 0.15)", class = "Intercept"),  
                   set_prior("normal(0, 0.15)", class = "b")))
summary(edge_brm1)

check_collinearity(edge_brm1)

#Posterior predictive checks
pp_check(edge_brm1, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(edge_brm1)
launch_shinystan(edge_brm1)

#Posterior distribution
as_draws_df(edge_brm1)
mcmc_areas(edge_brm1)
mcmc_intervals(edge_brm1)

#Evaluation and interpretation
loo(edge_brm1, moment_match = TRUE)
fitted(edge_brm1, scale = "response")
conditional_effects(edge_brm1)
conditional_effects(edge_brm1,effects = "HKLK_no_dual_events_mean:knowledge_comb2")

bayes_R2(edge_brm1)
emmeans(edge_brm1, ~ knowledge_comb2 * HKLK_no_dual_events_mean, type = "response") |> pairs()

#Sensitivity analysis and prior checks 
prior_summary(edge_brm1)

#Extract predictions
fitted(edge_brm1)






ce <- conditional_effects(displace_brm2, effects = "knowledge_comb:quality", re_formula = NA)
newdat <- as.data.frame(ce[["knowledge_comb:quality"]]) %>%
  dplyr::select(knowledge_comb, quality)
posterior <- posterior_epred(displace_brm2, newdata = newdat, re_formula = NA)

# Convert to long format
posterior_long <- as.data.frame(posterior) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "draw", values_to = "pred") %>%
  mutate(row = as.numeric(gsub("V", "", draw))) %>%
  left_join(newdat %>% mutate(row = row_number()), by = "row")

posterior_long_LK <- subset(posterior_long, quality == "low")
posterior_long_HK <- subset(posterior_long, quality == "high")

displace_plot2 <- ggplot(posterior_long, aes(x = knowledge_comb, y = pred, fill = quality)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values =c("black", "grey")) +
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8)) +
  labs(
    x = "Knowledge combination",
    y = "Predicted probability of displacement",
    fill = "Feeding station quality") +
  scale_x_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displace_plot2

displace_plot_LK <- ggplot(posterior_long_LK, aes(x = knowledge_comb, y = pred, fill = quality)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values =c("grey")) +
  ylim(c(0, 0.5)) +
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8)) +
  labs(
    x = "Knowledge combination",
    y = "Predicted probability of displacement",
    fill = "Feeding station quality") +
  scale_x_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displace_plot_LK

displace_plot_HK <- ggplot(posterior_long_HK, aes(x = knowledge_comb, y = pred, fill = quality)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values =c("black")) +
  ylim(c(0, 0.5)) +
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.8)) +
  labs(
    x = "Knowledge combination",
    y = "Predicted probability of displacement",
    fill = "Feeding station quality") +
  scale_x_discrete(labels=c('informed informed', 'informed naive', 'naive informed', 'naive naive'))
displace_plot_HK

ggarrange(displace_plot_HK, displace_plot_LK, widths = 1)

ggplot(posterior_long, aes(x = knowledge_comb, y = pred, fill = quality)) +
  geom_violin(
    alpha = 0.6,
    position = position_dodge(width = 0.7),
    width = 0.6,
    trim = TRUE,
    color = NA
  ) +
  # Posterior mean + interval markers
  stat_summary(fun = median, geom = "point",
               position = position_dodge(width = 0.7),
               size = 2.5, shape = 21, fill = "white", color = "black") +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.8),
               geom = "linerange",
               position = position_dodge(width = 0.7),
               linewidth = 1) +
  labs(
    x = "Factor 1",
    y = "Predicted probability",
    fill = "Factor 2",
    title = "Posterior predicted probabilities",
    subtitle = "Half-eye style inspired visualization"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 12, face = "bold")
  )


posterior_draw_knowledge_comb <- posterior_marginal %>%
  group_by(.draw, knowledge, quality) %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_draw_knowledge_comb_HK <- posterior_marginal %>%
  group_by(.draw, knowledge_comb) %>%
  filter(quality == "high") %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

posterior_draw_knowledge_comb_LK <- posterior_marginal %>%
  group_by(.draw, knowledge_comb) %>%
  filter(quality == "low") %>%
  summarise(
    epred = mean(.epred),
    .groups = "drop"
  )

pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  group_by(.draw, knowledge_comb, quality) %>%
  compare_levels(epred, by = knowledge_comb) %>%
  rename(
    contrast = knowledge_comb,
    diff = epred
  )

pairwise_contrasts <- posterior_draw_knowledge_comb %>%
  group_by(.draw, knowledge_comb, quality) %>%
  compare_levels(epred, by = knowledge_comb, quality)

pairwise_contrasts <- pairwise_contrasts %>%
  rename(
    contrast = knowledge_comb,
    diff = epred
  )
