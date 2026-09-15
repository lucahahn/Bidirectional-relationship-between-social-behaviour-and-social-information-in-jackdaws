
#PhD Ch 4 Social Bonds and Social Learning - Patch Discovery Experiment
#Author: Luca Hahn
#Last update: 15/12/2024

# (1) IMPORT DATA  ----

#Load packages 
install.packages("asnipe")
install.packages("car")
install.packages("carData")
install.packages("chisq.posthoc.test")
install.packages("ClusterR")
install.packages("corrplot")
install.packages("data.table")
install.packages("DHARMa")
install.packages("dplyr")
install.packages("emmeans")
install.packages("ggdist")
install.packages("ggplot2")
install.packages("glmmTMB")
install.packages("hms")
install.packages("igraph")
install.packages("lme4")
install.packages("lmerTest")
install.packages("lubridate")
install.packages("MASS")
install.packages("multcomp")
install.packages("RColorBrewer")
install.packages("reshape2")
install.packages("rptR")
install.packages("scales")
install.packages("stringi")
install.packages("stringr")
install.packages("svMisc")
install.packages("tidybayes")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("assortnet")

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
library(STbayes)
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

if (!require("devtools")) install.packages("devtools")
devtools::install_github("michaelchimento/STbayes", force = TRUE)
install.packages("STbayes")
install.packages("glue")
install.packages("posterior")
library(posterior)

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
RT_files_patch <- list.files("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Patches", pattern = "RT", recursive = TRUE)
RT_paths_patch <- paste("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Patches",RT_files_patch, sep = "/")

#Load saved life history csv file 
LH <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/LH20240419.csv", header = T, stringsAsFactors = F)
LH <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/LH_20241117.csv", header = T, stringsAsFactors = F)

LH$DATE <- strptime(LH$DATE,format= "%d/%m/%Y")
LH$DATE <- as.Date(LH$DATE, format = "%d/%m/%Y") # convert to date

LH_ring <- subset(LH, !LH$COMBINATION == "")

LH_tarsus <- subset(LH, !LH$TARSUS == "")

LH_tarsus <- LH_tarsus %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

#LH_body <- subset(LH, !LH$TARSUS == "")
LH_body <- subset(LH_body, !LH_body$WEIGHT == "")
LH_body$BODY_COND <- resid(lm(WEIGHT ~ TARSUS, data = LH_body))

LH_body <- LH_body %>%
  group_by(ID) %>% 
  arrange(desc(DATE)) %>% 
  slice(1:1)

#calculate mean tarsus per individual 

mean_tarsus <- as.data.frame(aggregate(LH$TARSUS, by = list(LH$ID), mean, na.rm = TRUE))
mean_tarsus$JID <- mean_tarsus$Group.1
mean_tarsus$tarsus <- mean_tarsus$x
mean_tarsus$Group.1 <- NULL
mean_tarsus$x <- NULL

#most recent body condition using average tarsus
recent_body_cond <- mean_tarsus
recent_body_cond$weight <- LH_body$WEIGHT[match(recent_body_cond$JID, LH_body$ID)]
recent_body_cond <- subset(recent_body_cond, !is.na(recent_body_cond$weight))
recent_body_cond$recent_body_cond <- resid(lm(weight ~ tarsus, data = recent_body_cond))

#Temperament data
temperament <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/temperament.csv", header = T, stringsAsFactors = F)

#Create sub-strings that contain feeder ID (e.g. "Y1.1") and date
day_arrays_patch <- unique(substr(RT_files_patch,6,16))

#Create empty list to place data into as we go
collapsed_list_patch <- list()

#Run through RT files in list 
for(i in 1:length(day_arrays_patch)) {
  progress(i, max.value = length(day_arrays_patch))
  
  day_array_list_patch <- list()
  
  for (j in 1:length(RT_files_patch[which(substr(RT_files_patch,6,16) == day_arrays_patch[i])])) {
    temp_day_file_patch <- (read.delim(RT_paths_patch[which(substr(RT_files_patch,6,16) == day_arrays_patch[i])][j], header = T, stringsAsFactors = F))[,1:13]
    temp_day_file_patch$feeder <- substr(RT_files_patch[which(substr(RT_files_patch,6,16) == day_arrays_patch[i])][j],6,9)
    day_array_list_patch[[j]] <- temp_day_file_patch
  }
  
  temp_RT_patch <- do.call(rbind, day_array_list_patch)
  
  temp_RT_patch$Time <- strptime(paste(temp_RT_patch$Date,stri_sub(temp_RT_patch$Hmsec/1024,2,5), sep = ""), "%Y-%m-%d %H:%M:%OS")  # Add in miliseconds (1024 in a second) and format time
  
  temp_RT_patch %>% filter(nchar(TagID_hex) == 10) -> temp_tags_patch  #remove times when no tag
  
  
  if(dim(temp_tags_patch)[1] >0){  
    
    visits_patch <- data.frame(Event = temp_tags_patch$Event, Start = temp_tags_patch$Time, End = temp_tags_patch$Time+(0.5*(temp_tags_patch$Reps -1)), tag = temp_tags_patch$TagID_hex, feeder = temp_tags_patch$feeder)  # adds reps to visit length (0.5 seconds for every extra detection as that was resampling speed)
    visits_patch <- arrange(visits_patch, Start)
    visits_patch <- arrange(visits_patch, feeder)
    within_errors_patch <- which(visits_patch$End < lag(visits_patch$End) & visits_patch$tag == lag(visits_patch$tag) & visits_patch$feeder == lag(visits_patch$feeder))
    if(length(within_errors_patch) > 0){
      visits_patch <- visits_patch[-which(visits_patch$End < lag(visits_patch$End) & visits_patch$tag == lag(visits_patch$tag) & visits_patch$feeder == lag(visits_patch$feeder)),]  ## Get rid of reads within bouts - almost always erroneous single reads that are repeated later in the datastream
    }
    visits_patch <- arrange(visits_patch, Start)
    visits_patch$count <- sapply(1:nrow(visits_patch),function(x)sum(visits_patch$tag[x]==visits_patch$tag[1:x]))  ## add individual visit counter - if two bouts don't have sequential counts then bird seen elsewhere in between
    visits_patch <- arrange(visits_patch, feeder)
    
    visits_patch$collapse <- "0"  # Temp column to tell me if this bout is to be collapsed
    
    visit_time_patch <- 10 ## How long between detections before a new bout is classed
    
    # Is the last visit ending within 'visit_time' seconds of this one starting, and with the same tag & in sequence?
    visits_patch$collapse[ which (visits_patch$Start-lag(visits_patch$End) < visit_time_patch & lead(visits_patch$Start) - visits_patch$End < visit_time_patch & visits_patch$feeder == lag(visits_patch$feeder) & visits_patch$feeder == lead(visits_patch$feeder) & visits_patch$tag == lag(visits_patch$tag) & visits_patch$tag == lead(visits_patch$tag) & (visits_patch$count-lag(visits_patch$count)) == 1 & (lead(visits_patch$count)-visits_patch$count) == 1 )] <- "1"  ## What about times they hop in between?!
    
    visits_patch$collapse[which(visits_patch$collapse == 0 & lag(visits_patch$collapse == 1))] <- "End"  ## Can work out the end based on the 0s and 1s
    visits_patch$collapse[which(visits_patch$collapse == 0 & lead(visits_patch$collapse == 1))] <- "Start" ## Can work out the start based on the 0s and 1s
    
    ## Mop up those that only have two potential detections (so no middle values to get assigned 1 above)
    visits_patch$collapse[which(visits_patch$collapse == 0 & lead(visits_patch$Start) - visits_patch$End < visit_time_patch & lead(visits_patch$Start) - visits_patch$End > -2 & visits_patch$tag == lead(visits_patch$tag) & lead(visits_patch$count) - visits_patch$count == 1)] <- "Start"
    visits_patch$collapse[which(visits_patch$collapse == 0 & visits_patch$Start - lag(visits_patch$End) < visit_time_patch & lag(visits_patch$feeder) == visits_patch$feeder & visits_patch$tag == lag(visits_patch$tag) & lag(visits_patch$count) - visits_patch$count == -1)] <- "End"
    
    # Add this file's data to the list
    collapsed_list_patch[[i]] <- data.frame(start = visits_patch$Start[which(visits_patch$collapse %in% c("Start","0"))], end = visits_patch$End[which(visits_patch$collapse %in% c("End","0"))], tag = visits_patch$tag[which(visits_patch$collapse %in% c("Start","0"))], event = visits_patch$Event[which(visits_patch$collapse %in% c("Start","0"))], feeder_patch =  visits_patch$feeder[which(visits_patch$collapse %in% c("Start","0"))], array_patch = rep(stri_sub(day_arrays[i], 8,11), length(which(visits_patch$collapse %in% c("Start","0")))))  }
}

#Turn the list into a data frame
visit_data_patch <- do.call(rbind,collapsed_list_patch)

#Find instances in which time is still NA
visit_data_NA_patch <- subset(visit_data_patch, is.na(visit_data_patch$start))

#Filter instances in which time is not NA
visit_data_patch <- subset(visit_data_patch, !is.na(visit_data_patch$start))

#Remove instances where JID = NA and test tags
visit_data_patch$JID <- LH$ID[match(visit_data_patch$tag,LH$RFID)]
visit_data_patch <- subset(visit_data_patch, JID != "NA")

# (2) ADD DATA ----

#Adding information about visit duration
visit_data_patch$interval <- interval(visit_data_patch$start,visit_data_patch$end)
visit_data_patch$visit_duration <- as.duration(visit_data_patch$interval)

#Adding information about site: X, Y, Z
visit_data_patch[substr(visit_data_patch$feeder_patch,1,2)=="YY","site"]<-"YY"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,2)=="YZ","site"]<-"YZ"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,2)=="ZZ","site"]<-"ZZ"

#Adding information about feeder position 
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YYP1","position"]<-"YYP1"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YYP2","position"]<-"YYP2"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YYP3","position"]<-"YYP3"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YYP4","position"]<-"YYP4"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YYP5","position"]<-"YYP5"

visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YZP1","position"]<-"YZP1"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YZP2","position"]<-"YZP2"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YZP3","position"]<-"YZP3"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YZP4","position"]<-"YZP4"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="YZP5","position"]<-"YZP5"

visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="ZZP1","position"]<-"ZZP1"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="ZZP2","position"]<-"ZZP2"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="ZZP3","position"]<-"ZZP3"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="ZZP4","position"]<-"ZZP4"
visit_data_patch[substr(visit_data_patch$feeder_patch,1,4)=="ZZP5","position"]<-"ZZP5"

table(visit_data_patch$site)
table(visit_data_patch$position)

#Adding day of year and day of study period 
visit_data_patch$day <- yday(visit_data_patch$start) #day of year
visit_data_patch$study_day <- visit_data_patch$day - 143 #day of study period

#Remove "array" column
visit_data_patch <- subset(visit_data_patch, select = -c(array_patch))

#Adding information about time of the day 
visit_data_patch$time <- as_hms(visit_data_patch$start)
visit_data_patch$hour <- hour(visit_data_patch$start)

write.csv(visit_data_patch,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/visit_data_patch.csv", row.names = FALSE)

table(visit_data$JID, visit_data$feeder)
individuals_patches <- as.data.frame(table(visit_data$JID, visit_data$feeder))

#Average visit duration per individual data set
visit_duration_patch <- as.data.frame(aggregate(visit_data_patch$visit_duration, by = list(visit_data_patch$JID), FUN = "mean", na.rm = TRUE))
visit_duration_patch$JID <- visit_duration_patch$Group.1
visit_duration_patch$visit_duration <- visit_duration_patch$x
visit_duration_patch <- subset(visit_duration_patch, select = -c(Group.1, x))

#Total visit duration per individual data set
visit_duration_total_patch <- as.data.frame(aggregate(visit_data_patch$visit_duration, by = list(visit_data_patch$JID), FUN = "sum", na.rm = TRUE))
visit_duration_total_patch$JID <- visit_duration_total_patch$Group.1
visit_duration_total_patch$visit_duration <- visit_duration_total_patch$x
visit_duration_total_patch <- subset(visit_duration_total_patch, select = -c(Group.1, x))

#Visit number per individual data set
visit_number_patch <- as.data.frame(table(visit_data_patch$JID))
visit_number_patch$JID <- visit_number_patch$Var1
visit_number_patch$visit_number <- visit_number_patch$Freq
visit_number_patch <- subset(visit_number_patch, select = -c(Var1, Freq))

visit_number_patch$ring <- LH_ring$COMBINATION[match(visit_number_patch$JID, LH_ring$ID)]

visit_data_recent_patch <- visit_data_patch %>%
  group_by(JID) %>% 
  arrange(desc(start)) %>% 
  slice(1:1)

visit_number_patch$recent_contact <- visit_data_recent_patch$start[match(visit_number_patch$JID, visit_data_recent_patch$JID)]

visit_number_patch$recent_contact <- as.character(visit_number_patch$recent_contact)

visit_number_patch$recent_contact <- ymd_hms(visit_number_patch$recent_contact)

write.csv(visit_number,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 5/Data/PhD Ch 5 Social Bonds and Animal Culture/visits_JID.csv", row.names = FALSE)

#visit number and duration per individual data set
visits_per_indiv_patch <- merge(x = visit_number_patch, y = visit_duration_patch, by = "JID", all = TRUE)

visits_per_indiv_patch$location_ringed <- LH_ring$LOCATION[match(visits_per_indiv_patch$JID, LH_ring$ID)]

#Visits per individual per site and per feeder, preferred sites
perindivpersite_patch <- visit_data_patch  %>%  count(JID, site)
perindivpersite2_patch <- reshape(perindivpersite_patch, idvar = "JID", timevar = "site", direction = "wide")
perindivpersite3_patch <- perindivpersite2_patch
perindivpersite3_patch[is.na(perindivpersite3_patch)] <- 0.99
perindivpersite3_patch$pref <- perindivpersite3_patch$n.Y / perindivpersite3_patch$n.Z
perindivpersite3 <- perindivpersite3 %>% 
  as_tibble() %>% 
  mutate(pref.site = if_else(pref > 0.5,"Y", "Z"))

table(perindivpersite3$pref.site)

#Individuals that were detected at different sites
perindivperpos_patch <- visit_data_patch %>% count(JID, position)
indiv_1site_patch <- perindivpersite_patch %>% distinct(JID, .keep_all=TRUE)
indiv_2sites_patch <- perindivpersite_patch[perindivpersite_patch$JID %in% perindivpersite_patch$JID[duplicated(perindivpersite_patch$JID)],]
length(table(indiv_2sites_patch$JID))

patchperindiv <- as.data.frame(table(perindivperpos_patch$JID))
patchperindiv$JID <- patchperindiv$Var1
patchperindiv$patch_number <- patchperindiv$Freq
patchperindiv$Var1 <- NULL
patchperindiv$Freq <- NULL

#Individuals' preferred position
perindivperpos2_patch <- dcast(setDT(perindivperpos_patch), JID ~ position, value.var = "n")
perindivperpos2_patch$pref_pos <- colnames(perindivperpos2_patch)[apply(perindivperpos2_patch,1,which.max)]
table(perindivperpos2_patch$pref_pos)

#number of individuals in study
all_individuals_patch <- visits_per_indiv_patch

#Adding individual age
Current_year <- 2024  #Set reference year

LH$DATE = as.Date(LH$DATE, format="%d-%m-%Y")
LH$year <- as.numeric(format(LH$DATE,"%Y"))   #  #Get year from the date
Ringed <- LH %>% filter(CODE == "RINGED")   #Get only records of when birds were ringed for the first time

Ringed$known_age <- 0   #binary 0/1 do we know the exact age (e.g. birds ringed as a 6 are 0)
Ringed$min_age <- 0  #Either actual age (if known_age = 1), or minimum age (if known_age = 0) - currently as number of new years crossed.

for (i in  1:nrow(Ringed)) {
  if(Ringed[i,12] == "4"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2)
  }
  else if(Ringed[i,12] %in% c("1","1J","3","3J")){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i]+1) 
  }
  else if(Ringed[i,12] == "5"){
    Ringed$known_age[i] = 1
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +2) 
  }
  else if (Ringed[i,12] == "6"){
    Ringed$known_age[i] = 0
    Ringed$min_age[i] = (Current_year - Ringed$year[i] +3) 
  }
  else { Ringed$known_age[i] = 0
  Ringed$min_age[i] = NA }
}

sum(Ringed$known_age)
length(Ringed$known_age)

visits_per_indiv_patch$age <- Ringed$min_age[match(visits_per_indiv_patch$JID, Ringed$ID)]

all_individuals_patch$age <- Ringed$min_age[match(all_individuals_patch$JID, Ringed$ID)]
all_individuals_patch$sex <- LH_sex$SEX[match(all_individuals_patch$JID,LH_sex$ID)]
all_individuals_patch$pair_ID <- LH_pairs$pair_ID[match(all_individuals_patch$JID,LH_pairs$ID)]
all_individuals_patch$pref_pos <- perindivperpos2$pref_pos_patch[match(all_individuals_patch$JID,perindivperpos2_patch$JID)] 
all_individuals_patch$box <- LH_box_24$BOX[match(all_individuals_patch$JID, LH_box_24$ID)]
all_individuals_patch$patch_number <- patchperindiv$patch_number[match(all_individuals_patch$JID, patchperindiv$JID)]

all_individuals_patch2 <- all_individuals_patch[, c("JID", "visit_number", "visit_duration", "patch_number")]
all_individuals_patch2$visit_number_patch <- all_individuals_patch2$visit_number
all_individuals_patch2$visit_number <- NULL
all_individuals_patch2$visit_duration_patch <- all_individuals_patch2$visit_duration
all_individuals_patch2$visit_duration <- NULL

#Pre-patch data
visit_data_prepatch <- subset(visit_data, visit_data$day > 124 & visit_data$day < 141)
visit_data_prepatch <- subset(visit_data_prepatch, !visit_data_prepatch$feeder == "Y1.2")
visit_data_prepatch <- subset(visit_data_prepatch, !visit_data_prepatch$feeder == "Y3.1")
visit_data_prepatch <- subset(visit_data_prepatch, !visit_data_prepatch$feeder == "Y3.2")

visit_data_RT_prepatch <- read.csv("C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/visit_data_RT_prepatch.csv", header = T, stringsAsFactors = F)
visit_data_RT_prepatch <- subset(visit_data_RT_prepatch, !visit_data_RT_prepatch$feeder == "Y1.2")
visit_data_RT_prepatch <- subset(visit_data_RT_prepatch, !visit_data_RT_prepatch$feeder == "Y3.1")
visit_data_RT_prepatch <- subset(visit_data_RT_prepatch, !visit_data_RT_prepatch$feeder == "Y3.2")

visit_number_prepatch <- as.data.frame(table(visit_data_prepatch$JID))
visit_number_prepatch$JID <- visit_number_prepatch$Var1
visit_number_prepatch$visit_number <- visit_number_prepatch$Freq
visit_number_prepatch <- subset(visit_number_prepatch, select = -c(Var1, Freq))

#Average visit duration per individual data set
visit_duration_prepatch <- as.data.frame(aggregate(visit_data_prepatch$visit_duration, by = list(visit_data_prepatch$JID), FUN = "mean", na.rm = TRUE))
visit_duration_prepatch$JID <- visit_duration_prepatch$Group.1
visit_duration_prepatch$visit_duration <- visit_duration_prepatch$x
visit_duration_prepatch <- subset(visit_duration_prepatch, select = -c(Group.1, x))

#visit number and duration per individual data set
visits_per_indiv_prepatch <- merge(x = visit_number_prepatch, y = visit_duration_prepatch, by = "JID", all = TRUE)

#Visits per individual per site and per feeder, preferred sites
perindivpersite_prepatch <- visit_data_prepatch  %>%  count(JID, site)
perindivpersite2_prepatch <- reshape(perindivpersite_prepatch, idvar = "JID", timevar = "site", direction = "wide")
perindivpersite3_prepatch <- perindivpersite2_prepatch
perindivpersite3_prepatch[is.na(perindivpersite3_prepatch)] <- 0.99
perindivpersite3_prepatch$pref <- perindivpersite3_prepatch$n.Y / perindivpersite3_prepatch$n.Z
perindivpersite3_prepatch <- perindivpersite3_prepatch %>% 
  as_tibble() %>% 
  mutate(pref.site = if_else(pref > 0.5,"Y", "Z"))

perindivpersite2_prepatch <- perindivpersite2_prepatch %>% 
  mutate(inter_sites = rowSums(!is.na(across(n.Z:n.Y))))

table(perindivpersite3$pref.site)

#Individuals that were detected at different sites
perindivperpos_prepatch <- visit_data_prepatch %>% count(JID, position)
indiv_1site_prepatch <- perindivpersite_prepatch %>% distinct(JID, .keep_all=TRUE)
indiv_2sites_prepatch <- perindivpersite_prepatch[perindivpersite_prepatch$JID %in% perindivpersite_prepatch$JID[duplicated(perindivpersite_prepatch$JID)],]
length(table(indiv_2sites_prepatch$JID))

prepatchperindiv <- as.data.frame(table(perindivperpos_prepatch$JID))
prepatchperindiv$JID <- prepatchperindiv$Var1
prepatchperindiv$patch_number <- prepatchperindiv$Freq
prepatchperindiv$Var1 <- NULL
prepatchperindiv$Freq <- NULL

#Individuals' preferred position
perindivperpos2_prepatch <- dcast(setDT(perindivperpos_prepatch), JID ~ position, value.var = "n")
perindivperpos2_prepatch$pref_pos <- colnames(perindivperpos2_prepatch)[apply(perindivperpos2_prepatch,1,which.max)]
table(perindivperpos2_prepatch$pref_pos)
perindivperpos2_prepatch <- perindivperpos2_prepatch %>% 
  mutate(feeder_number = rowSums(!is.na(across(Y1:Z6))))

#All individuals pre-patch
all_individuals_prepatch <- visits_per_indiv_prepatch
all_individuals_prepatch$age <- Ringed$min_age[match(all_individuals_prepatch$JID, Ringed$ID)]
all_individuals_prepatch$sex <- LH_sex$SEX[match(all_individuals_prepatch$JID,LH_sex$ID)]

all_individuals_prepatch <- all_individuals_prepatch %>%
  mutate(sex =case_when(
    sex =="F" ~ 0,
    sex =="M" ~ 1,    
    is.na(sex) ~ 0.5))

all_individuals_prepatch$sex_num <- all_individuals_prepatch$sex

all_individuals_prepatch <- all_individuals_prepatch %>%
  mutate(sex_cat =case_when(
    sex == 0 ~ "female",
    sex == 1 ~ "male",    
    sex == 0.5 ~ "unsexed"))

all_individuals_prepatch$tarsus <- mean_tarsus$tarsus[match(all_individuals_prepatch$JID, mean_tarsus$JID)]

all_individuals_prepatch$body_cond <- recent_body_cond$body_cond[match(all_individuals_prepatch$JID, recent_body_cond$JID)]

all_individuals_prepatch <- merge(x = all_individuals_prepatch, y = all_individuals_patch2, by = "JID", all = TRUE)

all_individuals_prepatch$patch_binary <- ifelse(is.na(all_individuals_prepatch$patch_number), 0, 1)

all_individuals_prepatch <- all_individuals_prepatch %>%
  mutate_at(vars(patch_number, visit_number_patch, visit_duration_patch), ~replace_na(., 0))

all_individuals_prepatch$mean_OA <- visit_data_patch_OA$OA[match(all_individuals_prepatch$JID, visit_data_patch_OA$JID)]
all_individuals_prepatch$bite <- temperament$bite[match(all_individuals_prepatch$JID, temperament$JID)]
all_individuals_prepatch$feeder_number <- perindivperpos2_prepatch$feeder_number[match(all_individuals_prepatch$JID, perindivperpos2_prepatch$JID)] 
all_individuals_prepatch$inter_sites <- perindivpersite2_prepatch$inter_sites[match(all_individuals_prepatch$JID, perindivpersite2_prepatch$JID)]
all_individuals_prepatch_sub$inter_sites <- perindivpersite2_prepatch$inter_sites[match(all_individuals_prepatch_sub$JID, perindivpersite2_prepatch$JID)]
all_individuals_prepatch$visit_duration_total_patch <- visit_duration_total_patch$visit_duration[match(all_individuals_prepatch$JID, visit_duration_total_patch$JID)]
all_individuals_prepatch_sub$visit_duration_total_patch <- visit_duration_total_patch$visit_duration[match(all_individuals_prepatch_sub$JID, visit_duration_total_patch$JID)]

all_individuals_prepatch_sub$sex_num <- all_individuals_prepatch_sub$sex

all_individuals_prepatch_sub <- all_individuals_prepatch_sub %>%
  mutate(sex_cat =case_when(
    sex == 0 ~ "female",
    sex == 1 ~ "male",    
    sex == 0.5 ~ "unsexed"))

all_individuals_prepatch_sub$sex_cat <- as.factor(all_individuals_prepatch_sub$sex_cat)

sum(is.na(all_individuals_prepatch$age))
sum(is.na(all_individuals_prepatch$sex))
sum(is.na(all_individuals_prepatch$tarsus))
sum(is.na(all_individuals_prepatch$visit_number))
sum(is.na(all_individuals_prepatch$degree))
sum(is.na(all_individuals_prepatch$strength))
sum(is.na(all_individuals_prepatch$eigenvector))
sum(is.na(all_individuals_prepatch$betweenness))
sum(is.na(all_individuals_prepatch$closeness))
sum(is.na(all_individuals_prepatch$bite))

#Add centrality measures to all_individuals_prepatch first (see 3 SOCIAL NETWORK)

all_individuals_prepatch_sub <- subset(all_individuals_prepatch, !is.na(all_individuals_prepatch$degree))
all_individuals_prepatch_sub$mean_OA_noNA <- all_individuals_prepatch_sub$mean_OA
#all_individuals_prepatch_sub <- all_individuals_prepatch_sub %>% 
mutate(mean_OA_noNA = replace_na(mean_OA_noNA, 0))

#all_individuals_prepatch_sub$mean_OA_noNA <- ceiling(all_individuals_prepatch_sub$mean_OA_noNA)

write.csv(visit_data_prepatch,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/visit_data_prepatch.csv", row.names = FALSE)

all_individuals_prepatch_centrality <- all_individuals_prepatch_sub[, c("degree", "strength", "betweenness", "eigenvector")]
all_individuals_prepatch_ilv <- all_individuals_prepatch_sub[, c("age_z", "sex_cat", "tarsus_z", "body_cond_z", "visit_number_z", "pref_site")]

# (3) SOCIAL NETWORK ----

write.csv(all_individuals_prepatch_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/all_individuals_prepatch_sub.csv", row.names = FALSE)

#Gaussian Mixture Models 

#Prepare variables 
## (i) Time stamp: seconds since start of study period (1 s before first visit)
int <- interval(ymd_hms("2024-05-04 05:30:00 UTC"), ymd_hms(visit_data_prepatch$start))
visit_data_prepatch$time <- time_length(int, "second")

#int <- interval(ymd_hms("2024-05-04 05:30:00 UTC"), ymd_hms(visit_data_RT_prepatch$time))
#visit_data_RT_prepatch$time <- time_length(int, "second")

visit_data_RT_prepatch$pref_pos <- perindivperpos2_prepatch$pref_pos[match(visit_data_RT_prepatch$JID, perindivperpos2_prepatch$JID)]
visit_data_RT_prepatch$JID_pos <- paste(visit_data_RT_prepatch$JID, visit_data_RT_prepatch$position, sep = "_")
visit_data_RT_prepatch$JID_pref_pos <- paste(visit_data_RT_prepatch$JID, visit_data_RT_prepatch$pref_pos, sep = "_")
visit_data_RT_prepatch$pos_pref_pos <- ifelse(visit_data_RT_prepatch$JID_pos == visit_data_RT_prepatch$JID_pref_pos, 1, 0)
visit_data_RT_prepatch_pref <- subset(visit_data_RT_prepatch, visit_data_RT_prepatch$pos_pref_pos == 1)

## (ii) Identity 
visit_data_prepatch$JID

visit_data_RT_prepatch$JID

#Visit number per individual data set
global_ids <- unique(sort(visit_data_prepatch$JID))

global_ids <- unique(sort(visit_data_RT_prepatch$JID))

## (iii) Location 
visit_data_prepatch$position  

visit_data_RT_prepatch$position  

## (iv) Location in time
visit_data_prepatch$loc_date <-
  paste(visit_data_prepatch$position,
        visit_data_prepatch$day,sep="_")

visit_data_RT_prepatch$loc_date <-
  paste(visit_data_RT_prepatch$position,
        visit_data_RT_prepatch$day,sep="_")

visit_data_RT_prepatch_pref$loc_date <-
  paste(visit_data_RT_prepatch_pref$position,
        visit_data_RT_prepatch_pref$day,sep="_")

## Data set 
visit_data_gmm <- visit_data_prepatch[, c("time", "JID", "position", "loc_date")]
visit_data_gmm <- na.omit(visit_data_gmm)

visit_data_RT_gmm <- visit_data_RT_prepatch[, c("time", "JID", "position", "loc_date")]
visit_data_RT_gmm <- na.omit(visit_data_RT_gmm)

visit_data_RT_gmm <- visit_data_RT_prepatch_pref[, c("time", "JID", "position", "loc_date")]
visit_data_RT_gmm <- na.omit(visit_data_RT_gmm)

# Generate GMM data
gmm_data_prepatch <- gmmevents(time= visit_data_gmm$time,
                               identity=visit_data_gmm$JID,
                               location=visit_data_gmm$loc_date,
                               global_ids=global_ids)

gmm_data_RT_prepatch <- gmmevents(time= visit_data_RT_gmm$time,
                                  identity=visit_data_RT_gmm$JID,
                                  location=visit_data_RT_gmm$loc_date,
                                  global_ids=global_ids)
# Extract output
gbi <- gmm_data_prepatch$gbi
events <- gmm_data_prepatch$metadata
observations_per_event <- gmm_data_prepatch$B

gbi <- gmm_data_RT_prepatch$gbi
events <- gmm_data_RT_prepatch$metadata
observations_per_event <- gmm_data_RT_prepatch$B

# Can also subset gbi to only individuals observed
# in the dataset to give same answer as if
# global_ids had not been provided
gbi <- gbi[,which(colSums(gbi)>0)]


# Split up location and date data
tmp <- strsplit(events$Location,"_")
tmp <- do.call("rbind",tmp)
events$Location <- tmp[,1]
events$Date <- tmp[,2]

#Get the adjacency matrix from group by individual matrix 
am <- get_network(gbi, data_format = "GBI",
                  association_index = "SRI", identities = NULL,
                  which_identities = NULL, times = NULL, occurrences = NULL,
                  locations = NULL, which_locations = NULL, start_time = NULL,
                  end_time = NULL, classes = NULL, which_classes = NULL,
                  enter_time = NULL, exit_time = NULL)

am <- am[sort(rownames(am)), sort(colnames(am))]

write.csv(am,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/am.csv", row.names = FALSE)

#Get graph from the adjacency matrix
g <- graph_from_adjacency_matrix(am, mode= "undirected",weighted=TRUE,diag=FALSE)

#Get sub-graph
E(g)$weight
edge_weight <- E(g)$weight
mean(E(g)$weight) * 2
summary(edge_weight)
hist(E(g)$weight, breaks = 100)
quantile(edge_weight, probs = seq(0, 1, 1/100))
gs <- subgraph_from_edges(g, E(g)[E(g)$weight > 0.028571429], del=F)
gs <- subgraph_from_edges(g, E(g)[E(g)$weight > 0.049019608], del=F)

gs2 <- delete_vertices(gs, degree(gs)==0)

#Vertex attributes
g <- set_vertex_attr(g, "pref_pos", value = all_individuals_prepatch_sub$pref_pos)
gs2 <- set_vertex_attr(gs2, "pref_pos", value = all_individuals_prepatch_sub$pref_pos)

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

V(g)$colour <- ifelse(V(g)$pref_pos == "Y1", "black", ifelse(V(g)$pref_pos == "Y2", "red", ifelse(V(g)$pref_pos == "Y4", "yellow", ifelse(V(g)$pref_pos == "Z2", "green", ifelse(V(g)$pref_pos == "Z6", "blue", "white")))))
V(gs2)$colour <- ifelse(V(gs2)$pref_pos == "Y1", "blue", ifelse(V(gs2)$pref_pos == "Y2", "red", ifelse(V(gs2)$pref_pos == "Y3", "yellow", ifelse(V(gs2)$pref_pos == "Z2", "green", ifelse(V(gs2)$pref_pos == "Z5", "black", "white")))))

V(g)$colour <- ifelse(V(g)$patch_binary == 1, "black", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "black", "white")

V(g)$colour <- ifelse(V(g)$patch_binary == 1, "grey", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "grey", "white")
V(gs)$colour <- ifelse(V(gs)$patch_binary == 1, "darkgoldenrod1", "white")

V(g)$colour <- ifelse(V(g)$patch_number == 0, "#FFFFFF", ifelse(V(g)$patch_number == 1, "#CCCCCC", ifelse(V(g)$patch_number == 2, "#999999", ifelse(V(g)$patch_number == 3, "#666666", ifelse(V(g)$patch_number == 4, "#333333", "#000000")))))

V(g)$degree_scaled <- scales::rescale(V(g)$degree, to = c(0,1))
V(gs)$degree_scaled <- rescale(V(gs)$degree)

# Create a color palette

degree_order <- seq(min(V(gs)$degree), max(V(gs)$degree))
colours <- data.frame(colour = heat.colors(length(degree_order)), rev = T, levels = degree_order)

V(gs)$colour <- colours$colour[match(V(gs)$degree, colours$levels)]

colours <- data.frame(colorRampPalette(c("white", "black"), rev = T),levels = degree_order)

V(gs)$colour <- colorRampPalette(c("white", "black"))(V(gs)$degree)

num_colours <- length(unique(V(gs)$degree))
colour_palette <- colorRampPalette(c("white", "black"))(num_colors)

# Assign colors based on degree centrality
vertex_colors <- color_palette[rank(degree_centrality)]
coords <- layout_(g, nicely())
coords <- layout_(g, as_star())
coords <- layout_(g, in_circle())
coords <- layout_(g, as_bipartite())
coords <- layout_(g, as_tree())
coords <- layout_(g, on_grid())
coords <- layout_(g, on_sphere())
coords <- layout_(g, randomly())
coords <- layout_(g, with_dh())
coords <- layout_(g, with_fr())
coords <- layout_(g, with_gem())
coords <- layout_(g, with_graphopt())
coords <- layout_(g, with_kk())
coords <- layout_(g, with_lgl())
coords <- layout_(g, with_mds())
coords <- layout_(g, with_sugiyama())
coords <- layout_(g, merge_coords())
coords <- layout_(g, norm_coords())
coords <- layout_(g, normalize())

plot(g, layout = coords, vertex.size = 3, vertex.label = NA, vertex.color = "darkseagreen", edge.width = E(g)$weight * 15, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = "darkseagreen", edge.width = E(g)$weight * 15, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = "white", edge.width = E(g)$weight * 10, edge.curved = 0.35)
plot(g, vertex.size = 4, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

plot(g, vertex.size = 4, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

#Patch binary (colour) and scaled degree (size)
plot(g, layout = coords, vertex.size = V(g)$degree_scaled * 3 + 3, vertex.label = NA, vertex.color = V(g)$colour, edge.width = igraph::E(g)$weight * 15, edge.curved = 0.35)

#Patch number (colour) and scaled degree (size)
plot(g, layout = coords, vertex.size = V(g)$degree_scaled * 3 + 3, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

plot(gs, vertex.size = 4, vertex.label = NA, vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)
plot(gs, vertex.size = 6, vertex.label = NA, vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)
plot(g, vertex.size = 6, vertex.label = NA, vertex.color = V(g)$colour, edge.width = E(g)$weight * 15, edge.curved = 0.35)

plot(gs, vertex.size = V(gs)$degree * 0.3, vertex.label = NA, vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

plot(gs, vertex.size = V(gs)$patch_binary * 1.5 + 4, vertex.label = NA, vertex.color = colour_palette, edge.width = E(gs)$weight * 15, edge.curved = 0.35)
plot(gs, vertex.size = V(gs)$patch_binary * 1.5 + 4, vertex.label = NA, vertex.color = V(gs)$colour, edge.width = E(gs)$weight * 15, edge.curved = 0.35)

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

degree_prepatch <- as.data.frame(degree(g))
degree_prepatch$JID <- rownames(degree_prepatch)
degree_prepatch$degree <- degree_prepatch$`degree(g)`
degree_prepatch$`degree(g)` <- NULL

strength_prepatch <- as.data.frame(strength(g))
strength_prepatch$JID <- rownames(strength_prepatch)
strength_prepatch$strength <- strength_prepatch$`strength(g)`
strength_prepatch$`strength(g)` <- NULL

betweenness_prepatch <- as.data.frame(betweenness(g))
betweenness_prepatch$JID <- rownames(betweenness_prepatch)
betweenness_prepatch$betweenness <- betweenness_prepatch$`betweenness(g)`
betweenness_prepatch$`betweenness(g)` <- NULL

eigenvector_prepatch <- as.data.frame(eigen_centrality(g))
eigenvector_prepatch$JID <- rownames(eigenvector_prepatch)
eigenvector_prepatch$eigenvector <- eigenvector_prepatch$vector

closeness_prepatch <- as.data.frame(igraph::closeness(g, weights = 1 / E(g)$weight))
closeness_prepatch$JID <- rownames(closeness_prepatch)
closeness_prepatch$closeness <- closeness_prepatch$`igraph::closeness(g, weights = 1/E(g)$weight)`
closeness_prepatch$`igraph::closeness(g, weights = 1/E(g)$weight)` <- NULL

transitivity_prepatch <- as.data.frame(transitivity(g, type = "local", isolates = "zero"))
transitivity_prepatch$JID <- rownames(transitivity_prepatch)
transitivity_prepatch$transitivity <- transitivity_prepatch$`transitivity(g, type = "local", isolates = "zero")`
transitivity_prepatch$`transitivity(g, type = "local", isolates = "zero")` <- NULL

all_individuals_prepatch$degree <- degree_prepatch$degree[match(all_individuals_prepatch$JID, degree_prepatch$JID)] 
all_individuals_prepatch$strength <- strength_prepatch$strength[match(all_individuals_prepatch$JID, strength_prepatch$JID)] 
all_individuals_prepatch$betweenness <- betweenness_prepatch$betweenness[match(all_individuals_prepatch$JID, betweenness_prepatch$JID)]
all_individuals_prepatch$eigenvector <- eigenvector_prepatch$eigenvector[match(all_individuals_prepatch$JID, eigenvector_prepatch$JID)]
all_individuals_prepatch$closeness <- closeness_prepatch$closeness[match(all_individuals_prepatch$JID, closeness_prepatch$JID)]
all_individuals_prepatch$transitivity <- transitivity_prepatch$transitivity[match(all_individuals_prepatch$JID, transitivity_prepatch$JID)]

all_individuals_prepatch_sub$degree_1pos <- degree_prepatch$degree[match(all_individuals_prepatch_sub$JID, degree_prepatch$JID)] 
all_individuals_prepatch_sub$strength_1pos <- strength_prepatch$strength[match(all_individuals_prepatch_sub$JID, strength_prepatch$JID)] 
all_individuals_prepatch_sub$betweenness_1pos <- betweenness_prepatch$betweenness[match(all_individuals_prepatch_sub$JID, betweenness_prepatch$JID)]
all_individuals_prepatch_sub$eigenvector_1pos <- eigenvector_prepatch$eigenvector[match(all_individuals_prepatch_sub$JID, eigenvector_prepatch$JID)]

all_individuals_prepatch_sub$mean_strength <- all_individuals_prepatch_sub$strength/all_individuals_prepatch_sub$degree
all_individuals_prepatch_sub$mean_strength_z <- scale(all_individuals_prepatch_sub$mean_strength)

am_data <- as.data.frame(am)
am_data$mean_strength <- rowMeans(am_data)
am_data$sd_strength <- apply(am_data, 1, sd, na.rm=TRUE)
am_data$social_diff <- am_data$sd_strength / am_data$mean_strength
am_data$JID <- row.names(am_data)

all_individuals_prepatch_sub$social_diff <- am_data$social_diff[match(all_individuals_prepatch_sub$JID, am_data$JID)]
all_individuals_prepatch_sub$social_diff_1pos <- am_data$social_diff[match(all_individuals_prepatch_sub$JID, am_data$JID)]

all_individuals_prepatch_sub$pref_pos <- perindivperpos2_prepatch$pref_pos[match(all_individuals_prepatch_sub$JID, perindivperpos2_prepatch$JID)]
all_individuals_prepatch_sub$pref_site <- perindivpersite3_prepatch$pref.site[match(all_individuals_prepatch_sub$JID, perindivpersite3_prepatch$JID)]

#group size
group_size <- rowSums(gbi)
mean(group_size)
median(group_size)
sd(group_size)

individual_mean_group_size <- apply(gbi, 2, function(x) {
  # x is a vector of 0/1 for one individual across all groups
  # Only include groups where the individual was present
  if (sum(x) == 0) return(NA)  # avoid division by zero
  mean(group_size[x == 1])
})

individual_mean_group_size

all_individuals_prepatch_sub$mean_group_size <- individual_mean_group_size

average_experienced_group_size <- mean(individual_mean_group_size, na.rm = TRUE)
average_experienced_group_size

unique_groups <- unique(gbi)
unique_compositions_per_individual <- colSums(unique_groups)
unique_compositions_per_individual

all_individuals_prepatch_sub$unique_group_comp <- unique_compositions_per_individual

#Datasets for manuscript ----

write.csv(all_individuals_prepatch_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_prepatch_sub.csv", row.names = FALSE)
write.csv(visit_data_patch_all_OA,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/visit_data_patch_all_OA.csv", row.names = FALSE)
write.csv(event_data,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/event_data.csv", row.names = FALSE)
write.csv(edge_list,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/edge_list.csv", row.names = FALSE)

all_individuals_prepatch_sub <- read.csv(all_individuals_prepatch_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/all_individuals_prepatch_sub.csv", row.names = FALSE)
visit_data_patch_all_OA <- read.csv(visit_data_patch_all_OA,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/visit_data_patch_all_OA.csv", row.names = FALSE)
event_data <- read.csv(event_data,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/event_data.csv", row.names = FALSE)
edge_list <- read.csv(edge_list,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/Data for manuscript/edge_list.csv", row.names = FALSE)


#(4) FIRST EXPLORATION ----

#boxplot(all_individuals_prepatch$patch_number ~ all_individuals_prepatch$sex)
plot(all_individuals_prepatch$visit_number, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$age, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$tarsus, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$degree, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$strength, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$betweenness, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$eigenvector, all_individuals_prepatch$patch_number)
plot(all_individuals_prepatch$closeness, all_individuals_prepatch$patch_number)

boxplot(all_individuals_prepatch$visit_number_patch ~ all_individuals_prepatch$sex)
plot(all_individuals_prepatch$age, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$tarsus, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$degree, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$strength, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$betweenness, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$eigenvector, all_individuals_prepatch$visit_number_patch)
plot(all_individuals_prepatch$closeness, all_individuals_prepatch$visit_number_patch)

boxplot(all_individuals_prepatch$patch_binary ~ all_individuals_prepatch$sex)
plot(all_individuals_prepatch$age, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$tarsus, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$degree, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$strength, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$betweenness, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$eigenvector, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$closeness, all_individuals_prepatch$patch_binary)
plot(all_individuals_prepatch$strength, all_individuals_prepatch$degree)

plot(all_individuals_prepatch_sub$tarsus, all_individuals_prepatch_sub$degree)
cor.test(all_individuals_prepatch_sub$tarsus, all_individuals_prepatch_sub$degree)

table(all_individuals_prepatch_sub$age)
table(all_individuals_prepatch_sub$sex)

individuals_per_patch <- as.data.frame(table(all_individuals_all_patches$patch_ID[all_individuals_all_patches$discovered_binary == 1]))
individuals_per_patch[15,] <- NA
individuals_per_patch[15,"Freq"] <- 0 
mean(individuals_per_patch$Freq)
median(individuals_per_patch$Freq)
sd(individuals_per_patch$Freq)

#(5) STATISTICAL MODELLING ----

#Rescale numeric variables (z-transformation)
all_individuals_prepatch_sub$degree_z <- as.numeric(scale(all_individuals_prepatch_sub$degree))
all_individuals_prepatch_sub$strength_z <- as.numeric(scale(all_individuals_prepatch_sub$strength))
all_individuals_prepatch_sub$eigenvector_z <- as.numeric(scale(all_individuals_prepatch_sub$eigenvector))
all_individuals_prepatch_sub$betweenness_z <- as.numeric(scale(all_individuals_prepatch_sub$betweenness))
all_individuals_prepatch_sub$closeness_z <- as.numeric(scale(all_individuals_prepatch_sub$closeness))
all_individuals_prepatch_sub$closeness_z <- ifelse(!is.na(all_individuals_prepatch_sub$closeness_z), all_individuals_prepatch_sub$closeness_z, 0)
all_individuals_prepatch_sub$tarsus_z <- as.numeric(scale(all_individuals_prepatch_sub$tarsus))
all_individuals_prepatch_sub$body_cond_z <- as.numeric(scale(all_individuals_prepatch_sub$body_cond))
all_individuals_prepatch_sub$age_z <- as.numeric(scale(all_individuals_prepatch_sub$age))
all_individuals_prepatch_sub$visit_number_z <- as.numeric(scale(all_individuals_prepatch_sub$visit_number))
all_individuals_prepatch_sub$feeder_number_z <- as.numeric(scale(all_individuals_prepatch_sub$feeder_number))
all_individuals_prepatch_sub$feeder_number_binary <- ifelse(all_individuals_prepatch_sub$feeder_number == 1, 1, 2)
all_individuals_prepatch_sub$mean_strength_z <- as.numeric(scale(all_individuals_prepatch_sub$strength))
all_individuals_prepatch_sub$social_diff_z <- as.numeric(scale(all_individuals_prepatch_sub$social_diff))
all_individuals_prepatch_sub$social_diff_z <- ifelse(!is.na(all_individuals_prepatch_sub$social_diff_z), all_individuals_prepatch_sub$social_diff_z, 0)
all_individuals_prepatch_sub$transitivity_z <- as.numeric(scale(all_individuals_prepatch_sub$transitivity))

all_individuals_prepatch_sub$inter_sites <- as.factor(all_individuals_prepatch_sub$inter_sites)
visit_data_patch_all_OA$inter_sites <- as.factor(visit_data_patch_all_OA$inter_sites)

#Subset data based on visit number (remove those that have not been sampled often)
all_individuals_prepatch_sub_vn5 <- subset(all_individuals_prepatch_sub, all_individuals_prepatch_sub$visit_number > 4)
all_individuals_prepatch_sub_vn10 <- subset(all_individuals_prepatch_sub, all_individuals_prepatch_sub$visit_number > 9)

#PCA  for general feeder use pre experiment
#normalise data
norm_data <- scale(all_individuals_prepatch_sub[, c("visit_number", "feeder_number", "inter_sites")])

#correlation matrix
corr_matrix <- cor(norm_data)
corr_matrix

#correlation plot
ggcorrplot(corr_matrix)

#PCA diagnostics
eigen(corr_matrix)
cortest.bartlett(corr_matrix, n= 167) #correlation
cortest.mat(corr_matrix, n1= 167)
KMO(corr_matrix)

factors(corr_matrix, n.obs = 167)
fa.parallel(corr_matrix, n.obs = 167)

#PCA
data.pca <- princomp(norm_data)
data.pca <- prcomp(norm_data)

#Create dataframe with PCs
data.pca <- prcomp(norm_data, scale = FALSE)$x
data.pca <- as.data.frame(data.pca)

summary(data.pca)
print(data.pca)
data.pca$scores

#PCA loadings
data.pca$loadings[, 1:2]
data.pca$PC1
data.pca$JID <- all_individuals_prepatch_sub$JID

all_individuals_prepatch_sub$pc1 <- data.pca$PC1[match(all_individuals_prepatch_sub$JID, data.pca$JID)]

fCORT_22_sna$pc1 <- all_individuals2$pc1[match(fCORT_22_sna$JID, all_individuals2$JID)]
fCORT_22_sna$pc2 <- all_individuals2$pc2[match(fCORT_22_sna$JID, all_individuals2$JID)]

#Explore correlations among predictor variables
cor(all_individuals_prepatch_centrality)
cor(all_individuals_prepatch_ilv)
cor.test(all_individuals_prepatch_ilv$sex, all_individuals_prepatch_ilv$tarsus)
plot(all_individuals_prepatch_ilv$sex, all_individuals_prepatch_ilv$tarsus)

cor(all_individuals_prepatch_sub[, c("degree", "strength", "mean_strength", "eigenvector", "betweenness","closeness", "social_diff", "visit_number")], use = "pairwise.complete.obs")
cor(all_individuals_prepatch_sub[, c("degree", "tarsus", "age", "visit_number", "body_cond")], use = "pairwise.complete.obs")
cor(all_individuals_prepatch_sub[, c("degree", "strength", "eigenvector", "betweenness", "mean_strength", "social_diff", "visit_number", "feeder_number_binary", "inter_sites")], use = "pairwise.complete.obs")

plot(all_individuals_prepatch_sub$visit_number, all_individuals_prepatch_sub$degree)
plot(all_individuals_prepatch_sub$visit_number, all_individuals_prepatch_sub$strength)
plot(all_individuals_prepatch_sub$visit_number, all_individuals_prepatch_sub$eigenvector)
plot(all_individuals_prepatch_sub$visit_number, all_individuals_prepatch_sub$betweenness)

cor(all_individuals_prepatch_sub$feeder_number, all_individuals_prepatch_sub$visit_number)
cor(all_individuals_prepatch_sub$degree, all_individuals_prepatch_sub$visit_number)
cor(all_individuals_prepatch_sub$degree, all_individuals_prepatch_sub$feeder_number)

#Residualised degree
degree_resid <- brm(
  degree_z ~ s(visit_number_z, k = 3),
  data = all_individuals_prepatch_sub,
  iter = 6000, control = list(adapt_delta = 0.98, max_treedepth =15))

all_individuals_prepatch_sub$degree_resid <- residuals(degree_resid)[, "Estimate"]

#Residualised strength
strength_resid <- brm(
  strength_z ~ s(visit_number_z, k = 3),
  data = all_individuals_prepatch_sub,
  iter = 6000, control = list(adapt_delta = 0.98, max_treedepth =15))

all_individuals_prepatch_sub$strength_resid <- residuals(strength_resid)[, "Estimate"]

#Residualised eigenvector
eigenvector_resid <- brm(
  eigenvector_z ~ s(visit_number_z, k = 3),
  data = all_individuals_prepatch_sub,
  iter = 6000, control = list(adapt_delta = 0.98, max_treedepth =15))

all_individuals_prepatch_sub$eigenvector_resid <- residuals(eigenvector_resid)[, "Estimate"]

#Residualised betweenness
betweenness_resid <- brm(
  betweenness_z ~ s(visit_number_z, k = 3),
  data = all_individuals_prepatch_sub,
  iter = 6000, control = list(adapt_delta = 0.98, max_treedepth =15))

all_individuals_prepatch_sub$betweenness_resid <- residuals(betweenness_resid)[, "Estimate"]

write.csv(all_individuals_prepatch_sub,"C:/Users/lh868/OneDrive - University of Exeter/Corvid Connections PhD/Chapter 4/Data/PhD Ch 4 Social Bonds and Social Learning/Data/all_individuals_prepatch_sub.csv", row.names = FALSE)

#(5.1) Patch binary ---- 

#Probability to discover at least one patch (binary data)

table(all_individuals_prepatch_sub$patch_binary)
mean(all_individuals_prepatch_sub$patch_binary)
sd(all_individuals_prepatch_sub$patch_binary)

plot(all_individuals_prepatch_sub$degree, all_individuals_prepatch_sub$patch_binary)
plot(all_individuals_prepatch_sub$strength, all_individuals_prepatch_sub$patch_binary)
plot(all_individuals_prepatch_sub$eigenvector, all_individuals_prepatch_sub$patch_binary)
plot(all_individuals_prepatch_sub$betweenness, all_individuals_prepatch_sub$patch_binary)
plot(all_individuals_prepatch_sub$visit_number, all_individuals_prepatch_sub$patch_binary)

all_individuals_prepatch_sub_nool <- subset(all_individuals_prepatch_sub, all_individuals_prepatch_sub$visit_number < 600)

all_individuals_prepatch_sub1 <- subset(all_individuals_prepatch_sub, all_individuals_prepatch_sub$visit_number > 10)

#Patch binary - degree
aggregate(all_individuals_prepatch_sub$degree, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$degree, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm1_prior <- brm(patch_binary ~ degree_z + age_z + sex_cat + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm1_prior, ndraws = 1000)

#Model including degree
patch_binary_brm1 <- brm(patch_binary ~ degree_z + age_z + sex_cat + tarsus_z + 
                         body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

patch_binary_brm1 <- brm(patch_binary ~ degree_z + age_z + sex_cat + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub1, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

summary(patch_binary_brm1)

check_collinearity(patch_binary_brm1)

#Posterior predictive checks
pp_check(patch_binary_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm1)
launch_shinystan(patch_binary_brm1)

#Posterior distribution
as_draws_df(patch_binary_brm1)
mcmc_areas(patch_binary_brm1)
mcmc_intervals(patch_binary_brm1)

#Evaluation and interpretation
loo(patch_binary_brm1, moment_match = TRUE)
fitted(patch_binary_brm1, scale = "response")
conditional_effects(patch_binary_brm1)
bayes_R2(patch_binary_brm1)
hypothesis(patch_binary_brm1, "degree_z > 0")
hypothesis(patch_binary_brm1, "inter_sites2 > 0")
hypothesis(patch_binary_brm1, "visit_number_z > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm1)

#Extract predictions
fitted(patch_binary_brm1)

patch_binary_brm1 %>%
  spread_draws(b_Intercept, b_degree_z) %>%
  ggplot(aes(y = , x = b_degree_z)) +
  theme_classic(base_size = 28) +
  theme(legend.position="none", text = element_text(size = 28, family = "Garamond")) +
  labs(x = "Est (Degree)", y = "Density") +
  stat_halfeye()

cond <- conditional_effects(patch_binary_brm1, effect = "degree_z")
p <- plot(cond, points=F)
patch_binary_degree_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = degree_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Degree centrality", y = "Probability of patch discovery")
patch_binary_degree_plot

newdata1 <- all_individuals_prepatch_sub
newdata1$degree_z <- 0
pp1 <- posterior_epred(patch_binary_brm1, newdata = newdata1)
pop_draws1 <- rowMeans(pp1)
pop_draws1
quantile(pop_draws1, c(0.025, 0.5, 0.975))

newdata2 <- all_individuals_prepatch_sub
newdata2$degree_z <- 1
pp2 <- posterior_epred(patch_binary_brm1, newdata = newdata2)
pop_draws2 <- rowMeans(pp2)
pop_draws2
quantile(pop_draws2, c(0.025, 0.5, 0.975))

newdata3 <- all_individuals_prepatch_sub
newdata3$degree_z <- 2
pp3 <- posterior_epred(patch_binary_brm1, newdata = newdata3)
pop_draws3 <- rowMeans(pp3)
pop_draws3
quantile(pop_draws3, c(0.025, 0.5, 0.975))

newdata4 <- all_individuals_prepatch_sub
newdata4$degree_z <- 3
pp4 <- posterior_epred(patch_binary_brm1, newdata = newdata4)
pop_draws4 <- rowMeans(pp4)
pop_draws4
quantile(pop_draws4, c(0.025, 0.5, 0.975))

diff <- rowMeans(pp4 - pp1)
quantile(diff, c(0.025, 0.5, 0.975))

#Patch binary - strength
aggregate(all_individuals_prepatch_sub$strength, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$strength, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm2_prior <- brm(patch_binary ~ strength_z + age_z + sex_cat + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm2_prior, ndraws = 1000)

#Model including strength
patch_binary_brm2 <- brm(patch_binary ~ strength_z + age_z + sex_cat + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm2)

#Posterior predictive checks
pp_check(patch_binary_brm2, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm2)
launch_shinystan(patch_binary_brm2)

#Posterior distribution
as_draws_df(patch_binary_brm2)
mcmc_areas(patch_binary_brm2)
mcmc_intervals(patch_binary_brm2)

#Evaluation and interpretation
loo(patch_binary_brm2)
fitted(patch_binary_brm2, scale = "response")
conditional_effects(patch_binary_brm2)
bayes_R2(patch_binary_brm2)
hypothesis(patch_binary_brm2, "strength_z > 0")
hypothesis(patch_binary_brm2, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm2)

#Extract predictions
fitted(patch_binary_brm2)

cond <- conditional_effects(patch_binary_brm2, effect = "strength_z")
p <- plot(cond, points=F)
patch_binary_strength_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = strength_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Strength (weighted degree)", y = "Probability of patch discovery")
patch_binary_strength_plot

#Patch binary - eigenvector
aggregate(all_individuals_prepatch_sub$eigenvector, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$eigenvector, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm3_prior <- brm(patch_binary ~ eigenvector_z + age_z + sex_cat + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm3_prior, ndraws = 1000)

#Model including eigenvector
patch_binary_brm3 <- brm(patch_binary ~ eigenvector_z + age_z + sex_cat + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm3)

#Posterior predictive checks
pp_check(patch_binary_brm3, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm3)
launch_shinystan(patch_binary_brm3)

#Posterior distribution
as_draws_df(patch_binary_brm3)
mcmc_areas(patch_binary_brm3)
mcmc_intervals(patch_binary_brm3)

#Evaluation and interpretation
loo(patch_binary_brm3)
fitted(patch_binary_brm3, scale = "response")
conditional_effects(patch_binary_brm3)
bayes_R2(patch_binary_brm3)
hypothesis(patch_binary_brm3, "eigenvector_z > 0")
hypothesis(patch_binary_brm3, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm3)

#Extract predictions
fitted(patch_binary_brm3)

cond <- conditional_effects(patch_binary_brm3, effect = "eigenvector_z")
p <- plot(cond, points=F)
patch_binary_eigenvector_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = eigenvector_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Eigenvector centrality", y = "Probability of patch discovery")
patch_binary_eigenvector_plot

#Patch binary - betweenness 
aggregate(all_individuals_prepatch_sub$betweenness, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$betweenness, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm4_prior <- brm(patch_binary ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm4_prior, ndraws = 1000)

#Model including betweenness
patch_binary_brm4 <- brm(patch_binary ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm4)

#Posterior predictive checks
pp_check(patch_binary_brm4, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm4)
launch_shinystan(patch_binary_brm4)

#Posterior distribution
as_draws_df(patch_binary_brm4)
mcmc_areas(patch_binary_brm4)
mcmc_intervals(patch_binary_brm4)

#Evaluation and interpretation
loo(patch_binary_brm4)
fitted(patch_binary_brm4, scale = "response")
conditional_effects(patch_binary_brm4)
bayes_R2(patch_binary_brm4)
hypothesis(patch_binary_brm4, "betweenness_z > 0")
hypothesis(patch_binary_brm4, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm4)

#Extract predictions
fitted(patch_binary_brm4)

cond <- conditional_effects(patch_binary_brm4, effect = "betweenness_z")
p <- plot(cond, points=F)
patch_binary_betweenness_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = betweenness_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Betweenness centrality", y = "Probability of patch discovery")
patch_binary_betweenness_plot

#Patch binary - closeness
aggregate(all_individuals_prepatch_sub$closeness, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$closeness, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm5_prior <- brm(patch_binary ~ closeness_z + age_z + sex_cat + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm5_prior, ndraws = 1000)

#Model including closeness
patch_binary_brm5 <- brm(patch_binary ~ closeness_z + age_z + sex_cat + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.99), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm5)

#Posterior predictive checks
pp_check(patch_binary_brm5, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm5)
launch_shinystan(patch_binary_brm5)

#Posterior distribution
as_draws_df(patch_binary_brm5)
mcmc_areas(patch_binary_brm5)
mcmc_intervals(patch_binary_brm5)

#Evaluation and interpretation
loo(patch_binary_brm5)
fitted(patch_binary_brm5, scale = "response")
conditional_effects(patch_binary_brm5)
bayes_R2(patch_binary_brm5)
hypothesis(patch_binary_brm5, "closeness_z > 0")
hypothesis(patch_binary_brm5, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm5)

#Extract predictions
fitted(patch_binary_brm5)

cond <- conditional_effects(patch_binary_brm5, effect = "closeness_z")
p <- plot(cond, points=F)
patch_binary_closeness_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = closeness_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Closeness centrality", y = "Probability of patch discovery")
patch_binary_closeness_plot

#Patch binary - Social differentiation
aggregate(all_individuals_prepatch_sub$social_diff, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$social_diff, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm6_prior <- brm(patch_binary ~ social_diff_z + age_z + sex + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm6_prior, ndraws = 1000)

#Model including social differentiation
patch_binary_brm6 <- brm(patch_binary ~ social_diff_z + age_z + sex + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm6)

#Posterior predictive checks
pp_check(patch_binary_brm6, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm6)
launch_shinystan(patch_binary_brm6)

#Posterior distribution
as_draws_df(patch_binary_brm6)
mcmc_areas(patch_binary_brm6)
mcmc_intervals(patch_binary_brm6)

#Evaluation and interpretation
loo(patch_binary_brm6)
fitted(patch_binary_brm6, scale = "response")
conditional_effects(patch_binary_brm6)
bayes_R2(patch_binary_brm6)

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm6)

#Extract predictions
fitted(patch_binary_brm6)

cond <- conditional_effects(patch_binary_brm6, effect = "social_diff_z")
p <- plot(cond, points=F)
patch_binary_social_diff_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = social_diff_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Social differentiation", y = "Probability of patch discovery")
patch_binary_social_diff_plot

#Patch binary - transitivity
aggregate(all_individuals_prepatch_sub$transitivity, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "mean", na.rm = TRUE)
aggregate(all_individuals_prepatch_sub$transitivity, by = list(all_individuals_prepatch_sub$patch_binary), FUN = "sd", na.rm = TRUE)

#Model just with prior
default_prior()

patch_binary_brm7_prior <- brm(patch_binary ~ transitivity_z + age_z + sex + tarsus_z + 
                                 body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                               data = all_individuals_prepatch_sub, 
                               family = bernoulli(link = "logit"),
                               prior = c(
                                 set_prior("normal(0, 0.5)", class = "b"),        
                                 set_prior("normal(0, 1)", class = "Intercept")),
                               sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_binary_brm7_prior, ndraws = 1000)

#Model including transitivity
patch_binary_brm7 <- brm(patch_binary ~ transitivity_z + age_z + sex + tarsus_z + 
                           body_cond_z + s(visit_number_z) + inter_sites + (1|pref_pos), 
                         data = all_individuals_prepatch_sub, 
                         family = bernoulli(link = "logit"),
                         prior = c(
                           set_prior("normal(0, 0.5)", class = "b"),        
                           set_prior("normal(0, 1)", class = "Intercept")
                         ), control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE))

#Summary
summary(patch_binary_brm7)

#Posterior predictive checks
pp_check(patch_binary_brm7, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_binary_brm7)
launch_shinystan(patch_binary_brm7)

#Posterior distribution
as_draws_df(patch_binary_brm7)
mcmc_areas(patch_binary_brm7)
mcmc_intervals(patch_binary_brm7)

#Evaluation and interpretation
loo(patch_binary_brm7)
fitted(patch_binary_brm7, scale = "response")
conditional_effects(patch_binary_brm7)
bayes_R2(patch_binary_brm7)

#Sensitivity analysis and prior checks 
prior_summary(patch_binary_brm7)

#Extract predictions
fitted(patch_binary_brm7)
cond <- conditional_effects(patch_binary_brm7, effect = "transitivity_z")
p <- plot(cond, points=F)
patch_binary_transitivity_plot <- p[[1]] + 
  theme_base(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = social_diff_z, y = patch_binary), 
    data = all_individuals_prepatch_sub, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Social differentiation", y = "Probability of patch discovery")
patch_binary_social_diff_plot

summary(patch_binary_brm1)
summary(patch_binary_brm2)
summary(patch_binary_brm3)
summary(patch_binary_brm4)
summary(patch_binary_brm5)
summary(patch_binary_brm6)
summary(patch_binary_brm7)

patch_binary_loo1 <- loo(patch_binary_brm1)
patch_binary_loo2 <- loo(patch_binary_brm2)
patch_binary_loo3 <- loo(patch_binary_brm3)
patch_binary_loo4 <- loo(patch_binary_brm4)
patch_binary_loo5 <- loo(patch_binary_brm5)
patch_binary_loo6 <- loo(patch_binary_brm6)
patch_binary_loo7 <- loo(patch_binary_brm7)

loo_compare(patch_binary_loo1, patch_binary_loo2, patch_binary_loo3, patch_binary_loo4, patch_binary_loo5, patch_binary_loo6, patch_binary_loo7)
loo_compare(patch_binary_loo1, patch_binary_loo2, patch_binary_loo3, patch_binary_loo4, patch_binary_loo5)

print(patch_binary_loo1)
plot(patch_binary_loo1)
print(patch_binary_loo2)
plot(patch_binary_loo2)
print(patch_binary_loo3)
plot(patch_binary_loo3)
print(patch_binary_loo4)
plot(patch_binary_loo4)
print(patch_binary_loo5)
plot(patch_binary_loo5)
print(patch_binary_loo6)
plot(patch_binary_loo6)

pp_check(patch_binary_brm1, type="bars")
pp_check(patch_binary_brm1, type="dens")
pp_check(patch_binary_brm1, type="dens_overlay")
pp_check(patch_binary_brm1, type="hist")
pp_check(patch_binary_brm1, type="stat", stat = "mean")
pp_check(patch_binary_brm1, type="stat", stat = "sd")
pp_check(patch_binary_brm1, type="error_scatter_avg")
pp_check(patch_binary_brm1, type="intervals")
pp_check(patch_binary_brm1, type = "loo_pit_overlay")
pp_check(patch_binary_brm1, type = "loo_intervals")


#(5.2) Patch proportion ----

#Proportion of patches discovered

hist(all_individuals_prepatch_sub$patch_number)
mean(all_individuals_prepatch_sub$patch_number, na.rm = TRUE)
sd(all_individuals_prepatch_sub$patch_number, na.rm = TRUE)

all_individuals_prepatch_sub$patches <- 15
all_individuals_prepatch_sub_1patch <- subset(all_individuals_prepatch_sub, all_individuals_prepatch_sub$patch_number > 0)

#Patch proportion - degree

#Model just with prior
default_prior()
patch_prop_brm1_prior <- brm(patch_number| trials(patches) ~ degree_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm1_prior, ndraws = 1000)

#Model including degree
patch_prop_brm1 <- brm(patch_number| trials(patches) ~ degree_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.999))

patch_prop_brm1 <- brm(patch_number| trials(patches) ~ degree_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.9999))

summary(patch_prop_brm1)
check_collinearity(patch_prop_brm1)

#Posterior predictive checks
pp_check(patch_prop_brm1, type = "dens_overlay", ndraws = 100)

#Plot model
plot(patch_prop_brm1)

#Posterior distribution
as_draws_df(patch_prop_brm1)
mcmc_areas(patch_prop_brm1)
mcmc_intervals(patch_prop_brm1)

#Evaluation and interpretation
loo(patch_prop_brm1)
fitted(patch_prop_brm1, scale = "response")
conditional_effects(patch_prop_brm1)
bayes_R2(patch_prop_brm1)
hypothesis(patch_prop_brm1, "degree_z > 0")
hypothesis(patch_prop_brm1, "sex_catunsexed > 0")
hypothesis(patch_prop_brm1, "inter_sites2 > 0")
hypothesis(patch_prop_brm1, "tarsus_z > 0")
hypothesis(patch_prop_brm1, "body_cond_z > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm1)

#Extract predictions
fitted(patch_prop_brm1)

#Patch proportion - strength

#Model just with prior
default_prior()
patch_prop_brm2_prior <- brm(patch_number| trials(patches) ~ strength_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm2_prior, ndraws = 1000)

#Model including strength
patch_prop_brm2 <- brm(patch_number| trials(patches) ~ strength_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.90))

patch_prop_brm2 <- brm(patch_number| trials(patches) ~ strength_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.999))

summary(patch_prop_brm2)

check_collinearity(patch_prop_brm2)

#Posterior predictive checks
pp_check(patch_prop_brm2, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm2)

#Posterior distribution
as_draws_df(patch_prop_brm2)
mcmc_areas(patch_prop_brm2)
mcmc_intervals(patch_prop_brm2)

#Evaluation and interpretation
loo(patch_prop_brm2)
fitted(patch_prop_brm2, scale = "response")
conditional_effects(patch_prop_brm2)
bayes_R2(patch_prop_brm2)
hypothesis(patch_prop_brm2, "strength_z > 0")
hypothesis(patch_prop_brm2, "tarsus_z > 0")
hypothesis(patch_prop_brm2, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm2)

#Extract predictions
fitted(patch_prop_brm2)

#Patch proportion - eigenvector

#Model just with prior
default_prior()
patch_prop_brm3_prior <- brm(patch_number| trials(patches) ~ eigenvector_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm3_prior, ndraws = 1000)

#Model including eigenvector
patch_prop_brm3 <- brm(patch_number| trials(patches) ~ eigenvector_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

patch_prop_brm3 <- brm(patch_number| trials(patches) ~ eigenvector_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

summary(patch_prop_brm3)
check_collinearity(patch_prop_brm3)

#Posterior predictive checks
pp_check(patch_prop_brm3, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm3)

#Posterior distribution
as_draws_df(patch_prop_brm3)
mcmc_areas(patch_prop_brm3)
mcmc_intervals(patch_prop_brm3)

#Evaluation and interpretation
loo(patch_prop_brm3)
fitted(patch_prop_brm3, scale = "response")
conditional_effects(patch_prop_brm3)
bayes_R2(patch_prop_brm3)
hypothesis(patch_prop_brm3, "eigenvector_z > 0")
hypothesis(patch_prop_brm3, "tarsus_z > 0")
hypothesis(patch_prop_brm3, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm3)

#Extract predictions
fitted(patch_prop_brm3)

#Patch proportion - betweenness

#Model just with prior
default_prior()
patch_prop_brm4_prior <- brm(patch_number| trials(patches) ~ betweenness_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm4_prior, ndraws = 1000)

#Model including betweenness
patch_prop_brm4 <- brm(patch_number| trials(patches) ~ betweenness_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.90))

patch_prop_brm4 <- brm(patch_number| trials(patches) ~ betweenness_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

summary(patch_prop_brm4)
check_collinearity(patch_prop_brm4)

#Posterior predictive checks
pp_check(patch_prop_brm4, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm4)

#Posterior distribution
as_draws_df(patch_prop_brm4)
mcmc_areas(patch_prop_brm4)
mcmc_intervals(patch_prop_brm4)

#Evaluation and interpretation
loo(patch_prop_brm4)
fitted(patch_prop_brm4, scale = "response")
conditional_effects(patch_prop_brm4)
bayes_R2(patch_prop_brm4)
hypothesis(patch_prop_brm4, "betweenness_z > 0")
hypothesis(patch_prop_brm4, "tarsus_z > 0")
hypothesis(patch_prop_brm4, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm4)

#Extract predictions
fitted(patch_prop_brm4)


#Patch proportion - closeness

#Model just with prior
default_prior()
patch_prop_brm5_prior <- brm(patch_number| trials(patches) ~ closeness_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm5_prior, ndraws = 1000)

#Model including betweenness
patch_prop_brm5 <- brm(patch_number| trials(patches) ~ closeness_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

patch_prop_brm5 <- brm(patch_number| trials(patches) ~ closeness_z + age_z + 
                         sex_cat + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

summary(patch_prop_brm5)
check_collinearity(patch_prop_brm5)

#Posterior predictive checks
pp_check(patch_prop_brm5, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm5)

#Posterior distribution
as_draws_df(patch_prop_brm5)
mcmc_areas(patch_prop_brm5)
mcmc_intervals(patch_prop_brm5)

#Evaluation and interpretation
loo(patch_prop_brm5)
fitted(patch_prop_brm5, scale = "response")
conditional_effects(patch_prop_brm5)
bayes_R2(patch_prop_brm5)
hypothesis(patch_prop_brm5, "closeness_z > 0")
hypothesis(patch_prop_brm5, "tarsus_z > 0")
hypothesis(patch_prop_brm5, "inter_sites2 > 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm5)

#Extract predictions
fitted(patch_prop_brm5)

#Patch proportion - social differentiation

#Model just with prior
default_prior()
patch_prop_brm6_prior <- brm(patch_number| trials(patches) ~ social_diff_z + age_z + 
                               sex_cat + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm4_prior, ndraws = 1000)

#Model including social differentiation
patch_prop_brm6 <- brm(patch_number| trials(patches) ~ social_diff_z + age_z + 
                         sex + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.90))

patch_prop_brm6 <- brm(patch_number| trials(patches) ~ social_diff_z + age_z + 
                         sex + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

summary(patch_prop_brm6)
check_collinearity(patch_prop_brm6)

#Posterior predictive checks
pp_check(patch_prop_brm6, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm6)

#Posterior distribution
as_draws_df(patch_prop_brm6)
mcmc_areas(patch_prop_brm6)
mcmc_intervals(patch_prop_brm6)

#Evaluation and interpretation
loo(patch_prop_brm6)
fitted(patch_prop_brm6, scale = "response")
conditional_effects(patch_prop_brm6)
bayes_R2(patch_prop_brm6)

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm6)

#Extract predictions
fitted(patch_prop_brm6)

#Patch proportion - transitivity

#Model just with prior
default_prior()
patch_prop_brm7_prior <- brm(patch_number| trials(patches) ~ transitivity_z + age_z + 
                               sex + tarsus_z + body_cond_z + s(visit_number_z) +
                               inter_sites + (1|pref_pos), 
                             data = all_individuals_prepatch_sub, 
                             family = binomial(link = "logit"),
                             prior = c(
                               set_prior("normal(0, 2)", class = "b"),        
                               set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                             sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_prop_brm7_prior, ndraws = 1000)

#Model including social differentiation
patch_prop_brm7 <- brm(patch_number| trials(patches) ~ transitivity_z + age_z + 
                         sex + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.99))

patch_prop_brm7 <- brm(patch_number| trials(patches) ~ transitivity_z + age_z + 
                         sex + tarsus_z + body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|pref_pos), 
                       data = all_individuals_prepatch_sub_1patch, 
                       family = binomial(link = "logit"),
                       prior = c(
                         set_prior("normal(0, 2)", class = "b"),        
                         set_prior("normal(-2.5, 1.0)", class = "Intercept")),
                       control = list(adapt_delta = 0.999))

summary(patch_prop_brm7)
check_collinearity(patch_prop_brm7)

#Posterior predictive checks
pp_check(patch_prop_brm7, type = "dens_overlay", ndraws = 1000)

#Plot model
plot(patch_prop_brm7)

#Posterior distribution
as_draws_df(patch_prop_brm7)
mcmc_areas(patch_prop_brm7)
mcmc_intervals(patch_prop_brm7)

#Evaluation and interpretation
loo(patch_prop_brm7)
fitted(patch_prop_brm7, scale = "response")
conditional_effects(patch_prop_brm7)
bayes_R2(patch_prop_brm7)

#Sensitivity analysis and prior checks 
prior_summary(patch_prop_brm7)

#Extract predictions
fitted(patch_prop_brm7)

#Compare models 
summary(patch_prop_brm1)
summary(patch_prop_brm2)
summary(patch_prop_brm3)
summary(patch_prop_brm4)
summary(patch_prop_brm5)
summary(patch_prop_brm6)
summary(patch_prop_brm7)

patch_prop_loo1 <- loo(patch_prop_brm1)
patch_prop_loo2 <- loo(patch_prop_brm2)
patch_prop_loo3 <- loo(patch_prop_brm3)
patch_prop_loo4 <- loo(patch_prop_brm4)
patch_prop_loo5 <- loo(patch_prop_brm5)
patch_prop_loo6 <- loo(patch_prop_brm6)
patch_prop_loo7 <- loo(patch_prop_brm7)

loo_compare(patch_prop_loo1, patch_prop_loo2, patch_prop_loo3, patch_prop_loo4, patch_prop_loo5, patch_prop_loo6, patch_prop_loo7)

print(patch_prop_loo1)
plot(patch_prop_loo1)
print(patch_prop_loo2)
plot(patch_prop_loo2)
print(patch_prop_loo3)
plot(patch_prop_loo3)
print(patch_prop_loo4)
plot(patch_prop_loo4)
print(patch_prop_loo5)
plot(patch_prop_loo5)
print(patch_prop_loo6)
plot(patch_prop_loo6)
print(patch_prop_loo7)
plot(patch_prop_loo7)


#(5.3) Patch order/time ----

visit_data_patch_all_OA_sub1 <- subset(visit_data_patch_all_OA, visit_data_patch_all_OA$visit_number >10)

#Relative order of patch discovery

#Model just with prior
default_prior()

visit_data_patch_all_OA$OA_log <- log(visit_data_patch_all_OA$OA)

patch_order_brm1_prior <- brm(OA ~ degree_z + age_z + sex_cat + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_time_brm1_prior <- brm(time_tada_hr ~ degree_z + age_z + sex_cat + tarsus_z + 
                         body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|JID) + (1|position), 
                       data = visit_data_patch_all_OA, 
                       family = lognormal(),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "b"),        
                         set_prior("normal(3.5, 1)", class = "Intercept"),
                         set_prior("student_t(3, 0, 1)", class = "sd"),
                         set_prior("student_t(3, 0, 1)", class = "sigma")
                       ), control = list(adapt_delta = 0.99),
                       sample_prior = "only"
)

#Prior predictive checks 
pp_check(patch_order_brm1_prior, ndraws = 1000)

yrep <- posterior_predict(patch_order_brm1_prior)
summary(yrep)
quantile(yrep, probs = c(0.025, 0.5, 0.975))

#Model including degree

patch_order_brm1 <- brm(OA ~ degree_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

patch_time_brm1 <- brm(time_tada_hr ~ degree_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(3.5, 1)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 1)", class = "sigma")
                        ), control = list(adapt_delta = 0.99))

summary(patch_time_brm1)

check_collinearity(patch_time_brm1)

#Posterior predictive checks
pp_check(patch_time_brm1, ndraws = 1000)

pp_check(patch_time_brm1, type = "dens_overlay", ndraws = 1000) + scale_x_log10()
pp_check(patch_time_brm1, type="stat")

#Plot model
plot(patch_time_brm1)

#Posterior distribution
as_draws_df(patch_time_brm1)
mcmc_areas(patch_time_brm1)
mcmc_intervals(patch_time_brm1)

#Evaluation and interpretation
loo(patch_time_brm1)
fitted(patch_time_brm1, scale = "response")
conditional_effects(patch_time_brm1)
bayes_R2(patch_time_brm1)
hypothesis(patch_time_brm1, "degree_z < 0")
hypothesis(patch_time_brm1, "tarsus_z < 0")
hypothesis(patch_time_brm1, "body_cond_z < 0")
hypothesis(patch_time_brm1, "inter_sites < 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_time_brm1)

#Extract predictions
fitted(patch_time_brm1)

#Model including strength

patch_order_brm2_prior <- brm(OA ~ strength_z + age_z + sex_cat + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_time_brm2_prior <- brm(time_tada_hr ~ strength_z + age_z + sex_cat + tarsus_z + 
                               body_cond_z + s(visit_number_z) + 
                               inter_sites + (1|JID) + (1|position), 
                             data = visit_data_patch_all_OA, 
                             family = lognormal(),
                             prior = c(
                               set_prior("normal(0, 0.5)", class = "b"),        
                               set_prior("normal(3.5, 1)", class = "Intercept"),
                               set_prior("student_t(3, 0, 1)", class = "sd"),
                               set_prior("student_t(3, 0, 1)", class = "sigma")
                             ), control = list(adapt_delta = 0.99),
                             sample_prior = "only"
)

patch_order_brm2 <- brm(OA ~ strength_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

patch_time_brm2 <- brm(time_tada_hr ~ strength_z + age_z + sex_cat + tarsus_z + 
                         body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|JID) + (1|position), 
                       data = visit_data_patch_all_OA, 
                       family = lognormal(),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "b"),        
                         set_prior("normal(3.5, 1)", class = "Intercept"),
                         set_prior("student_t(3, 0, 1)", class = "sd"),
                         set_prior("student_t(3, 0, 1)", class = "sigma")
                       ), control = list(adapt_delta = 0.99))

summary(patch_time_brm2)

check_collinearity(patch_time_brm2)

#Posterior predictive checks
pp_check(patch_time_brm2, ndraws = 1000)
pp_check(patch_time_brm2, type = "dens_overlay", ndraws = 1000)
pp_check(patch_time_brm2, type="stat")

#Plot model
plot(patch_time_brm2)

#Posterior distribution
as_draws_df(patch_time_brm2)
mcmc_areas(patch_time_brm2)
mcmc_intervals(patch_time_brm2)

#Evaluation and interpretation
loo(patch_time_brm2)
fitted(patch_time_brm2, scale = "response")
conditional_effects(patch_time_brm2)
bayes_R2(patch_time_brm2)
hypothesis(patch_time_brm2, "strength_z < 0")
hypothesis(patch_time_brm2, "tarsus_z < 0")
hypothesis(patch_time_brm2, "body_cond_z < 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_time_brm2)

#Extract predictions
fitted(patch_time_brm2)

#Model including eigenvector

patch_order_brm3_prior <- brm(OA ~ strength_z + age_z + sex + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_time_brm3_prior <- brm(time_tada_hr ~ eigenvector_z + age_z + sex_cat + tarsus_z + 
                               body_cond_z + s(visit_number_z) + 
                               inter_sites + (1|JID) + (1|position), 
                             data = visit_data_patch_all_OA, 
                             family = lognormal(),
                             prior = c(
                               set_prior("normal(0, 0.5)", class = "b"),        
                               set_prior("normal(3.5, 1)", class = "Intercept"),
                               set_prior("student_t(3, 0, 1)", class = "sd"),
                               set_prior("student_t(3, 0, 1)", class = "sigma")
                             ), control = list(adapt_delta = 0.99),
                             sample_prior = "only"
)

patch_order_brm3 <- brm(OA ~ eigenvector_z + age_z + sex + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

patch_time_brm3 <- brm(time_tada_hr ~ eigenvector_z + age_z + sex_cat + tarsus_z + 
                         body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|JID) + (1|position), 
                       data = visit_data_patch_all_OA, 
                       family = lognormal(),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "b"),        
                         set_prior("normal(3.5, 1)", class = "Intercept"),
                         set_prior("student_t(3, 0, 1)", class = "sd"),
                         set_prior("student_t(3, 0, 1)", class = "sigma")
                       ), control = list(adapt_delta = 0.99))

summary(patch_time_brm3)

check_collinearity(patch_time_brm3)

#Posterior predictive checks
pp_check(patch_time_brm3, ndraws = 1000)
pp_check(patch_time_brm3, type = "dens_overlay", ndraws = 1000)
pp_check(patch_time_brm3, type="stat")

#Plot model
plot(patch_time_brm3)

#Posterior distribution
as_draws_df(patch_time_brm3)
mcmc_areas(patch_time_brm3)
mcmc_intervals(patch_time_brm3)

#Evaluation and interpretation
loo(patch_time_brm3)
fitted(patch_time_brm3, scale = "response")
conditional_effects(patch_time_brm3)
bayes_R2(patch_time_brm3)
hypothesis(patch_time_brm3, "eigenvector_z < 0")
hypothesis(patch_time_brm3, "tarsus_z < 0")
hypothesis(patch_time_brm3, "body_cond_z < 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_time_brm3)

#Extract predictions
fitted(patch_time_brm3)

#Model including betweenness

patch_order_brm4_prior <- brm(OA ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_time_brm4_prior <- brm(time_tada_hr ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                               body_cond_z + s(visit_number_z) + 
                               inter_sites + (1|JID) + (1|position), 
                             data = visit_data_patch_all_OA, 
                             family = lognormal(),
                             prior = c(
                               set_prior("normal(0, 0.5)", class = "b"),        
                               set_prior("normal(3.5, 1)", class = "Intercept"),
                               set_prior("student_t(3, 0, 1)", class = "sd"),
                               set_prior("student_t(3, 0, 1)", class = "sigma")
                             ), control = list(adapt_delta = 0.99),
                             sample_prior = "only"
)

patch_order_brm4 <- brm(OA ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

patch_time_brm4 <- brm(time_tada_hr ~ betweenness_z + age_z + sex_cat + tarsus_z + 
                         body_cond_z + s(visit_number_z) + 
                         inter_sites + (1|JID) + (1|position), 
                       data = visit_data_patch_all_OA, 
                       family = lognormal(),
                       prior = c(
                         set_prior("normal(0, 0.5)", class = "b"),        
                         set_prior("normal(3.5, 1)", class = "Intercept"),
                         set_prior("student_t(3, 0, 1)", class = "sd"),
                         set_prior("student_t(3, 0, 1)", class = "sigma")
                       ), control = list(adapt_delta = 0.99))

summary(patch_time_brm4)

check_collinearity(patch_time_brm4)

#Posterior predictive checks
pp_check(patch_time_brm4, ndraws = 1000)
pp_check(patch_time_brm4, type = "dens_overlay", ndraws = 1000)
pp_check(patch_time_brm4, type="stat")

#Plot model
plot(patch_time_brm4)

#Posterior distribution
as_draws_df(patch_time_brm4)
mcmc_areas(patch_time_brm4)
mcmc_intervals(patch_time_brm4)

#Evaluation and interpretation
loo(patch_time_brm4)
fitted(patch_time_brm4, scale = "response")
conditional_effects(patch_time_brm4)
bayes_R2(patch_time_brm4)
hypothesis(patch_time_brm4, "betweenness_z > 0")
hypothesis(patch_time_brm4, "tarsus_z < 0")
hypothesis(patch_time_brm4, "body_cond_z < 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_time_brm4)

#Extract predictions
fitted(patch_time_brm4)

#Model including closeness

patch_order_brm5_prior <- brm(OA ~ closeness_z + age_z + sex_cat + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_time_brm5_prior <- brm(time_tada_hr ~ closeness_z + age_z + sex_cat + tarsus_z + 
                               body_cond_z + s(visit_number_z) + 
                               inter_sites + (1|JID) + (1|position), 
                             data = visit_data_patch_all_OA, 
                             family = lognormal(),
                             prior = c(
                               set_prior("normal(0, 0.5)", class = "b"),        
                               set_prior("normal(3.5, 1)", class = "Intercept"),
                               set_prior("student_t(3, 0, 1)", class = "sd"),
                               set_prior("student_t(3, 0, 1)", class = "sigma")
                             ), control = list(adapt_delta = 0.99),
                             sample_prior = "only"
)

patch_order_brm5 <- brm(OA ~ closeness_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999990))

patch_time_brm5 <- brm(time_tada_hr ~ closeness_z + age_z + sex_cat + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(3.5, 1)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 1)", class = "sigma")
                        ), control = list(adapt_delta = 0.99))

summary(patch_time_brm5)

check_collinearity(patch_time_brm5)

#Posterior predictive checks
pp_check(patch_time_brm5, ndraws = 100)
pp_check(patch_time_brm5, type = "dens_overlay", ndraws = 1000)
pp_check(patch_time_brm5, type="stat")

#Plot model
plot(patch_time_brm5)

#Posterior distribution
as_draws_df(patch_time_brm5)
mcmc_areas(patch_time_brm5)
mcmc_intervals(patch_time_brm5)

#Evaluation and interpretation
loo(patch_time_brm5)
fitted(patch_time_brm5, scale = "response")
conditional_effects(patch_time_brm5)
bayes_R2(patch_time_brm5)
hypothesis(patch_time_brm5, "closeness_z < 0")
hypothesis(patch_time_brm5, "tarsus_z < 0")
hypothesis(patch_time_brm5, "body_cond_z < 0")
hypothesis(patch_time_brm5, "inter_sites2 < 0")

#Sensitivity analysis and prior checks 
prior_summary(patch_time_brm5)

#Extract predictions
fitted(patch_time_brm5)

patch_time_brm5 %>%
  spread_draws(b_Intercept, b_closeness_z) %>%
  ggplot(aes(y = , x = b_closeness_z)) +
  theme_classic(base_size = 28) +
  theme(legend.position="none", text = element_text(size = 28, family = "Garamond")) +
  labs(x = "Est (Closeness)", y = "Density") +
  stat_halfeye()

cond <- conditional_effects(patch_time_brm5, effect = "closeness_z")
p <- plot(cond, points=F)
patch_time_closeness_plot <- p[[1]] + 
  theme_few(base_size = 14) +
  theme(legend.position="none", text = element_text(size = 14, family = "Garamond")) +
  geom_point(
    aes(x = closeness_z, y = time_tada_hr), 
    data = visit_data_patch_all_OA, 
    color = "black",
    size = 2,
    alpha = 0.2,
    inherit.aes = FALSE) + 
  geom_line(color="black", size=1) + 
  labs(x = "Closeness centrality", y = "Time of patch discovery (hrs)")
patch_time_closeness_plot


#Model including social differentiation

patch_order_brm6_prior <- brm(OA ~ social_diff_z + age_z + sex + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_order_brm6 <- brm(OA ~ social_diff_z + age_z + sex + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

summary(patch_order_brm6)

check_collinearity(patch_order_brm6)

#Posterior predictive checks
pp_check(patch_order_brm6, ndraws = 1000)
pp_check(patch_order_brm6, type = "dens_overlay", ndraws = 1000)
pp_check(patch_order_brm6, type="stat")

#Plot model
plot(patch_order_brm6)

#Posterior distribution
as_draws_df(patch_order_brm6)
mcmc_areas(patch_order_brm6)
mcmc_intervals(patch_order_brm6)

#Evaluation and interpretation
loo(patch_order_brm6)
fitted(patch_order_brm6, scale = "response")
conditional_effects(patch_order_brm6)
bayes_R2(patch_order_brm6)

#Sensitivity analysis and prior checks 
prior_summary(patch_order_brm6)

#Extract predictions
fitted(patch_order_brm6)

#Model including transitivity

patch_order_brm7_prior <- brm(OA ~ transitivity_z + age_z + sex + tarsus_z + 
                                body_cond_z + s(visit_number_z) +
                                inter_sites + (1|JID) + (1|position), 
                              data = visit_data_patch_all_OA, 
                              family = lognormal(),
                              prior = c(
                                set_prior("normal(0, 0.5)", class = "b"),        
                                set_prior("normal(2, 0.5)", class = "Intercept"),
                                set_prior("student_t(3, 0, 1)", class = "sd"),
                                set_prior("student_t(3, 0, 0.5)", class = "sigma")),
                              sample_prior = "only"
)

patch_order_brm7 <- brm(OA ~ transitivity_z + age_z + sex + tarsus_z + 
                          body_cond_z + s(visit_number_z) + 
                          inter_sites + (1|JID) + (1|position), 
                        data = visit_data_patch_all_OA, 
                        family = lognormal(),
                        prior = c(
                          set_prior("normal(0, 0.5)", class = "b"),        
                          set_prior("normal(2, 0.5)", class = "Intercept"),
                          set_prior("student_t(3, 0, 1)", class = "sd"),
                          set_prior("student_t(3, 0, 0.5)", class = "sigma")
                        ), control = list(adapt_delta = 0.999))

summary(patch_order_brm7)

check_collinearity(patch_order_brm7)

#Posterior predictive checks
pp_check(patch_order_brm7, ndraws = 1000)
pp_check(patch_order_brm7, type = "dens_overlay", ndraws = 1000)
pp_check(patch_order_brm7, type="stat")

#Plot model
plot(patch_order_brm7)

#Posterior distribution
as_draws_df(patch_order_brm7)
mcmc_areas(patch_order_brm7)
mcmc_intervals(patch_order_brm7)

#Evaluation and interpretation
loo(patch_order_brm7)
fitted(patch_order_brm7, scale = "response")
conditional_effects(patch_order_brm7)
bayes_R2(patch_order_brm7)

#Sensitivity analysis and prior checks 
prior_summary(patch_order_brm7)

#Extract predictions
fitted(patch_order_brm7)

#Compare models
summary(patch_order_brm1)
summary(patch_order_brm2)
summary(patch_order_brm3)
summary(patch_order_brm4)
summary(patch_order_brm5)
summary(patch_order_brm6)
summary(patch_order_brm7)

patch_order_loo1 <- loo(patch_order_brm1)
patch_order_loo2 <- loo(patch_order_brm2)
patch_order_loo3 <- loo(patch_order_brm3)
patch_order_loo4 <- loo(patch_order_brm4)
patch_order_loo5 <- loo(patch_order_brm5)
patch_order_loo6 <- loo(patch_order_brm6)
patch_order_loo7 <- loo(patch_order_brm7)

loo_compare(patch_order_loo1, patch_order_loo2, patch_order_loo3, patch_order_loo4, patch_order_loo5, patch_order_loo6, patch_order_loo7)

print(patch_order_loo1)
plot(patch_order_loo1)
print(patch_order_loo2)
plot(patch_order_loo2)
print(patch_order_loo3)
plot(patch_order_loo3)
print(patch_order_loo4)
plot(patch_order_loo4)
print(patch_order_loo5)
plot(patch_order_loo5)
print(patch_order_loo6)
plot(patch_order_loo6)
print(patch_order_loo6)
plot(patch_order_loo6)
print(patch_order_loo7)
plot(patch_order_loo7)


#(6) NETWORK-BASED DIFFUSION ANALYSIS ----

#Patch visit data ----

all_individuals_prepatch_JID <- visit_number_prepatch
all_individuals_prepatch_JID <- as.data.frame(all_individuals_prepatch_JID$JID)
all_individuals_prepatch_JID$JID <- all_individuals_prepatch_JID$`all_individuals_prepatch_JID$JID`
all_individuals_prepatch_JID$`all_individuals_prepatch_JID$JID` <- NULL

all_individuals_prepatch_ilv <- all_individuals_prepatch_sub[, c("age", "sex_cat", "tarsus", "visit_number")]

all_individuals_prepatch_diffusions <- all_individuals_prepatch_JID

#Order of acquisition of YYP1
visit_data_patch_YYP1_OA <- subset(visit_data_patch, visit_data_patch$position == "YYP1")
visit_data_patch_YYP1_OA <- visit_data_patch_YYP1_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YYP1_OA <- merge(visit_data_patch_YYP1_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YYP1_OA <- subset(visit_data_patch_YYP1_OA, !is.na(visit_data_patch_YYP1_OA$start))
visit_data_patch_YYP1_OA <- visit_data_patch_YYP1_OA %>% arrange(start) 
visit_data_patch_YYP1_OA$OA <- ifelse(is.na(visit_data_patch_YYP1_OA$start), NA, row.names(visit_data_patch_YYP1_OA))  
#visit_data_patch_YYP1_OA$OA <- row.names(visit_data_patch_YYP1_OA)  
#visit_data_patch_YYP1_OA <- visit_data_patch_YYP1_OA %>% arrange(JID) 
visit_data_patch_YYP1_OA$OA <- as.numeric(visit_data_patch_YYP1_OA$OA)
visit_data_patch_YYP1_OA$am_number <- am_ind$number[match(visit_data_patch_YYP1_OA$JID, am_ind$JID)]
visit_data_patch_YYP1_OA$patch_deployed <- ymd_hms("2024-05-22 04:59:00")
visit_data_patch_YYP1_OA$patch_removed <- ymd_hms("2024-05-27 17:14:00")
visit_data_patch_YYP1_OA$patch_last_read <- ymd_hms("2024-05-25 19:58:00")
visit_data_patch_YYP1_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YYP1_OA$patch_deployed), ymd_hms(visit_data_patch_YYP1_OA$patch_last_read))
visit_data_patch_YYP1_OA$patch_active <- as.duration(visit_data_patch_YYP1_OA$patch_interval)
visit_data_patch_YYP1_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YYP1_OA$patch_deployed), ymd_hms(visit_data_patch_YYP1_OA$start))
visit_data_patch_YYP1_OA$time_tada <- time_length(visit_data_patch_YYP1_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada1 <- visit_data_patch_YYP1_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YYP1_OA$JID)]
all_individuals_prepatch_diffusions$tada1 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada1), all_individuals_prepatch_diffusions$tada1, visit_data_patch_YYP1_OA$patch_active + 3600)

#Order of acquisition of YYP2
visit_data_patch_YYP2_OA <- subset(visit_data_patch, visit_data_patch$position == "YYP2")
visit_data_patch_YYP2_OA <- visit_data_patch_YYP2_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YYP2_OA <- merge(visit_data_patch_YYP2_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YYP2_OA <- subset(visit_data_patch_YYP2_OA, !is.na(visit_data_patch_YYP2_OA$start))
visit_data_patch_YYP2_OA <- visit_data_patch_YYP2_OA %>% arrange(start) 
visit_data_patch_YYP2_OA$OA <- ifelse(is.na(visit_data_patch_YYP2_OA$start), NA, row.names(visit_data_patch_YYP2_OA))  
#visit_data_patch_YYP2_OA$OA <- row.names(visit_data_patch_YYP2_OA)  
#visit_data_patch_YYP2_OA <- visit_data_patch_YYP2_OA %>% arrange(JID) 
visit_data_patch_YYP2_OA$OA <- as.numeric(visit_data_patch_YYP2_OA$OA)
visit_data_patch_YYP2_OA$am_number <- am_ind$number[match(visit_data_patch_YYP2_OA$JID, am_ind$JID)]
visit_data_patch_YYP2_OA$patch_deployed <- ymd_hms("2024-05-28 06:41:00")
visit_data_patch_YYP2_OA$patch_removed <- ymd_hms("2024-06-02 06:35:00")
visit_data_patch_YYP2_OA$patch_last_read <- ymd_hms("2024-06-02 06:35:00")
visit_data_patch_YYP2_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YYP2_OA$patch_deployed), ymd_hms(visit_data_patch_YYP2_OA$patch_last_read))
visit_data_patch_YYP2_OA$patch_active <- as.duration(visit_data_patch_YYP2_OA$patch_interval)
visit_data_patch_YYP2_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YYP2_OA$patch_deployed), ymd_hms(visit_data_patch_YYP2_OA$start))
visit_data_patch_YYP2_OA$time_tada <- time_length(visit_data_patch_YYP2_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada2 <- visit_data_patch_YYP2_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YYP2_OA$JID)]
all_individuals_prepatch_diffusions$tada2 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada2), all_individuals_prepatch_diffusions$tada2, visit_data_patch_YYP2_OA$patch_active + 3600)

#Order of acquisition of YYP3
visit_data_patch_YYP3_OA <- subset(visit_data_patch, visit_data_patch$position == "YYP3")
visit_data_patch_YYP3_OA <- visit_data_patch_YYP3_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YYP3_OA <- merge(visit_data_patch_YYP3_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YYP3_OA <- subset(visit_data_patch_YYP3_OA, !is.na(visit_data_patch_YYP3_OA$start))
visit_data_patch_YYP3_OA <- visit_data_patch_YYP3_OA %>% arrange(start) 
visit_data_patch_YYP3_OA$OA <- ifelse(is.na(visit_data_patch_YYP3_OA$start), NA, row.names(visit_data_patch_YYP3_OA))  
#visit_data_patch_YYP3_OA$OA <- row.names(visit_data_patch_YYP2_OA)  
#visit_data_patch_YYP3_OA <- visit_data_patch_YYP3_OA %>% arrange(JID) 
visit_data_patch_YYP3_OA$OA <- as.numeric(visit_data_patch_YYP3_OA$OA)
visit_data_patch_YYP3_OA$am_number <- am_ind$number[match(visit_data_patch_YYP3_OA$JID, am_ind$JID)]
visit_data_patch_YYP3_OA$patch_deployed <- ymd_hms("2024-06-02 07:08:00")
visit_data_patch_YYP3_OA$patch_removed <- ymd_hms("2024-06-06 18:02:00")
visit_data_patch_YYP3_OA$patch_last_read <- ymd_hms("2024-06-06 18:02:00") #feeder functional upon arrival
visit_data_patch_YYP3_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YYP3_OA$patch_deployed), ymd_hms(visit_data_patch_YYP3_OA$patch_last_read))
visit_data_patch_YYP3_OA$patch_active <- as.duration(visit_data_patch_YYP3_OA$patch_interval)
visit_data_patch_YYP3_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YYP3_OA$patch_deployed), ymd_hms(visit_data_patch_YYP3_OA$start))
visit_data_patch_YYP3_OA$time_tada <- time_length(visit_data_patch_YYP3_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada3 <- visit_data_patch_YYP3_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YYP3_OA$JID)]
all_individuals_prepatch_diffusions$tada3 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada3), all_individuals_prepatch_diffusions$tada3, visit_data_patch_YYP3_OA$patch_active + 3600)

#Order of acquisition of YYP4
visit_data_patch_YYP4_OA <- subset(visit_data_patch, visit_data_patch$position == "YYP4")
visit_data_patch_YYP4_OA <- visit_data_patch_YYP4_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YYP4_OA <- merge(visit_data_patch_YYP4_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YYP4_OA <- subset(visit_data_patch_YYP4_OA, !is.na(visit_data_patch_YYP4_OA$start))
visit_data_patch_YYP4_OA <- visit_data_patch_YYP4_OA %>% arrange(start) 
visit_data_patch_YYP4_OA$OA <- ifelse(is.na(visit_data_patch_YYP4_OA$start), NA, row.names(visit_data_patch_YYP4_OA))  
#visit_data_patch_YYP4_OA$OA <- row.names(visit_data_patch_YYP4_OA)  
#visit_data_patch_YYP4_OA <- visit_data_patch_YYP4_OA %>% arrange(JID) 
visit_data_patch_YYP4_OA$OA <- as.numeric(visit_data_patch_YYP4_OA$OA)
visit_data_patch_YYP4_OA$am_number <- am_ind$number[match(visit_data_patch_YYP4_OA$JID, am_ind$JID)]
visit_data_patch_YYP4_OA$patch_deployed <- ymd_hms("2024-06-07 05:45:00")
visit_data_patch_YYP4_OA$patch_removed <- ymd_hms("2024-06-11 21:54:00")
visit_data_patch_YYP4_OA$patch_last_read <- ymd_hms("2024-06-11 21:54:00") #feeder functional upon arrival
visit_data_patch_YYP4_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YYP4_OA$patch_deployed), ymd_hms(visit_data_patch_YYP4_OA$patch_last_read))
visit_data_patch_YYP4_OA$patch_active <- as.duration(visit_data_patch_YYP4_OA$patch_interval)
visit_data_patch_YYP4_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YYP4_OA$patch_deployed), ymd_hms(visit_data_patch_YYP4_OA$start))
visit_data_patch_YYP4_OA$time_tada <- time_length(visit_data_patch_YYP4_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada4 <- visit_data_patch_YYP4_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YYP4_OA$JID)]
all_individuals_prepatch_diffusions$tada4 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada4), all_individuals_prepatch_diffusions$tada4, visit_data_patch_YYP4_OA$patch_active + 3600)

#Order of acquisition of YYP5
visit_data_patch_YYP5_OA <- subset(visit_data_patch, visit_data_patch$position == "YYP5")
visit_data_patch_YYP5_OA <- visit_data_patch_YYP5_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YYP5_OA <- merge(visit_data_patch_YYP5_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YYP5_OA <- subset(visit_data_patch_YYP5_OA, !is.na(visit_data_patch_YYP5_OA$start))
visit_data_patch_YYP5_OA <- visit_data_patch_YYP5_OA %>% arrange(start) 
visit_data_patch_YYP5_OA$OA <- ifelse(is.na(visit_data_patch_YYP5_OA$start), NA, row.names(visit_data_patch_YYP5_OA))  
#visit_data_patch_YYP5_OA$OA <- row.names(visit_data_patch_YYP5_OA)  
#visit_data_patch_YYP5_OA <- visit_data_patch_YYP5_OA %>% arrange(JID) 
visit_data_patch_YYP5_OA$OA <- as.numeric(visit_data_patch_YYP5_OA$OA)
visit_data_patch_YYP5_OA$am_number <- am_ind$number[match(visit_data_patch_YYP5_OA$JID, am_ind$JID)]
visit_data_patch_YYP5_OA$patch_deployed <- ymd_hms("2024-06-12 06:26:00")
visit_data_patch_YYP5_OA$patch_removed <- ymd_hms("2024-06-17 07:35:00")
visit_data_patch_YYP5_OA$patch_last_read <- ymd_hms("2024-06-17 07:35:00") #feeder functional upon arrival
visit_data_patch_YYP5_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YYP5_OA$patch_deployed), ymd_hms(visit_data_patch_YYP5_OA$patch_last_read))
visit_data_patch_YYP5_OA$patch_active <- as.duration(visit_data_patch_YYP5_OA$patch_interval)
visit_data_patch_YYP5_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YYP5_OA$patch_deployed), ymd_hms(visit_data_patch_YYP5_OA$start))
visit_data_patch_YYP5_OA$time_tada <- time_length(visit_data_patch_YYP5_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada5 <- visit_data_patch_YYP5_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YYP5_OA$JID)]
all_individuals_prepatch_diffusions$tada5 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada5), all_individuals_prepatch_diffusions$tada5, visit_data_patch_YYP5_OA$patch_active + 3600)

#Order of acquisition of YZP1
visit_data_patch_YZP1_OA <- subset(visit_data_patch, visit_data_patch$position == "YZP1")
visit_data_patch_YZP1_OA <- visit_data_patch_YZP1_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YZP1_OA <- merge(visit_data_patch_YZP1_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YZP1_OA <- subset(visit_data_patch_YZP1_OA, !is.na(visit_data_patch_YZP1_OA$start))
visit_data_patch_YZP1_OA <- visit_data_patch_YZP1_OA %>% arrange(start) 
visit_data_patch_YZP1_OA$OA <- ifelse(is.na(visit_data_patch_YZP1_OA$start), NA, row.names(visit_data_patch_YZP1_OA))  
#visit_data_patch_YZP1_OA$OA <- row.names(visit_data_patch_YZP1_OA)  
#visit_data_patch_YZP1_OA <- visit_data_patch_YZP1_OA %>% arrange(JID) 
visit_data_patch_YZP1_OA$OA <- as.numeric(visit_data_patch_YZP1_OA$OA)
visit_data_patch_YZP1_OA$am_number <- am_ind$number[match(visit_data_patch_YZP1_OA$JID, am_ind$JID)]
visit_data_patch_YZP1_OA$patch_deployed <- ymd_hms("2024-05-22 05:53:00")
visit_data_patch_YZP1_OA$patch_removed <- ymd_hms("2024-05-27 17:05:00")
visit_data_patch_YZP1_OA$patch_last_read <- ymd_hms("2024-05-27 17:05:00") #feeder functional upon arrival
visit_data_patch_YZP1_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YZP1_OA$patch_deployed), ymd_hms(visit_data_patch_YZP1_OA$patch_last_read))
visit_data_patch_YZP1_OA$patch_active <- as.duration(visit_data_patch_YZP1_OA$patch_interval)
visit_data_patch_YZP1_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YZP1_OA$patch_deployed), ymd_hms(visit_data_patch_YZP1_OA$start))
visit_data_patch_YZP1_OA$time_tada <- time_length(visit_data_patch_YZP1_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada6 <- visit_data_patch_YZP1_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YZP1_OA$JID)]
all_individuals_prepatch_diffusions$tada6 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada6), all_individuals_prepatch_diffusions$tada6, visit_data_patch_YZP1_OA$patch_active + 3600)

#Order of acquisition of YZP2
visit_data_patch_YZP2_OA <- subset(visit_data_patch, visit_data_patch$position == "YZP2")
visit_data_patch_YZP2_OA <- visit_data_patch_YZP2_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YZP2_OA <- merge(visit_data_patch_YZP2_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YZP2_OA <- subset(visit_data_patch_YZP2_OA, !is.na(visit_data_patch_YZP2_OA$start))
visit_data_patch_YZP2_OA <- visit_data_patch_YZP2_OA %>% arrange(start) 
visit_data_patch_YZP2_OA$OA <- ifelse(is.na(visit_data_patch_YZP2_OA$start), NA, row.names(visit_data_patch_YZP2_OA))  
#visit_data_patch_YZP2_OA$OA <- row.names(visit_data_patch_YZP2_OA)  
#visit_data_patch_YZP2_OA <- visit_data_patch_YZP2_OA %>% arrange(JID) 
visit_data_patch_YZP2_OA$OA <- as.numeric(visit_data_patch_YZP2_OA$OA)
visit_data_patch_YZP2_OA$am_number <- am_ind$number[match(visit_data_patch_YZP2_OA$JID, am_ind$JID)]
visit_data_patch_YZP2_OA$patch_deployed <- ymd_hms("2024-05-28 06:23:00")
visit_data_patch_YZP2_OA$patch_removed <- ymd_hms("2024-06-02 06:29:00")
visit_data_patch_YZP2_OA$patch_last_read <- ymd_hms("2024-06-02 06:29:00") #feeder functional upon arrival
visit_data_patch_YZP2_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YZP2_OA$patch_deployed), ymd_hms(visit_data_patch_YZP2_OA$patch_last_read))
visit_data_patch_YZP2_OA$patch_active <- as.duration(visit_data_patch_YZP2_OA$patch_interval)
visit_data_patch_YZP2_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YZP2_OA$patch_deployed), ymd_hms(visit_data_patch_YZP2_OA$start))
visit_data_patch_YZP2_OA$time_tada <- time_length(visit_data_patch_YZP2_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada7 <- visit_data_patch_YZP2_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YZP2_OA$JID)]
all_individuals_prepatch_diffusions$tada7 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada7), all_individuals_prepatch_diffusions$tada7, visit_data_patch_YZP2_OA$patch_active + 3600)

#Order of acquisition of YZP3 (not discovered)
visit_data_patch_YZP3_OA <- subset(visit_data_patch, visit_data_patch$position == "YZP3")
visit_data_patch_YZP3_OA <- add_row(visit_data_patch_YZP3_OA)
visit_data_patch_YZP3_OA <- visit_data_patch_YZP3_OA %>%
  group_by(JID) %>% 
  slice(1:1)
#visit_data_patch_YZP3_OA <- merge(visit_data_patch_YZP3_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
#visit_data_patch_YZP3_OA <- subset(visit_data_patch_YZP3_OA, !is.na(visit_data_patch_YZP3_OA$start))
visit_data_patch_YZP3_OA <- visit_data_patch_YZP3_OA %>% arrange(start) 
visit_data_patch_YZP3_OA$OA <- ifelse(is.na(visit_data_patch_YZP3_OA$start), NA, row.names(visit_data_patch_YZP3_OA))  
#visit_data_patch_YZP3_OA$OA <- row.names(visit_data_patch_YZP3_OA)  
#visit_data_patch_YZP3_OA <- visit_data_patch_YZP3_OA %>% arrange(JID) 
visit_data_patch_YZP3_OA$OA <- as.numeric(visit_data_patch_YZP3_OA$OA)
visit_data_patch_YZP3_OA$am_number <- am_ind$number[match(visit_data_patch_YZP3_OA$JID, am_ind$JID)]
visit_data_patch_YZP3_OA$patch_deployed <- ymd_hms("2024-06-02 07:51:00")
visit_data_patch_YZP3_OA$patch_removed <- ymd_hms("2024-06-06 18:06:00")
visit_data_patch_YZP3_OA$patch_last_read <- ymd_hms("2024-06-06 18:06:00") #feeder functional upon arrival
visit_data_patch_YZP3_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YZP3_OA$patch_deployed), ymd_hms(visit_data_patch_YZP3_OA$patch_last_read))
visit_data_patch_YZP3_OA$patch_active <- as.duration(visit_data_patch_YZP3_OA$patch_interval)
visit_data_patch_YZP3_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YZP3_OA$patch_deployed), ymd_hms(visit_data_patch_YZP3_OA$start))
visit_data_patch_YZP3_OA$time_tada <- time_length(visit_data_patch_YZP3_OA$interval_tada, "second")

#Order of acquisition of YZP4
visit_data_patch_YZP4_OA <- subset(visit_data_patch, visit_data_patch$position == "YZP4")
visit_data_patch_YZP4_OA <- visit_data_patch_YZP4_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YZP4_OA <- merge(visit_data_patch_YZP4_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YZP4_OA <- subset(visit_data_patch_YZP4_OA, !is.na(visit_data_patch_YZP4_OA$start))
visit_data_patch_YZP4_OA <- visit_data_patch_YZP4_OA %>% arrange(start) 
visit_data_patch_YZP4_OA$OA <- ifelse(is.na(visit_data_patch_YZP4_OA$start), NA, row.names(visit_data_patch_YZP4_OA))  
#visit_data_patch_YZP4_OA$OA <- row.names(visit_data_patch_YZP4_OA)  
#visit_data_patch_YZP4_OA <- visit_data_patch_YZP4_OA %>% arrange(JID) 
visit_data_patch_YZP4_OA$OA <- as.numeric(visit_data_patch_YZP4_OA$OA)
visit_data_patch_YZP4_OA$am_number <- am_ind$number[match(visit_data_patch_YZP4_OA$JID, am_ind$JID)]
visit_data_patch_YZP4_OA$patch_deployed <- ymd_hms("2024-06-07 06:43:00")
visit_data_patch_YZP4_OA$patch_removed <- ymd_hms("2024-06-11 21:22:00")
visit_data_patch_YZP4_OA$patch_last_read <- ymd_hms("2024-06-11 21:22:00") #feeder functional upon arrival
visit_data_patch_YZP4_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YZP4_OA$patch_deployed), ymd_hms(visit_data_patch_YZP4_OA$patch_last_read))
visit_data_patch_YZP4_OA$patch_active <- as.duration(visit_data_patch_YZP4_OA$patch_interval)
visit_data_patch_YZP4_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YZP4_OA$patch_deployed), ymd_hms(visit_data_patch_YZP4_OA$start))
visit_data_patch_YZP4_OA$time_tada <- time_length(visit_data_patch_YZP4_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada8 <- visit_data_patch_YZP4_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YZP4_OA$JID)]
all_individuals_prepatch_diffusions$tada8 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada8), all_individuals_prepatch_diffusions$tada8, visit_data_patch_YZP4_OA$patch_active + 3600)

#Order of acquisition of YZP5
visit_data_patch_YZP5_OA <- subset(visit_data_patch, visit_data_patch$position == "YZP5")
visit_data_patch_YZP5_OA <- visit_data_patch_YZP5_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_YZP5_OA <- merge(visit_data_patch_YZP5_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_YZP5_OA <- subset(visit_data_patch_YZP5_OA, !is.na(visit_data_patch_YZP5_OA$start))
visit_data_patch_YZP5_OA <- visit_data_patch_YZP5_OA %>% arrange(start) 
visit_data_patch_YZP5_OA$OA <- ifelse(is.na(visit_data_patch_YZP5_OA$start), NA, row.names(visit_data_patch_YZP5_OA))  
#visit_data_patch_YZP5_OA$OA <- row.names(visit_data_patch_YZP5_OA)  
#visit_data_patch_YZP5_OA <- visit_data_patch_YZP5_OA %>% arrange(JID) 
visit_data_patch_YZP5_OA$OA <- as.numeric(visit_data_patch_YZP5_OA$OA)
visit_data_patch_YZP5_OA$am_number <- am_ind$number[match(visit_data_patch_YZP5_OA$JID, am_ind$JID)]
visit_data_patch_YZP5_OA$patch_deployed <- ymd_hms("2024-06-12 06:57:00")
visit_data_patch_YZP5_OA$patch_removed <- ymd_hms("2024-06-17 07:49:00")
visit_data_patch_YZP5_OA$patch_last_read <- ymd_hms("2024-06-17 07:49:00") #feeder functional upon arrival
visit_data_patch_YZP5_OA$patch_interval <- interval(ymd_hms(visit_data_patch_YZP5_OA$patch_deployed), ymd_hms(visit_data_patch_YZP5_OA$patch_last_read))
visit_data_patch_YZP5_OA$patch_active <- as.duration(visit_data_patch_YZP5_OA$patch_interval)
visit_data_patch_YZP5_OA$interval_tada <- interval(ymd_hms(visit_data_patch_YZP5_OA$patch_deployed), ymd_hms(visit_data_patch_YZP5_OA$start))
visit_data_patch_YZP5_OA$time_tada <- time_length(visit_data_patch_YZP5_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada9 <- visit_data_patch_YZP5_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_YZP5_OA$JID)]
all_individuals_prepatch_diffusions$tada9 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada9), all_individuals_prepatch_diffusions$tada9, visit_data_patch_YZP5_OA$patch_active + 3600)

#Order of acquisition of ZZP1
visit_data_patch_ZZP1_OA <- subset(visit_data_patch, visit_data_patch$position == "ZZP1")
visit_data_patch_ZZP1_OA <- visit_data_patch_ZZP1_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_ZZP1_OA <- merge(visit_data_patch_ZZP1_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_ZZP1_OA <- subset(visit_data_patch_ZZP1_OA, !is.na(visit_data_patch_ZZP1_OA$start))
visit_data_patch_ZZP1_OA <- visit_data_patch_ZZP1_OA %>% arrange(start) 
visit_data_patch_ZZP1_OA$OA <- ifelse(is.na(visit_data_patch_ZZP1_OA$start), NA, row.names(visit_data_patch_ZZP1_OA))  
#visit_data_patch_ZZP1_OA$OA <- row.names(visit_data_patch_ZZP1_OA)  
#visit_data_patch_ZZP1_OA <- visit_data_patch_ZZP1_OA %>% arrange(JID) 
visit_data_patch_ZZP1_OA$OA <- as.numeric(visit_data_patch_ZZP1_OA$OA)
visit_data_patch_ZZP1_OA$am_number <- am_ind$number[match(visit_data_patch_ZZP1_OA$JID, am_ind$JID)]
visit_data_patch_ZZP1_OA$patch_deployed <- ymd_hms("2024-05-22 05:31:00")
visit_data_patch_ZZP1_OA$patch_removed <- ymd_hms("2024-05-27 16:01:00")
visit_data_patch_ZZP1_OA$patch_last_read <- ymd_hms("2024-05-27 16:01:00") #feeder functional upon arrival
visit_data_patch_ZZP1_OA$patch_interval <- interval(ymd_hms(visit_data_patch_ZZP1_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP1_OA$patch_last_read))
visit_data_patch_ZZP1_OA$patch_active <- as.duration(visit_data_patch_ZZP1_OA$patch_interval)
visit_data_patch_ZZP1_OA$interval_tada <- interval(ymd_hms(visit_data_patch_ZZP1_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP1_OA$start))
visit_data_patch_ZZP1_OA$time_tada <- time_length(visit_data_patch_ZZP1_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada10 <- visit_data_patch_ZZP1_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_ZZP1_OA$JID)]
all_individuals_prepatch_diffusions$tada10 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada10), all_individuals_prepatch_diffusions$tada10, visit_data_patch_ZZP1_OA$patch_active + 3600)

#Order of acquisition of ZZP2
visit_data_patch_ZZP2_OA <- subset(visit_data_patch, visit_data_patch$position == "ZZP2")
visit_data_patch_ZZP2_OA <- visit_data_patch_ZZP2_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_ZZP2_OA <- merge(visit_data_patch_ZZP2_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_ZZP2_OA <- subset(visit_data_patch_ZZP2_OA, !is.na(visit_data_patch_ZZP2_OA$start))
visit_data_patch_ZZP2_OA <- visit_data_patch_ZZP2_OA %>% arrange(start) 
visit_data_patch_ZZP2_OA$OA <- ifelse(is.na(visit_data_patch_ZZP2_OA$start), NA, row.names(visit_data_patch_ZZP2_OA))  
#visit_data_patch_ZZP2_OA$OA <- row.names(visit_data_patch_ZZP2_OA)  
#visit_data_patch_ZZP2_OA <- visit_data_patch_ZZP2_OA %>% arrange(JID) 
visit_data_patch_ZZP2_OA$OA <- as.numeric(visit_data_patch_ZZP2_OA$OA)
visit_data_patch_ZZP2_OA$am_number <- am_ind$number[match(visit_data_patch_ZZP2_OA$JID, am_ind$JID)]
visit_data_patch_ZZP2_OA$patch_deployed <- ymd_hms("2024-05-28 06:08:00")
visit_data_patch_ZZP2_OA$patch_removed <- ymd_hms("2024-05-30 09:00:00") #feeder had to be removed
visit_data_patch_ZZP2_OA$patch_last_read <- ymd_hms("2024-05-30 09:00:00") #feeder had to be removed
visit_data_patch_ZZP2_OA$patch_interval <- interval(ymd_hms(visit_data_patch_ZZP2_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP2_OA$patch_last_read))
visit_data_patch_ZZP2_OA$patch_active <- as.duration(visit_data_patch_ZZP2_OA$patch_interval)
visit_data_patch_ZZP2_OA$interval_tada <- interval(ymd_hms(visit_data_patch_ZZP2_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP2_OA$start))
visit_data_patch_ZZP2_OA$time_tada <- time_length(visit_data_patch_ZZP2_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada11 <- visit_data_patch_ZZP2_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_ZZP2_OA$JID)]
all_individuals_prepatch_diffusions$tada11 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada11), all_individuals_prepatch_diffusions$tada11, visit_data_patch_ZZP2_OA$patch_active + 3600)

#Order of acquisition of ZZP3
visit_data_patch_ZZP3_OA <- subset(visit_data_patch, visit_data_patch$position == "ZZP3")
visit_data_patch_ZZP3_OA <- visit_data_patch_ZZP3_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_ZZP3_OA <- merge(visit_data_patch_ZZP3_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_ZZP3_OA <- subset(visit_data_patch_ZZP3_OA, !is.na(visit_data_patch_ZZP3_OA$start))
visit_data_patch_ZZP3_OA <- visit_data_patch_ZZP3_OA %>% arrange(start) 
visit_data_patch_ZZP3_OA$OA <- ifelse(is.na(visit_data_patch_ZZP3_OA$start), NA, row.names(visit_data_patch_ZZP3_OA))  
#visit_data_patch_ZZP3_OA$OA <- row.names(visit_data_patch_ZZP3_OA)  
#visit_data_patch_ZZP3_OA <- visit_data_patch_ZZP3_OA %>% arrange(JID) 
visit_data_patch_ZZP3_OA$OA <- as.numeric(visit_data_patch_ZZP3_OA$OA)
visit_data_patch_ZZP3_OA$am_number <- am_ind$number[match(visit_data_patch_ZZP3_OA$JID, am_ind$JID)]
visit_data_patch_ZZP3_OA$patch_deployed <- ymd_hms("2024-06-02 07:31:00")
visit_data_patch_ZZP3_OA$patch_removed <- ymd_hms("2024-06-06 18:16:00")
visit_data_patch_ZZP3_OA$patch_last_read <- ymd_hms("2024-06-06 18:16:00") #feeder functional upon arrival
visit_data_patch_ZZP3_OA$patch_interval <- interval(ymd_hms(visit_data_patch_ZZP3_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP3_OA$patch_last_read))
visit_data_patch_ZZP3_OA$patch_active <- as.duration(visit_data_patch_ZZP3_OA$patch_interval)
visit_data_patch_ZZP3_OA$interval_tada <- interval(ymd_hms(visit_data_patch_ZZP3_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP3_OA$start))
visit_data_patch_ZZP3_OA$time_tada <- time_length(visit_data_patch_ZZP3_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada12 <- visit_data_patch_ZZP3_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_ZZP3_OA$JID)]
all_individuals_prepatch_diffusions$tada12 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada12), all_individuals_prepatch_diffusions$tada12, visit_data_patch_ZZP3_OA$patch_active + 3600)

#Order of acquisition of ZZP4
visit_data_patch_ZZP4_OA <- subset(visit_data_patch, visit_data_patch$position == "ZZP4")
visit_data_patch_ZZP4_OA <- visit_data_patch_ZZP4_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_ZZP4_OA <- merge(visit_data_patch_ZZP4_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_ZZP4_OA <- subset(visit_data_patch_ZZP4_OA, !is.na(visit_data_patch_ZZP4_OA$start))
visit_data_patch_ZZP4_OA <- visit_data_patch_ZZP4_OA %>% arrange(start) 
visit_data_patch_ZZP4_OA$OA <- ifelse(is.na(visit_data_patch_ZZP4_OA$start), NA, row.names(visit_data_patch_ZZP4_OA))  
#visit_data_patch_ZZP4_OA$OA <- row.names(visit_data_patch_ZZP4_OA)  
#visit_data_patch_ZZP4_OA <- visit_data_patch_ZZP4_OA %>% arrange(JID) 
visit_data_patch_ZZP4_OA$OA <- as.numeric(visit_data_patch_ZZP4_OA$OA)
visit_data_patch_ZZP4_OA$am_number <- am_ind$number[match(visit_data_patch_ZZP4_OA$JID, am_ind$JID)]
visit_data_patch_ZZP4_OA$patch_deployed <- ymd_hms("2024-06-07 06:20:00")
visit_data_patch_ZZP4_OA$patch_removed <- ymd_hms("2024-06-11 20:59:00")
visit_data_patch_ZZP4_OA$patch_last_read <- ymd_hms("2024-06-11 20:59:00") #feeder functional upon arrival
visit_data_patch_ZZP4_OA$patch_interval <- interval(ymd_hms(visit_data_patch_ZZP4_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP4_OA$patch_last_read))
visit_data_patch_ZZP4_OA$patch_active <- as.duration(visit_data_patch_ZZP4_OA$patch_interval)
visit_data_patch_ZZP4_OA$interval_tada <- interval(ymd_hms(visit_data_patch_ZZP4_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP4_OA$start))
visit_data_patch_ZZP4_OA$time_tada <- time_length(visit_data_patch_ZZP4_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada13 <- visit_data_patch_ZZP4_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_ZZP4_OA$JID)]
all_individuals_prepatch_diffusions$tada13 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada13), all_individuals_prepatch_diffusions$tada13, visit_data_patch_ZZP4_OA$patch_active + 3600)

#Order of acquisition of ZZP5
visit_data_patch_ZZP5_OA <- subset(visit_data_patch, visit_data_patch$position == "ZZP5")
visit_data_patch_ZZP5_OA <- visit_data_patch_ZZP5_OA %>%
  group_by(JID) %>% 
  slice(1:1)
visit_data_patch_ZZP5_OA <- merge(visit_data_patch_ZZP5_OA, all_individuals_prepatch_JID, by = "JID", all.y = TRUE)
visit_data_patch_ZZP5_OA <- subset(visit_data_patch_ZZP5_OA, !is.na(visit_data_patch_ZZP5_OA$start))
visit_data_patch_ZZP5_OA <- visit_data_patch_ZZP5_OA %>% arrange(start) 
visit_data_patch_ZZP5_OA$OA <- ifelse(is.na(visit_data_patch_ZZP5_OA$start), NA, row.names(visit_data_patch_ZZP5_OA))  
#visit_data_patch_ZZP5_OA$OA <- row.names(visit_data_patch_ZZP5_OA)  
#visit_data_patch_ZZP5_OA <- visit_data_patch_ZZP5_OA %>% arrange(JID) 
visit_data_patch_ZZP5_OA$OA <- as.numeric(visit_data_patch_ZZP5_OA$OA)
visit_data_patch_ZZP5_OA$am_number <- am_ind$number[match(visit_data_patch_ZZP5_OA$JID, am_ind$JID)]
visit_data_patch_ZZP5_OA$patch_deployed <- ymd_hms("2024-06-12 05:54:00")
visit_data_patch_ZZP5_OA$patch_removed <- ymd_hms("2024-06-17 08:01:00")
visit_data_patch_ZZP5_OA$patch_last_read <- ymd_hms("2024-06-13 17:42:00") #feeder functional upon arrival
visit_data_patch_ZZP5_OA$patch_interval <- interval(ymd_hms(visit_data_patch_ZZP5_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP5_OA$patch_last_read))
visit_data_patch_ZZP5_OA$patch_active <- as.duration(visit_data_patch_ZZP5_OA$patch_interval)
visit_data_patch_ZZP5_OA$interval_tada <- interval(ymd_hms(visit_data_patch_ZZP5_OA$patch_deployed), ymd_hms(visit_data_patch_ZZP5_OA$start))
visit_data_patch_ZZP5_OA$time_tada <- time_length(visit_data_patch_ZZP5_OA$interval_tada, "second")

all_individuals_prepatch_diffusions$tada14 <- visit_data_patch_ZZP5_OA$time_tada[match(all_individuals_prepatch_diffusions$JID, visit_data_patch_ZZP5_OA$JID)]
all_individuals_prepatch_diffusions$tada14 <- ifelse(!is.na(all_individuals_prepatch_diffusions$tada14), all_individuals_prepatch_diffusions$tada14, visit_data_patch_ZZP5_OA$patch_active + 3600)

#Dataset combining all patch discoveries
visit_data_patch_OA <- rbind(visit_data_patch_YYP1_OA, visit_data_patch_YYP2_OA, visit_data_patch_YYP3_OA, visit_data_patch_YYP4_OA, visit_data_patch_YYP5_OA,
                             visit_data_patch_YZP1_OA, visit_data_patch_YZP2_OA, visit_data_patch_YZP4_OA, visit_data_patch_YZP5_OA,
                             visit_data_patch_ZZP1_OA, visit_data_patch_ZZP2_OA, visit_data_patch_ZZP3_OA, visit_data_patch_ZZP4_OA, visit_data_patch_ZZP5_OA)

table(visit_data_patch_OA$feeder_patch)
table(visit_data_patch_OA$JID)

#Mean order of arrival of individuals that discovered at least one patch 
visit_data_patch_OA <- as.data.frame(aggregate(visit_data_patch_OA$OA, by = list(visit_data_patch_OA$JID), FUN = "mean", na.rm = TRUE))
visit_data_patch_OA$JID <- visit_data_patch_OA$Group.1
visit_data_patch_OA$OA <- visit_data_patch_OA$x
visit_data_patch_OA <- subset(visit_data_patch_OA, select = -c(Group.1, x))

#Dataset combining all patch discoveries
visit_data_patch_all_OA <- rbind(visit_data_patch_YYP1_OA, visit_data_patch_YYP2_OA,
                                 visit_data_patch_YYP3_OA, visit_data_patch_YYP4_OA,
                                 visit_data_patch_YYP5_OA, visit_data_patch_YZP1_OA,
                                 visit_data_patch_YZP2_OA, visit_data_patch_YZP4_OA,
                                 visit_data_patch_YZP5_OA, visit_data_patch_ZZP1_OA, 
                                 visit_data_patch_ZZP2_OA, visit_data_patch_ZZP3_OA, 
                                 visit_data_patch_ZZP4_OA, visit_data_patch_ZZP5_OA)

visit_data_patch_all_OA$JID_patch <- paste(visit_data_patch_all_OA$JID, visit_data_patch_all_OA$position, sep = "_")

visit_data_patch_all_OA$degree_z <- all_individuals_prepatch_sub$degree_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$strength_z <- all_individuals_prepatch_sub$strength_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$eigenvector_z <- all_individuals_prepatch_sub$eigenvector_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$betweenness_z <- all_individuals_prepatch_sub$betweenness_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$closeness_z <- all_individuals_prepatch_sub$closeness_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$social_diff_z <- all_individuals_prepatch_sub$social_diff_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$transitivity_z <- all_individuals_prepatch_sub$transitivity_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]

visit_data_patch_all_OA$age_z <- all_individuals_prepatch_sub$age_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$sex <- all_individuals_prepatch_sub$sex[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$sex_cat <- all_individuals_prepatch_sub$sex_cat[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$body_cond_z <- all_individuals_prepatch_sub$body_cond_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$tarsus_z <- all_individuals_prepatch_sub$tarsus_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$visit_number <- all_individuals_prepatch_sub$visit_number[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$visit_number_z <- all_individuals_prepatch_sub$visit_number_z[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$feeder_number_binary <- all_individuals_prepatch_sub$feeder_number_binary[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]
visit_data_patch_all_OA$inter_sites <- all_individuals_prepatch_sub$inter_sites[match(visit_data_patch_all_OA$JID, all_individuals_prepatch_sub$JID)]

visit_data_patch_all_OA$time_tada_hr <- visit_data_patch_all_OA$time_tada/60/60

patch_ID <- c("YYP1", "YYP2", "YYP3", "YYP4", "YYP5", "YZP1", "YZP2", "YZP3", "YZP4", "YZP5", "ZZP1", "ZZP2", "ZZP3", "ZZP4", "ZZP5")
patch_active <- c(1:15)
all_patches <- as.data.frame(cbind(patch_ID, patch_active))
all_patches$patch_active <- visit_data_patch_all_OA$patch_active[match(all_patches$patch_ID, visit_data_patch_all_OA$position)]
all_patches[8,2] <- visit_data_patch_YZP3_OA$patch_active
mean(all_patches$patch_active)/60/60/24
median(all_patches$patch_active)/60/60/24
sd(all_patches$patch_active)/60/60/24

#Data set combining all individuals and all patches, including individuals that discovered patches but were not in the social network
all_individuals_all_patches <- crossing(JID = all_individuals_prepatch_sub$JID, patch_ID = all_patches$patch_ID)
all_individuals_all_patches$JID_patch <- paste(all_individuals_all_patches$JID, all_individuals_all_patches$patch_ID, sep = "_")

all_individuals_all_patches <- merge(visit_data_patch_all_OA, all_individuals_all_patches, by = "JID_patch", all = TRUE)
all_individuals_all_patches$JID.x <- NULL
all_individuals_all_patches$JID <- all_individuals_all_patches$JID.y
all_individuals_all_patches$JID.y <- NULL
all_individuals_all_patches$discovered_binary <- ifelse(is.na(all_individuals_all_patches$time_tada), 0, 1)

table(all_individuals_all_patches$JID_patch)


#Event data ----

#Event data for STBayes NBDA 

#167 individuals * 14 patches (1 never discovered)
event_data <- pivot_longer(all_individuals_prepatch_diffusions[, c("JID", "tada1", "tada2", "tada3", "tada4", "tada5", "tada6", "tada7", "tada8", "tada9", "tada10", "tada11", "tada12", "tada13", "tada14")], cols = c("tada1", "tada2", "tada3", "tada4", "tada5", "tada6", "tada7", "tada8", "tada9", "tada10", "tada11", "tada12", "tada13", "tada14"), names_to = "trial", values_to = "time")
event_data   <- event_data  %>% rename(id = JID)

#Add end of each trial (hours after set-up)
event_data$t_end <- NA
event_data$t_end[event_data$trial == "tada1"] <- 86.9833
event_data$t_end[event_data$trial == "tada2"] <- 119.9
event_data$t_end[event_data$trial == "tada3"] <- 106.9
event_data$t_end[event_data$trial == "tada4"] <- 112.15
event_data$t_end[event_data$trial == "tada5"] <- 121.15
event_data$t_end[event_data$trial == "tada6"] <- 131.2
event_data$t_end[event_data$trial == "tada7"] <- 120.1
event_data$t_end[event_data$trial == "tada8"] <- 110.65
event_data$t_end[event_data$trial == "tada9"] <- 120.866667
event_data$t_end[event_data$trial == "tada10"] <- 130.5
event_data$t_end[event_data$trial == "tada11"] <- 50.866667
event_data$t_end[event_data$trial == "tada12"] <- 106.75
event_data$t_end[event_data$trial == "tada13"] <- 110.65
event_data$t_end[event_data$trial == "tada14"] <- 35.8

#cTADA so integer not relevant
#event_data$time <- as.integer(event_data$time)
#event_data$t_end <- as.integer(event_data$t_end)

#Edge list ----

#Create edge list for STBayes NBDA

#Edge list from graph object (prepatch social network)
g_dyads_prepatch <- igraph::ends(g, igraph::E(g), names=TRUE)
g_edge_weights_prepatch <- igraph::E(g)$weight

dyad_edge_prepatch <- data.frame(matrix(NA, nrow = 1515, ncol = 2))
dyad_edge_prepatch$X1 <- g_dyads_prepatch
dyad_edge_prepatch <- dyad_edge_prepatch %>% rename(trial = X2)
dyad_edge_prepatch$from <- dyad_edge_prepatch$X1[,1]
dyad_edge_prepatch$to <- dyad_edge_prepatch$X1[,2]
dyad_edge_prepatch$weight <- g_edge_weights_prepatch 
dyad_edge_prepatch <- dyad_edge_prepatch[, c("from", "to", "trial", "weight")]

edge_list <- dyad_edge_prepatch

#Add self-edges for individuals with 0 connections
edge_list[1516,1] <- "J1057"
edge_list[1517,1] <- "J2701"
edge_list[1518,1] <- "J4511"
edge_list[1519,1] <- "J4605"
edge_list[1520,1] <- "J558"

edge_list[1516,2] <- "J1057"
edge_list[1517,2] <- "J2701"
edge_list[1518,2] <- "J4511"
edge_list[1519,2] <- "J4605"
edge_list[1520,2] <- "J558"

edge_list_1 <- edge_list
edge_list_2 <- edge_list
edge_list_3 <- edge_list
edge_list_4 <- edge_list
edge_list_5 <- edge_list
edge_list_6 <- edge_list
edge_list_7 <- edge_list
edge_list_8 <- edge_list
edge_list_9 <- edge_list
edge_list_10 <- edge_list
edge_list_11 <- edge_list
edge_list_12 <- edge_list
edge_list_13 <- edge_list
edge_list_14 <- edge_list

edge_list_1$trial <- "tada1"
edge_list_2$trial <- "tada2"
edge_list_3$trial <- "tada3"
edge_list_4$trial <- "tada4"
edge_list_5$trial <- "tada5"
edge_list_6$trial <- "tada6"
edge_list_7$trial <- "tada7"
edge_list_8$trial <- "tada8"
edge_list_9$trial <- "tada9"
edge_list_10$trial <- "tada10"
edge_list_11$trial <- "tada11"
edge_list_12$trial <- "tada12"
edge_list_13$trial <- "tada13"
edge_list_14$trial <- "tada14"

edge_list <- rbind(edge_list_1, edge_list_2, edge_list_3, edge_list_4, edge_list_5, edge_list_6, edge_list_7, edge_list_8, edge_list_9, edge_list_10, edge_list_11, edge_list_12, edge_list_13, edge_list_14)

#Data list
data_list <- import_user_STb(event_data = event_data, 
                             networks = edge_list,
                             network_type = "undirected")

#First models ----

#Models without transmission weights, which are added in the final models

#Social model 
model_full <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = c("continuous"))

#Fit social model
full_fit <- fit_STb(data_list,
                    model_full,
                    parallel_chains = 4,
                    chains = 4,
                    cores = 4,
                    iter = 4000,
                    refresh=1000
)

STb_save(full_fit, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit, digits = 3)

data_list <- import_user_STb(event_data = event_data, 
                             networks = edge_list,
                             network_type = "undirected")

#Asocial model
model_asoc = generate_STb_model(data_list, model_type="asocial")

#Fit asocial model 
asocial_fit = fit_STb(data_list,
                      model_asoc,
                      parallel_chains = 4,
                      chains = 4,
                      cores = 4,
                      iter = 4000,
                      refresh=1000)

#Compare models
loo_output = STb_compare(full_fit, asocial_fit, method="loo-psis")

print(loo_output$comparison, simplify = FALSE)


#ILVs ----

#Individual-level variables

ILV_c <- data.frame(
  id = all_individuals_prepatch_JID$JID,
  age = age, 
  sex = sex, 
  tarsus = tarsus,
  body_cond = body_cond,
  feederuse = feederuse
)

data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  ILV_c = ILV_c,
  ILVi = c("tarsus", "body_cond", "feederuse", "age", "sex"),
  ILVs = c( "age", "sex")
)

model_full_ilv <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = c("continuous"))

full_fit_ilv <- fit_STb(data_list,
                        model_full_ilv,
                        parallel_chains = 4,
                        chains = 4,
                        cores = 4,
                        iter = 4000,
                        refresh=1000
)

STb_save(full_fit_ilv, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_ilv, digits = 3)

loo_output = STb_compare(full_fit_ilv, full_fit, asocial_fit, method="loo-psis")

print(loo_output$comparison, simplify = FALSE)

#Transmission weights ----

#Adding transmission weights
all_individuals_prepatch_diffusions <- cbind(all_individuals_prepatch_diffusions$JID, all_individuals_prepatch_diffusions[,2:15]/3600)
all_individuals_prepatch_diffusions <- all_individuals_prepatch_diffusions%>%
  rename(JID = `all_individuals_prepatch_diffusions$JID`)

#Visit number per individual per patch
all_individuals_prepatch_diffusions$tw1 <- visit_data_patch_visit_number[, c("YYP1")]
all_individuals_prepatch_diffusions$tw2 <- visit_data_patch_visit_number[, c("YYP2")]
all_individuals_prepatch_diffusions$tw3 <- visit_data_patch_visit_number[, c("YYP3")]
all_individuals_prepatch_diffusions$tw4 <- visit_data_patch_visit_number[, c("YYP4")]
all_individuals_prepatch_diffusions$tw5 <- visit_data_patch_visit_number[, c("YYP5")]
all_individuals_prepatch_diffusions$tw6 <- visit_data_patch_visit_number[, c("YZP1")]
all_individuals_prepatch_diffusions$tw7 <- visit_data_patch_visit_number[, c("YZP2")]
all_individuals_prepatch_diffusions$tw8 <- visit_data_patch_visit_number[, c("YZP4")]
all_individuals_prepatch_diffusions$tw9 <- visit_data_patch_visit_number[, c("YZP5")]
all_individuals_prepatch_diffusions$tw10 <- visit_data_patch_visit_number[, c("ZZP1")]
all_individuals_prepatch_diffusions$tw11 <- visit_data_patch_visit_number[, c("ZZP2")]
all_individuals_prepatch_diffusions$tw12 <- visit_data_patch_visit_number[, c("ZZP3")]
all_individuals_prepatch_diffusions$tw13 <- visit_data_patch_visit_number[, c("ZZP4")]
all_individuals_prepatch_diffusions$tw14 <- visit_data_patch_visit_number[, c("ZZP5")]

#Trial end (hours)
all_individuals_prepatch_diffusions$t_end1 <- 86.9833
all_individuals_prepatch_diffusions$t_end2 <- 119.9
all_individuals_prepatch_diffusions$t_end3 <- 106.9
all_individuals_prepatch_diffusions$t_end4 <- 112.15
all_individuals_prepatch_diffusions$t_end5 <- 121.15
all_individuals_prepatch_diffusions$t_end6 <- 131.2
all_individuals_prepatch_diffusions$t_end7 <- 120.1
all_individuals_prepatch_diffusions$t_end8 <- 110.65
all_individuals_prepatch_diffusions$t_end9 <- 120.866667
all_individuals_prepatch_diffusions$t_end10 <- 130.5
all_individuals_prepatch_diffusions$t_end11 <- 50.866667
all_individuals_prepatch_diffusions$t_end12 <- 106.75
all_individuals_prepatch_diffusions$t_end13 <- 110.65
all_individuals_prepatch_diffusions$t_end14 <- 35.8

#Transmission weight: visit rates for individual at each patch per hour of being knowledgeable
all_individuals_prepatch_diffusions$tww1 <- all_individuals_prepatch_diffusions$tw1/(all_individuals_prepatch_diffusions$t_end1 - all_individuals_prepatch_diffusions$tada1)
all_individuals_prepatch_diffusions$tww2 <- all_individuals_prepatch_diffusions$tw2/(all_individuals_prepatch_diffusions$t_end2 - all_individuals_prepatch_diffusions$tada2)
all_individuals_prepatch_diffusions$tww3 <- all_individuals_prepatch_diffusions$tw3/(all_individuals_prepatch_diffusions$t_end3 - all_individuals_prepatch_diffusions$tada3)
all_individuals_prepatch_diffusions$tww4 <- all_individuals_prepatch_diffusions$tw4/(all_individuals_prepatch_diffusions$t_end4 - all_individuals_prepatch_diffusions$tada4)
all_individuals_prepatch_diffusions$tww5 <- all_individuals_prepatch_diffusions$tw5/(all_individuals_prepatch_diffusions$t_end5 - all_individuals_prepatch_diffusions$tada5)
all_individuals_prepatch_diffusions$tww6 <- all_individuals_prepatch_diffusions$tw6/(all_individuals_prepatch_diffusions$t_end6 - all_individuals_prepatch_diffusions$tada6)
all_individuals_prepatch_diffusions$tww7 <- all_individuals_prepatch_diffusions$tw7/(all_individuals_prepatch_diffusions$t_end7 - all_individuals_prepatch_diffusions$tada7)
all_individuals_prepatch_diffusions$tww8 <- all_individuals_prepatch_diffusions$tw8/(all_individuals_prepatch_diffusions$t_end8 - all_individuals_prepatch_diffusions$tada8)
all_individuals_prepatch_diffusions$tww9 <- all_individuals_prepatch_diffusions$tw9/(all_individuals_prepatch_diffusions$t_end9 - all_individuals_prepatch_diffusions$tada9)
all_individuals_prepatch_diffusions$tww10 <- all_individuals_prepatch_diffusions$tw10/(all_individuals_prepatch_diffusions$t_end10 - all_individuals_prepatch_diffusions$tada10)
all_individuals_prepatch_diffusions$tww11 <- all_individuals_prepatch_diffusions$tw11/(all_individuals_prepatch_diffusions$t_end11 - all_individuals_prepatch_diffusions$tada11)
all_individuals_prepatch_diffusions$tww12 <- all_individuals_prepatch_diffusions$tw12/(all_individuals_prepatch_diffusions$t_end12 - all_individuals_prepatch_diffusions$tada12)
all_individuals_prepatch_diffusions$tww13 <- all_individuals_prepatch_diffusions$tw13/(all_individuals_prepatch_diffusions$t_end13 - all_individuals_prepatch_diffusions$tada13)
all_individuals_prepatch_diffusions$tww14 <- all_individuals_prepatch_diffusions$tw14/(all_individuals_prepatch_diffusions$t_end14 - all_individuals_prepatch_diffusions$tada14)

#Data set including transmission weights
t_weights_static <- all_individuals_prepatch_diffusions %>%
  pivot_longer(
    cols = starts_with("tada"),
    names_to = "trial",
    values_to = "tada_value"
  ) %>%
  pivot_longer(
    cols = starts_with("tww"),
    names_to = "trial_w",
    values_to = "t_weight"
  ) %>%
  # match tww columns to the corresponding tada rows
  filter(trial_w == sub("tada", "tww", trial)) %>%
  dplyr::select(JID, trial, t_weight)

t_weights_static <- t_weights_static  %>% rename(id = JID)

#Models with weights ----

#Data list
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static
)

#Full social model including transmission weights
model_full_weights <- generate_STb_model(data_list, 
                                         gq = T, 
                                         est_acqTime = F, 
                                         data_type = c("continuous"))

#Fit full social model including transmission weights
full_fit_weights <- fit_STb(data_list,
                            model_full_weights,
                            parallel_chains = 4,
                            chains = 4,
                            cores = 4,
                            iter = 4000,
                            refresh=1000
)

STb_save(full_fit_weights, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_weights, digits = 3)

model_asoc_weights = generate_STb_model(data_list, model_type="asocial")

asocial_fit_weights = fit_STb(data_list,
                              model_asoc_weights,
                              parallel_chains = 4,
                              chains = 4,
                              cores = 4,
                              iter = 4000,
                              refresh=1000)

STb_save(asocial_fit_weights, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(asocial_fit_weights, digits = 3)

#Compare social and asocial models
loo_output = STb_compare(full_fit_weights, asocial_fit_weights, method="loo-psis")

print(loo_output$comparison, simplify = FALSE)

#Models with weights and ILVs ----

#Weights and ILVs

data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static,
  ILV_c = ILV_c,
  ILVi = c("tarsus", "body_cond", "feederuse", "age", "sex"),
  ILVs = c( "age", "sex")
)

model_full_weights_ilv <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = c("continuous"))

full_fit_weights_ilv <- fit_STb(data_list,
                                model_full_ilv,
                                parallel_chains = 4,
                                chains = 4,
                                cores = 4,
                                iter = 4000,
                                refresh=1000
)

STb_save(full_fit_weights_ilv, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_weights_ilv, digits = 3)

data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static,
  ILV_c = ILV_c,
  ILVi = c("tarsus", "body_cond", "feederuse", "age", "sex"),
  ILVs = c( "age", "sex")
)

model_asoc = generate_STb_model(data_list, model_type="asocial")

asocial_fit_weights_ilv = fit_STb(data_list,
                                  model_asoc,
                                  parallel_chains = 4,
                                  chains = 4,
                                  cores = 4,
                                  iter = 4000,
                                  refresh=1000)

STb_save(full_fit_weights_ilv, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_weights_ilv, digits = 3)


#Model comparison ----

#Comparing 4 models: social and asocial with and without ILV; all with transmission weights
loo_output = STb_compare(full_fit_weights, full_fit_weights_ilv, asocial_fit_weights, asocial_fit_weights_ilv, method="loo-psis")
print(loo_output$comparison, simplify = FALSE)

#Comparing 2 models without ILVs
loo_output = STb_compare(full_fit_weights, asocial_fit_weights, method="loo-psis")
print(loo_output$comparison, simplify = FALSE)

#Comparing 2 social models 
loo_output = STb_compare(full_fit_weights, full_fit_weights_ilv, method="loo-psis")
print(loo_output$comparison, simplify = FALSE)

#Tied individuals ----
event_data_ties <- event_data %>% arrange(trial, time)
event_data <- event_data %>% arrange(id, trial, time)

threshold <- 0.25

event_data_ties <- event_data_ties %>%
  arrange(trial, time) %>%
  group_by(trial) %>%
  mutate(
    time_adj = time
  ) %>%
  group_modify(~{
    x <- .x$time
    for (i in 2:length(x)) {
      if (x[i] - x[i-1] <= threshold) {
        x[i] <- x[i-1]
      }
    }
    .x$time_adj <- x
    .x
  }) %>%
  ungroup()

#Data list
data_list <- import_user_STb(
  event_data = event_data_ties,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static
)

#Full social model including transmission weights
model_full_weights_ties <- generate_STb_model(data_list, 
                                         gq = T, 
                                         est_acqTime = F, 
                                         data_type = c("continuous"))

#Fit full social model including transmission weights
full_fit_weights_ties <- fit_STb(data_list,
                            model_full_weights_ties,
                            parallel_chains = 4,
                            chains = 4,
                            cores = 4,
                            iter = 4000,
                            refresh=1000
)

STb_save(full_fit_weights_ties, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_weights_ties, digits = 3)

model_asoc_weights_ties = generate_STb_model(data_list, model_type="asocial")

asocial_fit_weights_ties = fit_STb(data_list,
                              model_asoc_weights_ties,
                              parallel_chains = 4,
                              chains = 4,
                              cores = 4,
                              iter = 4000,
                              refresh=1000)

STb_save(asocial_fit_weights, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(asocial_fit_weights, digits = 3)

#Compare social and asocial models
loo_output = STb_compare(full_fit_weights_ties, asocial_fit_weights_ties, method="loo-psis")

print(loo_output$comparison, simplify = FALSE)

#Posterior predictive checks ----

#Cumulative diffusion curves 

# we need to store num inds per trial to refer to later
event_data <- event_data %>%
  group_by(trial) %>%
  mutate(n_trial = n())

plsot_data_obs <- event_data %>%
  filter(time > 0, time <= t_end) %>% # exclude demonstrators (time == 0) and censored (time > t_end)
  group_by(trial) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    cum_prop = row_number() / n_trial, # this denominator needs to be the number of individuals per trial
    type = "observed"
  ) %>%
  dplyr::select(trial, time, cum_prop, type) %>%
  ungroup()

# add in 0,0 starting point for diffusions w/o demonstrators
starting_points <- plot_data_obs %>%
  distinct(trial) %>%
  anti_join(
    plot_data_obs %>% filter(time == 0) %>% distinct(trial),
    by = "trial"
  ) %>%
  mutate(time = 0, cum_prop = 0, type = "observed")

plot_data_obs <- bind_rows(plot_data_obs, starting_points) %>%
  arrange(trial, time)

draws_df <- as_draws_df(full_fit_weights$draws(variables = "acquisition_time", inc_warmup = FALSE))

# pivot longer
ppc_long <- draws_df %>%
  dplyr::select(starts_with("acquisition_time[")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("trial", "ind"),
    names_pattern = "acquisition_time\\[(\\d+),(\\d+)\\]",
    values_to = "time"
  ) %>%
  mutate(
    trial = as.integer(trial),
    ind = as.integer(ind),
    draw = rep(1:(nrow(draws_df)), 
               each = length(unique(.$trial)) * length(unique(.$ind)))
  )
#> Warning: Dropping 'draws_df' class as required metadata was removed.

# thin sample for plotting
sample_idx <- sample(c(1:max(ppc_long$draw)), 100)
ppc_long <- ppc_long %>% filter(draw %in% sample_idx)

# same as before, we need a way to reference the number of individuals in each trial
ppc_long <- ppc_long %>%
  group_by(draw, trial) %>%
  mutate(n_trial = n())

# we also need to remove individuals predicted as censored, which
# have value of -1 in predicted data
ppc_long <- ppc_long %>%
  filter(time > -1)

# build cumulative curves per draw
plot_data_ppc <- ppc_long %>%
  group_by(draw, trial, time) %>%
  summarise(n = n(), n_trial = first(n_trial), .groups = "drop") %>%
  group_by(draw, trial) %>%
  arrange(time) %>%
  mutate(cum_prop = cumsum(n) / n_trial)

# add in 0,0 starting point when no demos, similar to above
starting_points_ppc <- plot_data_ppc %>%
  distinct(trial, draw) %>%
  anti_join(
    plot_data_ppc %>%
      filter(time == 0) %>%
      distinct(trial, draw),
    by = c("trial", "draw")
  ) %>%
  mutate(time = 0, cum_prop = 0, type = "ppc")

plot_data_ppc <- bind_rows(plot_data_ppc, starting_points_ppc) %>%
  arrange(trial, draw, time)

# plot it
ggplot() +
  geom_line(data = plot_data_ppc, 
            aes(x = time, y = cum_prop, 
                group = interaction(draw, trial)), alpha = .1) +
  geom_line(data = plot_data_obs, aes(x = time, y = cum_prop), linewidth = 1) +
  labs(x = "Time", y = "Cumulative proportion informed", color = "Trial") +
  theme_minimal()

#Estimated vs observed acquisition times
acqdata = extract_acqTime(full_fit_weights, data_list)

ggplot(acqdata, aes(x = observed_time, y = median_time)) +
  geom_segment(
    aes(x = observed_time, xend = observed_time, 
        y = median_time, yend = observed_time),
    color = "red",
    alpha = 0.3) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(x = "Observed time", y = "Estimated time") +
  theme_minimal()

#Model comparison
comparison_df <- as.data.frame(loo_output$comparison)
comparison_df$model <- rownames(comparison_df)

ggplot(comparison_df, aes(x = reorder(model, elpd_diff), y = elpd_diff)) +
  geom_point(size = 3) + #elpd_diff
  geom_errorbar(aes(ymin = elpd_diff - se_diff, 
                    ymax = elpd_diff + se_diff), width = 0.2) + #SE of elpd diff
  coord_flip() +
  labs(x = "Model", y = "ELPD Difference", title = "Model Comparison") +
  theme_minimal()

pareto_df = as.data.frame(loo_output$pareto_diagnostics)
ggplot(pareto_df, aes(x=observation, y=pareto_k, color=model))+
  geom_point() +
  scale_color_viridis_d(begin=0.2, end=0.7)+
  geom_hline(yintercept = 0.7, linetype="dashed", color="orange")+
  geom_hline(yintercept = 1, linetype="dashed", color="red")+
  labs(x="Observation", y="Pareto-k value", title="Pareto-k diagnostics")+
  theme_minimal()


#(7) FIGURES ----
colours <- c("#FFFFFF", "#999999", "#333333")

patch_binary <- ggplot(all_individuals_prepatch, aes(x = degree, y = patch_binary)) +
  geom_point(aes(fill = sex), shape = 21, size = 3) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x="Degree", y = "Patch discovered (binary)") +
  theme_bw(base_size = 18)
patch_binary

patch_number <- ggplot(all_individuals_prepatch, aes(x = degree, y = patch_number)) +
  geom_point(aes(fill = age), shape = 21, size = 3) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(x="Degree", y = "Patch discovered (binary)") +
  theme_bw(base_size = 18)
patch_number

patch_number <- ggplot(all_individuals_prepatch, aes(x = sex, y = patch_number)) + 
  geom_violin(width = 1) +
  labs(x="Sex", y = "Patch discovery") +
  theme_bw(base_size = 18) +
  stat_summary(fun=mean, geom="point", shape=8, size=3) +  geom_jitter(shape=16, position= jitter, size =3) + theme(legend.position="none")
patch_number


#Data list
data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static
)

#Full social model including transmission weights
model_full_weights <- generate_STb_model(data_list, 
                                         gq = T, 
                                         est_acqTime = F, 
                                         data_type = c("continuous"))

#Fit full social model including transmission weights
full_fit_weights <- fit_STb(data_list,
                            model_full,
                            parallel_chains = 4,
                            chains = 4,
                            cores = 4,
                            iter = 4000,
                            refresh=1000
)

STb_save(full_fit_weights, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(full_fit_weights, digits = 3)

data_list <- import_user_STb(
  event_data = event_data,
  networks = edge_list,
  network_type = "undirected",
  t_weights = t_weights_static
)

model_asoc = generate_STb_model(data_list, model_type="asocial")

asocial_fit_weights = fit_STb(data_list,
                              model_asoc,
                              parallel_chains = 4,
                              chains = 4,
                              cores = 4,
                              iter = 4000,
                              refresh=1000)

STb_save(asocial_fit_weights, output_dir = "cmdstan_saves", name="my_first_fit")

STb_summary(asocial_fit_weights, digits = 3)

#Compare social and asocial models
loo_output = STb_compare(full_fit_weights, asocial_fit_weights, method="loo-psis")

print(loo_output$comparison, simplify = FALSE)

