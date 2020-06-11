library(seqinr)
library(data.table)
library(ggplot2)
library(grid)
library(Hmisc)
library(memisc)


#########################################################
#Load 
#########################################################
#see ReadMe for sample data formats

# NHSS data
# variables: "newnum", "rsh_county_name", "hivyr", "hivmo", "self_rpt_first_pos_yr", "self_rpt_first_pos_mo",
#            "birthRegion", "gender", "age_diag", "transmission", "first_cd4_rnd", "first_vl_rnd" 
nhss <- read.csv("nhss.csv")

# partner services data
# variables: "newnum", "interviewed", "entryUS_date", "HIV_1stpos_date", "HIV_lastneg_date", 
#            "TestBeforeUS", "TestBeforeUSResult", ???
ps <- read.csv("ps.csv")

# MHS sequence-level metadata (note: multiple sequences per person)
# variables: "seqID", "newnum", "type", "KC", "subtype"
metadata <- read.csv("metadata.csv")

# alignment and person-level cluster data from clusters.R
# variables: "newnum", "clusterID", "is.clustered.tn02"
alignment <- read.fasta("alignment.fasta")
clusters <- read.csv("TN93.02.csv") #



#########################################################
#Combine 
#########################################################
#combine data by person-level ID, "newnum"
nhss <- merge(nhss, ps, by="newnum", all.x=TRUE)
nhss <- merge(nhss, clusters, by="newnum", all.x=TRUE)

#aggregate from sequence-level to person-level
metadata$n <- 1
metadata <- with(metadata, aggregate(n, FUN=sum, by=list(newnum, subtype))) 
names(metadata) <- c("newnum", "subtype", "n")
#tag and drop duplicates,keep subtypes that occurs most frequently
metadata <- metadata[order(metadata$newnum, -metadata$n),]
metadata$dup <- duplicated(metadata$newnum) 
metadata <- subset(metadata, dup==FALSE)
#merge person-level subtype data
nhss <- merge(nhss, metadata[,c("newnum", "subtype")], by="newnum", all.x=TRUE)


#########################################################
#New Variables-Dates
#########################################################
#this is complex due to incomplete dates
nhss$hiv_dx_date_complete <- with(nhss, ifelse(hivmo=="" & hivyr=="", "no date",
                                        ifelse(hivmo=="", "no month",
                                        ifelse(hivyr=="", "no year", "complete"))))
#tag differences in hiv date and prev. diag. date
nhss$self_report_different <- with(nhss, ifelse(self_rpt_first_pos_yr=="", "missing",
                                         ifelse(self_rpt_first_pos_yr==hivyr & self_rpt_first_pos_mo=="", "same year/no self reported month", 
                                         ifelse(self_rpt_first_pos_yr==hivyr & self_rpt_first_pos_mo==hivmo, "same date", 
                                         ifelse(self_rpt_first_pos_yr==hivyr & self_rpt_first_pos_mo!=hivmo, "same year/different month", 
                                         ifelse(self_rpt_first_pos_yr!=hivyr, "different year", NA))))))
#create a new variables: *dx_year*, *dx_month*, and *time_prev_to_current_dx*
#case1: keep hivyr and hivmo
nhss$dx_year <- with(nhss, ifelse(self_report_different=="missing" | self_report_different=="same year/no self reported month" | self_report_different=="same date", hivyr, NA))
nhss$dx_month <- with(nhss, ifelse(self_report_different=="missing" | self_report_different=="same year/no self reported month" | self_report_different=="same date", hivmo, NA))
nhss$time_prev_to_current_dx <- with(nhss, ifelse(self_report_different=="missing" | self_report_different=="same year/no self reported month" | self_report_different=="same date", 0, NA))
#case2: keep hivyr, calculate minimum month
nhss$dx_year  <- with(nhss, ifelse(self_report_different=="same year/different month", hivyr, dx_year))
nhss$dx_month <- with(nhss, ifelse(self_report_different=="same year/different month", min(hivmo, self_rpt_first_pos_mo), dx_month))
nhss$time_prev_to_current_dx <- with(nhss, ifelse(self_report_different=="same year/different month", abs(hivmo-self_rpt_first_pos_mo), NA))
#case3: calculate earliest year 
nhss$dx_year <- with(nhss, ifelse(self_report_different=="different year", ifelse(hivyr<self_rpt_first_pos_yr, hivyr, self_rpt_first_pos_yr), dx_year))
nhss$dx_month <- with(nhss, ifelse(self_report_different=="different year", ifelse(hivyr<self_rpt_first_pos_yr, hivmo, self_rpt_first_pos_mo), dx_month))
nhss$time_prev_to_current_dx <- with(nhss, ifelse(self_report_different=="different year", (abs(hivyr - dx_year) - 1)*12 + hivmo_0 + (12-self_rpt_first_pos_mo_0), time_prev_to_current_dx))
#case4: same date
nhss$time_prev_to_current_dx <- with(nhss, ifelse(self_report_different %in% c("same year/no self reported month", "same date"), 0, time_prev_to_current_dx))


#########################################################
#New Variables-Other
#########################################################
nhss$has.ps <- with(nhss, ifelse(is.na(interviewed), 0, interviewed)) #1=yes,0=no
nhss$has.seq <- (nhss$newnum %in% metadata$newnum) #TRUE/FALSE
nhss$has.PRRTseq <- (nhss$newnum %in% metadata[metadata$seqID %in% names(alignment),]$newnum) #TRUE/FALSE
nhss$is.clustered.tn02 <- with(nhss, ifelse(is.clustered.tn02==1, 1, ifelse(has.PRRTseq==TRUE, 0, NA)))

nhss$is.USA    <- with(nhss, ifelse(birthRegion=="USA", "USA", ifelse(birthRegion=="Unknown", "Unknown", "Foreign")))
nhss$is.IDU    <- with(nhss, ifelse(transmission %in% c('MSM/IDU',"IDU", "TRANS-SM/IDU"), 1, 0)) #TRANS-SM = transgender, sex with men
nhss$is.MSM    <- with(nhss, ifelse(transmission %in% c("MSM", "MSM/IDU"), 1, 0))
nhss$is.HETERO <- with(nhss, ifelse(transmission %in% c("HSX_F", "HSX_M"), 1, 0))
nhss$is.OTHER  <- with(nhss, ifelse(transmission %in% c("HSX_TRANS","PERINAT","TRANS-SM","UNKNOWN"), 1, 0))

nhss$age_diag_cat <- with(nhss, ifelse(age_diag <= 24, "<25",
                                ifelse(age_diag >= 25  & age_diag <= 34, "25-34",
                                ifelse(age_diag >= 35  & age_diag <= 44, "35-44",
                                ifelse(age_diag >= 45, "45+", NA)))))
nhss$cd4_cat <- with(nhss, ifelse(first_cd4_rnd<200, "<200", 
                           ifelse(first_cd4_rnd>=200 & first_cd4_rnd<500, "200-500", 
                           ifelse(first_cd4_rnd>=500, ">500", NA))))
nhss$vl_cat <- with(nhss, ifelse(first_vl_rnd<200, "<200", 
                          ifelse(first_vl_rnd>=200 & first_vl_rnd<1000, "200-1,000", 
                          ifelse(first_vl_rnd>=1000 & first_vl_rnd<10000, "1,000-10,000",
                          ifelse(first_vl_rnd>=10000 & first_vl_rnd<100000, "10,000-100,000",
                          ifelse(first_vl_rnd>=100000, ">100,000", NA))))))

nhss[is.na(nhss$subtype),]$subtype <- "No Sequence Available"
nhss$is.B <- with(nhss, ifelse(subtype=="No Sequence Available", "No Sequence Available",
                        ifelse(subtype=="B", "B", "notB")))
nhss$subtype_cat <- with(nhss, ifelse(subtype=="No Sequence Available", "No Sequence Available",
                               ifelse(subtype=="B", "B", 
                               ifelse(subtype=="A1", "A1",
                               ifelse(subtype=="C", "C",
                               ifelse(subtype=="01_AE", "01_AE",
                               ifelse(subtype=="02_AG", "02_AG", "Other")))))))

#########################################################
# Limit Dataset 
#########################################################

nhss <- subset(nhss, rsh_county_name=="KING CO.") #diagnosed with HIV in KC
nhss <- subset(nhss, hivyr>=2010 & hivyr<=2018) #diagnosed 2010-2018
saveRDS(nhss, file = "processeddata.rds")

