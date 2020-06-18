library(ape)
library(data.table)
library(ggplot2)
library(grid)
library(Hmisc)
library(memisc)
library(sp)
library(stringr)
library(viridis)
library(seqinr)
library(tableone)
library(phytools)
library(diversitree)

nhss <- readRDS("processeddata.rds")

#########################################################
# NHSS & PHSKC cases
#########################################################

#identify individuals with a previous diagnosis & those categorized as a "new" or "local" case by PHSKC
#note: nhss$prev_diag=="previously diagnosed" are excluded from primary analysis
nhss$prev_diag <- with(nhss, ifelse( is.na(time_prev_to_current_dx) | time_prev_to_current_dx==0, 
                                     "new PHSKC case", "previously diagnosed"))
nhss.new <- subset(nhss, prev_diag=="new PHSKC case")


#########################################################
# Partner Services data
#########################################################
#see Figure 1 in publication and supplement for PS questions asked

#limit to people with an interview
nhss.new.has.ps <- subset(nhss.new, has.ps)

#ever tested for HIV?
nhss.new.has.ps$nevertested <- with(nhss.new.has.ps, ifelse(HIV_lastneg_date=="" | is.na(HIV_lastneg_date), 1, 0)) #if no negative date available
nhss.new.has.ps$nevertested <- with(nhss.new.has.ps, ifelse(TestBeforeUS==1, 0, nevertested)) #replace with FALSE if they had tested pos/neg before US

#date of arrival to US known?
nhss.new.has.ps$nodate_US <- with(nhss.new.has.ps, ifelse(has.ps==1, ifelse(is.na(entryUS_date), 1, 0), NA))

#date of last negative test occured AFTER arrival to the US?
#note: TestBeforeUSResult results are 0=not tested, 1=positive, 2=negative
nhss.new.has.ps$negative_after_US <- with(nhss.new.has.ps, ifelse(is.na(entryUS_date) | is.na(HIV_lastneg_date) | TestBeforeUSResult==1, NA,
                                                           ifelse(year(entryUS_date)<year(HIV_lastneg_date), 1, 0)))
nhss.new.has.ps$no_previous_negative <- with(nhss.new.has.ps, ifelse(is.na(HIV_lastneg_date), 1, 0))

#reported sex while traveling abroad? 
nhss.new.has.ps$homecountry_sex <- with(nhss.new.has.ps, ifelse(HomeCountryReturnSex==1, 1, NA))

#tested positive for HIV after having sex abroad?
nhss.new.has.ps$positive_after_homecountry_sex <- with(nhss.new.has.ps, ifelse(HomeCountryReturnSex==1 & NegHomeCountrySex==1, 0,
                                                                        ifelse(HomeCountryReturnSex==1 & NegHomeCountrySex==0, 1, NA)))

#inference
nhss.new.has.ps$ps.acquisition <- with(nhss.new.has.ps, ifelse(nevertested==1, "Unsure - Never Tested",
                                                        ifelse(negative_after_US==1 & is.na(homecountry_sex), "US", 
                                                        ifelse(positive_after_homecountry_sex==0, "US", 
                                                        ifelse(positive_after_homecountry_sex==1, "Abroad", 
                                                        ifelse(no_previous_negative==0, "Unsure - Noninformative Testing History", 
                                                        ifelse(no_previous_negative==1, "Unsure - Never Tested", NA)))))))

#merge onto main dataset
nhss.new <- merge(nhss.new, nhss.new.has.ps[,c("newnum","ps.acquisition")], by="newnum", all.x=TRUE)


#########################################################
# Molecular Epi
#########################################################
nhss.new$molecularepi.acquisition <- with(nhss.new, ifelse(is.clustered.tn02==1 & is.B=="B", "Subtype B, Clustered", 
                                                    ifelse(is.clustered.tn02==1 & is.B=="notB", "Subtype non-B, Clustered", 
                                                    ifelse(is.clustered.tn02==0 & is.B=="B", "Subtype B, Not Clustered", 
                                                    ifelse(is.clustered.tn02==0 & is.B=="notB", "Subtype non-B, Not Clustered", NA)))))

#########################################################
# Combined Analysis
#########################################################
nhss.new$combined.acquisition <- with(nhss.new, 
                                      ifelse(ps.acquisition %in% c("Abroad", "US"), ps.acquisition, 
                                      ifelse(molecularepi.acquisition=="Subtype non-B, Not Clustered", "Abroad",
                                      ifelse(molecularepi.acquisition %in% c("Subtype B, Clustered", "Subtype non-B, Clustered"), "US",
                                      ifelse(molecularepi.acquisition=="Subtype B, Not Clustered" & birthRegion=="Sub-Saharan Africa", "US",
                                      ifelse(molecularepi.acquisition=="Subtype B, Not Clustered" & birthRegion %in% c("Asia", "Europe/Canada", "Latin/Southern America", "Middle East/North Africa", "Oceania"), "Non-informative",
                                      "No Data"))))))
nhss.new$combined.acquisition <- with(nhss.new, ifelse(birthRegion=="USA", NA, combined.acquisition))


nhss <- merge(nhss, nhss.new[,c("newnum", "ps.acquisition", "molecularepi.acquisition", "combined.acquisition")], by="newnum", all.x=TRUE)
nhss$combined.acquisition <- with(nhss, ifelse(prev_diag=="previously diagnosed", "previously diagnosed (not a new PHSKC case)", combined.acquisition))
saveRDS(nhss, file = "localHIV.rds")



