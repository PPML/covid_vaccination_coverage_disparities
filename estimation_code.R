#########################################################################################
## Joint KFF-Stanford Project to Estimate COVID-19 Vaccination Coverage by Race/Ethnicity
## Contact: mreitsma@stanford.edu
#########################################################################################

rm(list = ls())

## Load packages
library(data.table)
library(ggplot2)
library(DBI)
library(RSQLite)
library(lubridate)

setwd("~/COVID/vaccination_analysis/KFF Analysis/Final/")

## State FIPS codes for merging
fips <- as.data.table(tidycensus::fips_codes)
fips <- unique(fips[,.(state_name, state_code)])
fips <- fips[, state_code:=as.numeric(state_code)]

#### Prepare population data from 2019 Single-Year ACS
df <- fread("./data/usa_00026.csv")
df <- df[,.(SERIAL, STATEFIP, GQ, PERNUM, PERWT, AGE, RACE, RACED, HISPAN, HISPAND, PUMA)]
df <- df[GQ%in%c(1, 2, 5)] ## Remove institutionalized populations
df <- df[AGE>=12] ## Restrict to eligible population

## Assign race/ethnicity
df <- df[, ethnicity:=ifelse(HISPAN%in%c(1,2,3,4), "Hispanic", "Not Hispanic")]
df <- df[RACE==1, race:="White"]
df <- df[RACE==2, race:="Black"]
df <- df[RACE==3, race:="American Indian or Alaska Native"]
df <- df[RACE%in%c(4,5), race:="Asian"]
df <- df[RACE%in%c(8, 9), race:="Other"]
df <- df[RACE%in%c(7), race:="NEC"]
df <- df[RACED%in%c(630, 680, 681, 682, 685, 689, 690, 699), race:="Native Hawaiian or Other Pacific Islander"]
df <- df[is.na(race), race:="Asian"]

## Account for race == "NEC" based on census approach
race_breakdowns <- copy(df)
race_breakdowns <- race_breakdowns[race!="NEC"]
race_breakdowns <- race_breakdowns[, tot_race_eth:=sum(PERWT, na.rm=T), by = c("ethnicity", "race", "STATEFIP")]
race_breakdowns <- race_breakdowns[, tot_race:=sum(PERWT, na.rm=T), by = c("ethnicity", "STATEFIP")]
race_breakdowns <- race_breakdowns[, pct_race:=tot_race_eth/tot_race]
race_breakdowns <- unique(race_breakdowns[,.(STATEFIP, ethnicity, race, pct_race)])

nec_expand <- df[race=="NEC"]
nec_expand <- nec_expand[, race:=NULL]
nec_expand <- merge(nec_expand, race_breakdowns, by = c("STATEFIP", "ethnicity"), allow.cartesian=T)
nec_expand <- nec_expand[, PERWT:=PERWT*pct_race]

df <- df[race!="NEC"]
df <- rbind(df, nec_expand[, pct_race:=NULL])

df <- df[ethnicity=="Hispanic", race_grp:="Hispanic"]
df <- df[is.na(race_grp) & race == "White", race_grp:="White"]
df <- df[is.na(race_grp) & race == "Black", race_grp:="Black"]
df <- df[is.na(race_grp) & race == "American Indian or Alaska Native", race_grp:="American Indian or Alaska Native"]
df <- df[is.na(race_grp) & race == "Asian", race_grp:="Asian"]
df <- df[is.na(race_grp) & race == "Other", race_grp:="Other"]
df <- df[is.na(race_grp) & race == "Native Hawaiian or Other Pacific Islander", race_grp:="Native Hawaiian or Other Pacific Islander"]

## Hispanic population breakdown by race
df <- merge(df, fips, by.y = "state_code", by.x="STATEFIP")
hisp_breakdown <- copy(df)
hisp_breakdown <- hisp_breakdown[, tot_race_eth:=sum(PERWT, na.rm=T), by = c("race", "ethnicity", "state_name")]
hisp_breakdown <- hisp_breakdown[, tot_eth:=sum(PERWT, na.rm=T), by = c("ethnicity", "state_name")]
hisp_breakdown <- hisp_breakdown[, tot_race:=sum(PERWT, na.rm=T), by = c("race", "state_name")]
hisp_breakdown <- hisp_breakdown[, pct_race:=tot_race_eth/tot_eth]
hisp_breakdown <- hisp_breakdown[, pct_eth:=tot_race_eth/tot_race]

race_eth <- copy(hisp_breakdown)
race_eth <- unique(race_eth[,.(race, ethnicity, state_name, tot_eth, tot_race)])
race_eth <- race_eth[, pct_eth:=tot_eth/sum(tot_eth, na.rm=T), by = c("state_name", "race")]
race_eth <- race_eth[, pct_race:=tot_race/sum(tot_race, na.rm=T), by = c("state_name", "ethnicity")]
race_eth <- race_eth[ethnicity == "Hispanic"]
setnames(race_eth, "race", "race_grp")
race_eth <- race_eth[, .(state_name, race_grp, ethnicity, pct_eth, pct_race)]
race_eth <- melt(race_eth, id.vars = c("race_grp", "state_name", "ethnicity"), value.var = c("pct_eth", "pct_race"))
race_eth <- race_eth[variable=="pct_eth", race_grp:="Hispanic"]
race_eth <- race_eth[, c("variable", "ethnicity"):=NULL]
setnames(race_eth, "value", "pct_pop")
race_eth <- unique(race_eth)

race_breakdown_by_ethnicity <- copy(hisp_breakdown)
race_breakdown_by_ethnicity <- unique(race_breakdown_by_ethnicity[, .(tot_race_eth, state_name, race, ethnicity)])
race_breakdown_by_ethnicity <- race_breakdown_by_ethnicity[, pct_race_eth:=tot_race_eth/sum(tot_race_eth, na.rm=T), by = c("state_name")]

hisp_breakdown <- hisp_breakdown[ethnicity=="Hispanic"]
hisp_breakdown <- unique(hisp_breakdown[,.(race, state_name, pct_race)])
setnames(hisp_breakdown, "race", "race_grp")

## Merge state name
df <- df[AGE>=12, pop_12_allrace:=sum(PERWT, na.rm=T), by = c("state_name")]
df <- df[AGE>=12, pop_12:=sum(PERWT, na.rm=T), by = c("race_grp", "state_name")]
df <- df[!is.na(pop_12)]
df <- unique(df[,.(pop_12, pop_12_allrace, state_name, race_grp)])

#### Read in KFF extractions
kff_full <- NULL
for (i in c("May 3, 2021", "April 26, 2021", "April 19, 2021", "April 12, 2021", "April 5, 2021", "March 29, 2021", "March 15, 2021")) {
  kff <- as.data.table(readxl::read_xlsx("~/COVID/vaccination_analysis/KFF Analysis/COVID-19 Vaccinations by Race_Ethnicity March 1 - May 17.xlsx", sheet = paste0("as of ", i)))
  setnames(kff, "...1", "state_name")
  kff <- kff[!is.na(state_name)]
  kff <- kff[, c("Total People Vaccinated", "Total People Vaccinated with known Race data", "Total People Vaccinated with known ethnicity"):=NULL]
  num_cols <- colnames(kff)[colnames(kff)%like%"%"]
  kff <- kff[, c(num_cols):=lapply(.SD, as.numeric), .SDcols = c(num_cols), with = F]
  kff <- kff[, c("% of Vaccinations with Known Race", "% of Vaccinations with Known Ethnicity"):=NULL]
  
  kff <- melt(kff, id.vars = c("state_name", "Race Categories Include Hispanic Individuals"))
  kff <- kff[, race_grp:=tstrsplit(x = variable, split = " % of Vaccinations")]
  kff <- kff[, variable:=NULL]
  kff <- kff[race_grp=="Hispanic" & is.na(value), flip_toggle:=1]
  kff <- kff[, flip_toggle:=mean(flip_toggle, na.rm=T), by = "state_name"]
  kff <- kff[flip_toggle==1, `Race Categories Include Hispanic Individuals`:=NA]
  kff <- kff[, flip_toggle:=NULL]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", value_adj:=value]
  
  kff <- kff[race_grp!="% of Vaccinations with Unknown Race" & race_grp!="% of Vaccinations with Unknown Ethnicity"]
  kff <- kff[race_grp=="Other", value:=NA] # Other unrealistically high, consider these data unknown for computing uptake rates
  
  kff <- merge(kff, race_eth, by = c("race_grp", "state_name"), all.x=T)
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic", value:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", fixed_hisp:=value_adj]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", fixed_hisp:=mean(fixed_hisp, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic" & !is.na(value), pct_pop:=pct_pop/sum(pct_pop, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_race:=value/pct_pop]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", rr_eth:=rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_hisp:=mean(rr_eth, na.rm=T), by = "state_name"]
  
  kff <- merge(kff, race_breakdown_by_ethnicity[ethnicity=="Hispanic"], by.x = c("state_name", "race_grp"), by.y=c("state_name", "race"), all.x=T)
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint:=pct_race_eth*rr_hisp*rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint_hisp_sum:=sum(temp_joint, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", final_joint:=temp_joint/(temp_joint_hisp_sum/fixed_hisp)]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & !is.na(final_joint), value_adj:=value-final_joint]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="" | is.na(`Race Categories Include Hispanic Individuals`), value_adj:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[, Date:=mdy(paste0(i))]
  kff_full <- rbind(kff_full, kff, fill = T)
}

for (i in c("05-10-2021", "05-17-2021", "05-24-2021", "06-07-2021")) {
  kff <- as.data.table(readxl::read_xlsx("~/COVID/vaccination_analysis/KFF Analysis/Vaccine Distribution 5_10-6_7.xlsx", sheet = paste0(i)))
  setnames(kff, "...1", "state_name")
  kff <- kff[!is.na(state_name)]
  kff <- kff[, c("Total People Vaccinated", "Total People Vaccinated with known Race data", "Total People Vaccinated with known ethnicity"):=NULL]
  num_cols <- colnames(kff)[colnames(kff)%like%"%"]
  kff <- kff[, c(num_cols):=lapply(.SD, as.numeric), .SDcols = c(num_cols), with = F]
  kff <- kff[, c("% of Vaccinations with Known Race", "% of Vaccinations with Known Ethnicity"):=NULL]
  
  kff <- melt(kff, id.vars = c("state_name", "Race Categories Include Hispanic Individuals"))
  kff <- kff[, race_grp:=tstrsplit(x = variable, split = " % of Vaccinations")]
  kff <- kff[, variable:=NULL]
  kff <- kff[race_grp=="Hispanic" & is.na(value), flip_toggle:=1]
  kff <- kff[, flip_toggle:=mean(flip_toggle, na.rm=T), by = "state_name"]
  kff <- kff[flip_toggle==1, `Race Categories Include Hispanic Individuals`:=NA]
  kff <- kff[, flip_toggle:=NULL]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", value_adj:=value]
  
  kff <- kff[race_grp!="% of Vaccinations with Unknown Race" & race_grp!="% of Vaccinations with Unknown Ethnicity"]
  kff <- kff[race_grp=="Other", value:=NA] # Other unrealistically high, consider these data unknown for computing uptake rates
  
  kff <- merge(kff, race_eth, by = c("race_grp", "state_name"), all.x=T)
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic", value:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", fixed_hisp:=value_adj]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", fixed_hisp:=mean(fixed_hisp, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic" & !is.na(value), pct_pop:=pct_pop/sum(pct_pop, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_race:=value/pct_pop]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", rr_eth:=rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_hisp:=mean(rr_eth, na.rm=T), by = "state_name"]
  
  kff <- merge(kff, race_breakdown_by_ethnicity[ethnicity=="Hispanic"], by.x = c("state_name", "race_grp"), by.y=c("state_name", "race"), all.x=T)
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint:=pct_race_eth*rr_hisp*rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint_hisp_sum:=sum(temp_joint, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", final_joint:=temp_joint/(temp_joint_hisp_sum/fixed_hisp)]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & !is.na(final_joint), value_adj:=value-final_joint]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="" | is.na(`Race Categories Include Hispanic Individuals`), value_adj:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[, Date:=mdy(paste0(i))]
  kff_full <- rbind(kff_full, kff, fill = T)
}

kff_full <- merge(kff_full, df, by = c("race_grp", "state_name"), all=T)
kff_full <- kff_full[!is.na(value_adj), pop_pct:=pop_12/sum(pop_12, na.rm=T), by = c("state_name", "Date")]

#### Read in and clean CDC state vaccination time series
vax_stats <- fread("~/Downloads/COVID-19_Vaccinations_in_the_United_States_Jurisdiction (1).csv") #https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc
vax_stats <- vax_stats[, Date:=mdy(Date)]
setnames(vax_stats, "Location", "state")
fips <- as.data.table(tidycensus::fips_codes)
fips <- unique(fips[,.(state_name, state, state_code)])
vax_stats <- merge(vax_stats, fips, by = c("state"))

# Impute 12+ for full time series
vax_stats <- vax_stats[, pop12:=Administered_Dose1_Recip_12Plus/(Administered_Dose1_Recip_12PlusPop_Pct/100)]
vax_stats <- vax_stats[, pop12:=mean(pop12, na.rm=T), by = "state_name"]
vax_stats <- vax_stats[Administered_Dose1_Recip_12Plus==0, Administered_Dose1_Recip_12Plus:=Administered_Dose1_Recip]
vax_stats <- vax_stats[Administered_Dose1_Recip_12PlusPop_Pct==0, Administered_Dose1_Recip_12PlusPop_Pct:=(Administered_Dose1_Recip_12Plus/pop12)*100]

vax_stats <- vax_stats[,.(Date, state_name, Administered_Dose1_Recip_12Plus, Administered_Dose1_Recip_12PlusPop_Pct)]
vax_stats <- vax_stats[, lapply(.SD, as.numeric), by = c("state_name", "Date")]
# vax_stats <- vax_stats[state_name=="New York State", state_name:="New York"]

vax_stats <- merge(vax_stats, kff_full, by = c("state_name", "Date"))
vax_stats <- vax_stats[, pop_share_agg:=pop_12/sum(pop_12, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[is.na(value_adj), missing:=1]
vax_stats <- vax_stats[!is.na(value_adj), missing:=0]
vax_stats <- vax_stats[is.na(value_adj), pct_missing:=pop_share_agg]
vax_stats <- vax_stats[, pct_missing:=sum(pct_missing, na.rm=T), by = c("state_name", "Date")]

vax_stats <- vax_stats[, value_adj:=value_adj*(1-pct_missing)]
vax_stats <- vax_stats[is.na(value_adj), value_adj:=pop_share_agg]
vax_stats <- vax_stats[, value_adj:=value_adj/sum(value_adj, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[, doses_adj:=(Administered_Dose1_Recip_12PlusPop_Pct/100)*pop_12_allrace]
vax_stats <- vax_stats[, race_doses:=doses_adj*value_adj]

## Deal with 100+ round 1
vax_stats <- vax_stats[race_doses>pop_12, add_doses:=race_doses-pop_12]
vax_stats <- vax_stats[, add_doses:=sum(add_doses, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[race_doses>pop_12, exceeded:=1]
vax_stats <- vax_stats[race_doses>pop_12, race_doses:=pop_12]
# vax_stats <- vax_stats[is.na(exceeded), second_pct:=pop_share_agg/sum(pop_share_agg), by = c("state_name", "Date")]
vax_stats <- vax_stats[is.na(exceeded), second_pct:=value_adj/sum(value_adj), by = c("state_name", "Date")]
vax_stats <- vax_stats[is.na(exceeded), race_doses:=race_doses+(second_pct*add_doses)]
# 
# ## Deal with 100+ round 2
# vax_stats <- vax_stats[, c("add_doses"):=NULL]
# vax_stats <- vax_stats[race_doses>pop_12, add_doses:=race_doses-pop_12]
# vax_stats <- vax_stats[, add_doses:=sum(add_doses, na.rm=T), by = c("state_name", "Date")]
# vax_stats <- vax_stats[race_doses>pop_12, exceeded:=1]
# vax_stats <- vax_stats[race_doses>pop_12, race_doses:=pop_12]
# vax_stats <- vax_stats[is.na(exceeded), second_pct:=value_adj/sum(value_adj), by = c("state_name", "Date")]
# vax_stats <- vax_stats[is.na(exceeded), race_doses:=race_doses+(second_pct*add_doses)]

vax_stats <- vax_stats[, vaccinated_pct_12:=(race_doses/pop_12)*100]

pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/historical_coverage_6.7.21.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & missing!=1], aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & missing!=1], aes(y = (doses_adj/pop_12_allrace)*100), size = 1.2, color = "dark gray") +
    geom_point(size = 2) + geom_line( size = 1.3) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"), values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage", color = "Race/Ethnicity", title = paste0(i)) +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 70, linetype = "dashed")+
    scale_y_continuous(limits = c(0, 100))
  print(plot)
}
dev.off()

## PER UNVACCINATED ##
backup <- copy(vax_stats)

## NEBRASKA MISSING 6/7
vax_stats <- copy(backup)
vax_stats <- vax_stats[state_name%in%"Nebraska"]
vax_stats <- vax_stats[order(state_name, race_grp, Date)]
vax_stats <- vax_stats[, l2:=shift(race_doses, n = 2, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff:=race_doses-l2]
vax_stats <- vax_stats[, diff_rate:=(diff/(pop_12*14))]

vax_stats <- vax_stats[, all_race:=doses_adj]
vax_stats <- vax_stats[, l2_all:=shift(all_race, n = 2, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff_all:=all_race-l2_all]
vax_stats <- vax_stats[, diff_all_rate:=(diff_all/(pop_12_allrace*14))]
vax_stats_out <- vax_stats

## Idaho/Tennessee missing 5/24
vax_stats <- copy(backup)
# vax_stats <- vax_stats[!(state_name%in%c("Nebraska"))]
vax_stats <- vax_stats[state_name%in%c("Idaho", "Tennessee", "New Hampshire")]
vax_stats <- vax_stats[order(state_name, race_grp, Date)]
vax_stats <- vax_stats[, l2:=shift(race_doses, n = 2, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff:=race_doses-l2]
vax_stats <- vax_stats[, diff_rate:=(diff/(pop_12*21))]

vax_stats <- vax_stats[, all_race:=doses_adj]
vax_stats <- vax_stats[, l2_all:=shift(all_race, n = 2, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff_all:=all_race-l2_all]
vax_stats <- vax_stats[, diff_all_rate:=(diff_all/(pop_12_allrace*21))]
vax_stats_out <- rbind(vax_stats_out, vax_stats)

## Remainder use difference between 5/24 and 6/7
vax_stats <- copy(backup)
vax_stats <- vax_stats[!(state_name%in%c("Idaho", "Tennessee", "Nebraska", "New Hampshire"))]
vax_stats <- vax_stats[order(state_name, race_grp, Date)]
vax_stats <- vax_stats[, l2:=shift(race_doses, n = 1, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff:=race_doses-l2]
vax_stats <- vax_stats[, diff_rate:=(diff/(pop_12*14))]

vax_stats <- vax_stats[, all_race:=doses_adj]
vax_stats <- vax_stats[, l2_all:=shift(all_race, n = 1, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff_all:=all_race-l2_all]
vax_stats <- vax_stats[, diff_all_rate:=(diff_all/(pop_12_allrace*14))]
vax_stats_out <- rbind(vax_stats_out, vax_stats)
vax_stats <- copy(vax_stats_out)

## Edge cases
vax_stats <- vax_stats[diff_rate<0, diff_rate:=0]
vax_stats <- vax_stats[diff_all_rate<0, diff_all_rate:=0]
vax_stats <- vax_stats[!(state_name%in%c("Nebraska")) |
                         (state_name == "Nebraska" & Date<="2021-05-24")]

backup <- copy(vax_stats)
vax_stats <- copy(backup)

for (i in ymd("2021-06-07"):ymd("2021-08-31")) {
  projection <- vax_stats[Date==i & !(state_name%in%c("Nebraska"))]
  # Per population
  projection <- projection[, race_doses:=race_doses+(diff_rate*(pop_12))]
  projection <- projection[, doses_adj:=doses_adj+(diff_all_rate*(pop_12_allrace))]
  
  projection <- projection[, Date:=Date+1]
  projection <- projection[, project:=1]
  vax_stats <- rbind(vax_stats, projection, fill = T)
}

for (i in ymd("2021-05-24"):ymd("2021-08-31")) {
  projection <- vax_stats[Date==i & (state_name%in%c("Nebraska"))]
  # Per population
  projection <- projection[, race_doses:=race_doses+(diff_rate*(pop_12))]
  projection <- projection[, doses_adj:=doses_adj+(diff_all_rate*(pop_12_allrace))]
  
  projection <- projection[, Date:=Date+1]
  projection <- projection[, project:=1]
  vax_stats <- rbind(vax_stats, projection, fill = T)
}

vax_stats <- vax_stats[, vaccinated_pct_12:=(race_doses/pop_12)*100]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), race_doses:=pop_12]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), vaccinated_pct_12:=100]

vax_stats <- vax_stats[, tot_doses:=sum(race_doses, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[, agg:=(tot_doses/pop_12_allrace)*100]

vax_stats_out <- copy(vax_stats)
vax_stats_out <- vax_stats_out[Date=="2021-07-04"]
vax_stats_out <- vax_stats_out[, upper:=vaccinated_pct_12+1.96*sqrt((vaccinated_pct_12*(100-vaccinated_pct_12))/pop_12)]
vax_stats_out <- vax_stats_out[, moe:=upper-vaccinated_pct_12]
vax_stats_out <- vax_stats_out[, vaccinated_pct_12_out:=as.character(round(vaccinated_pct_12))]
vax_stats_out <- vax_stats_out[, All:=as.character(round(agg))]
vax_stats_out <- vax_stats_out[missing == 1, vaccinated_pct_12_out:="NR"]
vax_stats_out <- vax_stats_out[missing == 0 & (moe>0.5 | race_doses<1500), vaccinated_pct_12_out:="NE"]

suppress_table <- copy(vax_stats_out)
suppress_table <- suppress_table[, nr:=ifelse(missing == 1, 1, 0)]
suppress_table <- suppress_table[, ne:=ifelse(missing == 0 & (moe>0.5 | race_doses<1500), 1, 0)]
suppress_table <- unique(suppress_table[,.(state_name, race_grp, nr, ne)])
suppress_table <- suppress_table[race_grp%in%c("Asian", "Black", "Hispanic", "White")]

vax_stats_out <- vax_stats_out[vaccinated_pct_12>=70, vaccinated_pct_12_out:="> 70%"]
vax_stats_out <- vax_stats_out[agg>=70, All:="> 70%"]


vax_stats_out <- vax_stats_out[,.(state_name, race_grp, All, vaccinated_pct_12_out)]
vax_stats_out <- vax_stats_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_stats_out <- dcast(vax_stats_out, state_name + All ~ race_grp, value.var = "vaccinated_pct_12_out")

write.csv(vax_stats_out, "~/COVID/vaccination_analysis/KFF Analysis/Final/output/coverage_july_4.csv", na = "", row.names = F)

vax_date_out <- copy(vax_stats)
vax_date_out <- vax_date_out[,.(state_name, Date, race_grp, vaccinated_pct_12)]
vax_date_out <- vax_date_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_dates_square <- as.data.table(expand.grid(Date = seq(ymd("2021-03-15"), ymd("2021-09-01"), by = "day"), race_grp = c("Asian", "Black", "Hispanic", "White"), state_name = unique(vax_date_out$state_name)))
vax_date_out <- merge(vax_date_out, vax_dates_square, by = c("Date", "race_grp", "state_name"), all = T)
vax_date_out <- vax_date_out[order(state_name, race_grp, Date)]
vax_date_out <- vax_date_out[, vaccinated_pct_12_approx:=na.approx(vaccinated_pct_12), by = c("state_name", "race_grp")]

vax_date_out <- vax_date_out[vaccinated_pct_12_approx>=70, date_:=min(Date, na.rm=T), by = c("state_name", "race_grp")]
vax_date_out <- vax_date_out[, date_:=mean(date_, na.rm=T), by = c("state_name", "race_grp")]
vax_date_out <- vax_date_out[, date_:=as.character(date_)]
vax_date_out <- vax_date_out[is.na(date_), date_:="After September 1"]

vax_date_out <- unique(vax_date_out[,.(state_name, race_grp, date_)])
vax_date_out <- vax_date_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_date_out <- merge(vax_date_out, suppress_table, by = c("state_name", "race_grp"))
vax_date_out <- vax_date_out[nr==1, date_:="NR"]
vax_date_out <- vax_date_out[ne==1, date_:="NE"]
vax_date_out <- vax_date_out[, c("nr", "ne"):=NULL]
vax_date_out <- dcast(vax_date_out, state_name ~ race_grp, value.var = "date_")
vax_date_out <- vax_date_out[,.(state_name, White, Black, Hispanic, Asian)]

write.csv(vax_date_out, "~/COVID/vaccination_analysis/KFF Analysis/Final/output/date_70_state.csv", na = "", row.names = F)

# vax_stats <- merge(vax_stats, suppress_table, by = c("state_name", "race_grp"))

pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/projections_state_6.7.21_twoweek.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & is.na(project) & nr == 0 & ne == 0 & missing==0],
                 aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & is.na(project)], aes(y = agg), size = 1, color = "dark gray") +
    geom_point(size = 2) + geom_line( size = 1.1) +
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & (project == 1 | Date == "2021-06-07") & Date <="2021-07-04"], aes(y = agg), size = 1.2, color = "dark gray",
              linetype = "dashed", alpha = 0.9) +
    geom_line(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & (project == 1 | Date == "2021-06-07") & Date <="2021-07-04"  & nr == 0 & ne == 0],
              size = 1.3,
              linetype = "dashed", alpha = 0.9) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"),
                       values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage Among Eligible (12+)", color = "Race/Ethnicity", title = paste0(i),
         caption = "Gray line shows coverage among whole eligible (12+) population.\nSolid lines show historical (observed) data, dashed lines show projections.") +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 70, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_shape_manual(values = c(19, 1))
  print(plot)
}
dev.off()

vax_stats_nat <- copy(vax_stats)
vax_stats_nat <- vax_stats_nat[, nat_doses_race:=sum(race_doses, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_race:=sum(pop_12, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_doses_agg:=sum(race_doses, na.rm=T), by = c("Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_agg:=sum(pop_12, na.rm=T), by = c("Date")]

vax_stats_nat <- vax_stats_nat[, nat_race:=(nat_doses_race/nat_pop_race)*100]
vax_stats_nat <- vax_stats_nat[, nat_race_agg:=(nat_doses_agg/nat_pop_agg)*100]

vax_stats_nat <- unique(vax_stats_nat[,.(race_grp, Date, nat_race, nat_race_agg, project)])
vax_stats_nat <- vax_stats_nat[(Date>"2021-06-07" & project==1) | is.na(project)]

write.csv(vax_stats_nat[Date=="2021-07-04"], "~/COVID/vaccination_analysis/KFF Analysis/Final/output/national_current.csv", na = "", row.names = F)
write.csv(vax_stats_nat[Date<="2021-07-04"], "~/COVID/vaccination_analysis/KFF Analysis/Final/output/national_time_series.csv", na = "", row.names = F)

vax_stats_out <- copy(vax_stats_nat)
vax_stats_out <- vax_stats_out[Date=="2021-07-04"]
vax_stats_out <- vax_stats_out[,.(race_grp, nat_race, nat_race_agg)]
vax_stats_out <- vax_stats_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_stats_out <- dcast(vax_stats_out, nat_race_agg ~ race_grp, value.var = "nat_race")
setnames(vax_stats_out, "nat_race_agg", "All")
vax_stats_out <- vax_stats_out[, state_name:="United States"]

write.csv(vax_stats_out, "~/COVID/vaccination_analysis/KFF Analysis/Final/output/coverage_july_4_national.csv", na = "", row.names = F)

vax_date_out <- copy(vax_stats_nat)
vax_dates_square <- as.data.table(expand.grid(Date = seq(ymd("2021-03-15"), ymd("2021-09-01"), by = "day"), race_grp = c("Asian", "Black", "Hispanic", "White")))
vax_date_out <- merge(vax_date_out, vax_dates_square, by = c("Date", "race_grp"), all = T)
vax_date_out <- vax_date_out[order(race_grp, Date)]
vax_date_out <- vax_date_out[, nat_race_approx:=na.approx(nat_race), by = c("race_grp")]
vax_date_out <- vax_date_out[, nat_agg_approx:=na.approx(nat_race_agg), by = c("race_grp")]


vax_date_out <- vax_date_out[nat_race_approx>=70, date_:=min(Date, na.rm=T), by = c( "race_grp")]
vax_date_out <- vax_date_out[, date_:=mean(date_, na.rm=T), by = c("race_grp")]
vax_date_out <- vax_date_out[, date_:=as.character(date_)]
vax_date_out <- vax_date_out[is.na(date_), date_:="After September 1"]
vax_date_out <- vax_date_out[nat_agg_approx>=70, date_agg:=min(Date, na.rm=T), by = c("race_grp")]
vax_date_out <- vax_date_out[, date_agg:=mean(date_agg, na.rm=T), by = c("race_grp")]
vax_date_out <- vax_date_out[, date_agg:=as.character(date_agg)]
vax_date_out <- vax_date_out[is.na(date_agg), date_agg:="After September 1"]

vax_date_out <- vax_date_out[, state_name:="United States"]
vax_date_out <- unique(vax_date_out[,.(state_name, race_grp, date_, date_agg)])
vax_date_out <- vax_date_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_date_out <- dcast(vax_date_out, state_name + date_agg ~ race_grp, value.var = "date_")
setnames(vax_date_out, "date_agg", "All")

write.csv(vax_date_out, "~/COVID/vaccination_analysis/KFF Analysis/Final/output/date_70_national.csv", na = "", row.names = F)



pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/national_scale_up_share_6.7.21.pdf", width =8, height = 6)
plot <- ggplot(data = vax_stats_nat[race_grp%in%c("Hispanic", "White", "Black", "Asian") & is.na(project)],
               aes(x = Date, y = nat_race, color = race_grp)) + 
  geom_line(data = vax_stats_nat[race_grp=="White" & is.na(project)], aes(y = nat_race_agg), size = 1, color = "dark gray") +
  geom_point(size = 2) + geom_line( size = 1.1) +
  geom_line(data = vax_stats_nat[race_grp=="White" & (project == 1 | Date == "2021-06-07") & Date>="2021-06-07" & Date <="2021-07-04"], aes(y = nat_race_agg), size = 1.2, color = "dark gray",
            linetype = "dashed", alpha = 0.9) +
  geom_line(data = vax_stats_nat[race_grp%in%c("Hispanic", "White", "Black", "Asian") & 
                                   (project == 1 | Date == "2021-06-07") & Date>="2021-06-07" & Date <="2021-07-04"], size = 1.3,
            linetype = "dashed", alpha = 0.9) +
  theme_bw() +
  scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"), values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
  labs(x = "", y = "Coverage Among Eligible (12+)", color = "Race/Ethnicity", title = "United States",
       caption = "Gray line shows coverage among whole eligible (12+) population.\nSolid lines show historical (observed) data, dashed lines show projections.") +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 70, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 90)) +
  scale_shape_manual(values = c(19, 1))
print(plot)

plot <- ggplot(data = vax_stats_nat[race_grp%in%c("Hispanic", "White", "Black", "Asian") & is.na(project)],
               aes(x = Date, y = nat_race, color = race_grp)) + 
  geom_line(data = vax_stats_nat[race_grp=="White" & is.na(project)], aes(y = nat_race_agg), size = 1, color = "dark gray") +
  geom_point(size = 2) + geom_line( size = 1.1) +
  geom_line(data = vax_stats_nat[race_grp=="White" & (project == 1 | Date == "2021-06-07") & Date>="2021-06-07" & Date <="2021-08-31"], aes(y = nat_race_agg), size = 1.2, color = "dark gray",
            linetype = "dashed", alpha = 0.9) +
  geom_line(data = vax_stats_nat[race_grp%in%c("Hispanic", "White", "Black", "Asian") & 
                                   (project == 1 | Date == "2021-06-07") & Date>="2021-06-07" & Date <="2021-08-31"], size = 1.3,
            linetype = "dashed", alpha = 0.9) +
  theme_bw() +
  scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"), values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
  labs(x = "", y = "Coverage Among Eligible (12+)", color = "Race/Ethnicity", title = "United States",
       caption = "Gray line shows coverage among whole eligible (12+) population.\nSolid lines show historical (observed) data, dashed lines show projections.") +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 70, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 90)) +
  scale_shape_manual(values = c(19, 1))
print(plot)
dev.off()

vax_stats <- vax_stats[, missing_data:=ifelse(is.na(value), 1, 0)]

plot_data <- copy(vax_stats)
plot_data <- merge(plot_data, suppress_table, by = c("state_name", "race_grp"))

plot_data <- plot_data[Date=="2021-07-04"]
plot_data <- plot_data[, binary:=ifelse(vaccinated_pct_12>=70, 1, 0)]
all <- plot_data[race_grp=="White"]
all <- all[, race_grp:="All"]
plot_data <- rbind(plot_data, all)
plot_data[race_grp=="All", binary:=ifelse(agg>=70, 1, 0)]
plot_data[race_grp=="All", vaccinated_pct_12:=agg]

plot_data <- plot_data[vaccinated_pct_12<50 , cov_cat:="30-49"]
plot_data <- plot_data[vaccinated_pct_12<60 & vaccinated_pct_12>=50, cov_cat:="50-59"]
plot_data <- plot_data[vaccinated_pct_12<70 & vaccinated_pct_12>=60, cov_cat:="60-69"]
plot_data <- plot_data[vaccinated_pct_12>=70, cov_cat:="70+"]
plot_data <- plot_data[(missing_data==1 & race_grp!="All") | ne==1 | nr ==1, cov_cat:="Data or Estimates\nNot Reported"]

plot_data <- plot_data[, cov_cat:=factor(cov_cat, levels = c("30-49", "50-59", "60-69", "70+", "Data or Estimates\nNot Reported"))]

pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/tile_plot.pdf", width = 8, height = 12)
ggplot(data = plot_data[race_grp%in%c("All", "Black", "White", "Asian", "Hispanic")],
       aes(x = race_grp, y = state_name)) +
  geom_tile(aes(fill = cov_cat)) +
  scale_fill_manual(values = c('#ffa500', '#ffbf6b', '#ffd7ac', '#6d87ba', "#eeeeee")) +
  theme_bw() +
  labs(x = "", y = "", fill = "Coverage Among\nEligible (12+)",
       title = "Projected Coverage Among Eligible (12+) on July 4, 2021",
       caption = "Estimates not reported if margin of binomial error exceeded 0.5\npercentage points or a group received fewer than 1500 vaccinations") +
  theme(legend.position = "bottom")+
  theme(legend.key.width=unit(1.5,"cm")) +
  guides(fill=guide_legend(nrow=4,byrow=TRUE)) +
  scale_y_discrete(limits = rev(levels(as.factor(plot_data$state_name))), expand = c(0,0)) +
  theme(text = element_text(size = 14)) +
  scale_x_discrete(expand = c(0,0))
dev.off()  

plot_data <- plot_data[, state_name_ordered:=factor(state_name, levels = unique(plot_data$state_name[order(plot_data$agg)]))]
plot_data <- plot_data[agg>80, agg:=80]
plot_data <- plot_data[vaccinated_pct_12>80, vaccinated_pct_12:=80]

pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/dots_july4.pdf", width = 8, height = 10)

plot_data <- plot_data[pop_12<500000, pop_cat:="0.0-0.4"]
plot_data <- plot_data[pop_12<1000000 & pop_12>=500000, pop_cat:="0.5-0.9"]
plot_data <- plot_data[pop_12>=1000000 & pop_12<2000000, pop_cat:="1.0-1.9"]
plot_data <- plot_data[pop_12>=2000000, pop_cat:="2.0+"]
plot_data <- plot_data[, pop_cat:=factor(pop_cat, levels = c("0.0-0.4", "0.5-0.9", "1.0-1.9", "2.0+"))]

ggplot(data = plot_data[race_grp%in%c("Black", "White", "Asian", "Hispanic") & ne!=1 & nr!=1],
       aes(x = vaccinated_pct_12, y = state_name_ordered)) +
  geom_point(data = plot_data[race_grp=="White"], aes(x = agg), color = "#222222", fill = "#555555", size = 3, shape = 23, alpha = 0) +
  geom_point(data = plot_data[race_grp%in%c("Black", "White", "Asian", "Hispanic") & missing == 0 & ne!=1 & nr!=1],
             aes(color = race_grp, size = pop_cat), alpha = 0.9) +
  geom_point(data = plot_data[race_grp=="White"], aes(x = agg), color = "#222222", fill = "#555555", size = 3, shape = 23, alpha = .8) +
  scale_fill_manual(values = c('#ffa500', '#ffbf6b', '#ffd7ac', '#6d87ba', "#eeeeee")) +
  theme_bw() +
  labs(x = "Coverage Among Eligible (12+)", y = "", color = "Race/Ethnicity", size = "Population Size (Millions)",
       title = "Projected Coverage Among Eligible (12+) on July 4, 2021") +
  theme(legend.position = "bottom")+
  # theme(legend.key.width=unit(1.5,"cm")) +
  guides(color=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size = 4))) +
  guides(size=guide_legend(nrow=4,byrow=TRUE)) +
  # scale_y_discrete(limits = rev(levels(as.factor(plot_data$state_name)))) +
  theme(text = element_text(size = 14)) +
  scale_size_manual(values = c(1, 2, 4, 6)) +
  scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"),
                     values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
  geom_vline(xintercept = 70, linetype = "dashed")+
  scale_x_continuous(labels = c("30%", "40%", "50%", "60%", "70%", "80%+"), breaks = c(30, 40, 50, 60, 70, 80), limits = c(30, 80))
dev.off()

## HEX MAP
plot_data <- copy(vax_stats)

plot_data <- plot_data[Date=="2021-07-04"]
plot_data <- plot_data[, binary:=ifelse(vaccinated_pct_12>=70, 1, 0)]
all <- plot_data[race_grp=="White"]
all <- all[, race_grp:="All"]
plot_data <- rbind(plot_data, all)
plot_data[race_grp=="All", binary:=ifelse(agg>=70, 1, 0)]
plot_data[race_grp=="All", vaccinated_pct_12:=agg]

plot_data <- merge(plot_data, suppress_table, by = c("state_name", "race_grp"))

plot_data <- plot_data[race_grp%in%c("Black", "White", "Hispanic", "Asian")]
plot_data <- plot_data[binary==1, col_cat:="On Track to 70%"]
plot_data <- plot_data[binary==0, col_cat:="Not on Track to 70%"]
plot_data <- plot_data[ne == 1 | nr == 1 | missing == 1, col_cat:="Data or Estimates\nNot Reported"]
plot_data <- plot_data[,.(race_grp, col_cat, state_name, binary)]


library(tidyverse)
library(geojsonio)
library(RColorBrewer)
library(rgdal)
# Load this file. (Note: I stored in a folder called DATA)
spdf <- geojson_read("~/Downloads/us_states_hexgrid.geojson",  what = "sp")

# Bit of reformating
spdf@data = spdf@data %>%
  mutate(google_name = gsub(" \\(United States\\)", "", google_name))

# Show it
plot(spdf)

library(broom)
spdf@data = spdf@data %>% mutate(google_name = gsub(" \\(United States\\)", "", google_name))
spdf_fortified <- tidy(spdf, region = "google_name")

# Calculate the centroid of each hexagon to add the label:
library(rgeos)
centers <- as.data.table(cbind.data.frame(data.frame(gCentroid(spdf, byid=TRUE), id=spdf@data$iso3166_2)))

spdf_fortified <- merge(spdf_fortified, plot_data, by.x = "id", by.y="state_name")
spdf_fortified <- as.data.table(spdf_fortified)

# Now I can plot this shape easily as described before:
pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/hex_maps_binary.pdf", width = 14, height = 10)
plot <- ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group, fill=col_cat), color="white") +
  # geom_text(data=centers, aes(x=x, y=y, label=id), color = "white") +
  geom_text(data=centers, aes(x=x, y=y, label=id), color = "white") +
  theme_void() +
  coord_map() +
  scale_fill_manual(values = c('#C1C1C1', '#ffa500', '#00579f')) +
  # scale_fill_viridis_c(option = "A", na.value = "gray") +
  # scale_fill_viridis_c(option = "A", breaks = c(-17, -15, -10, -5, 0, 5, 10, 12), na.value = "gray") +
  theme(legend.position = "bottom", legend.key.width=unit(1.8,"cm")) +
  guides(fill=guide_legend(nrow=1)) +
  labs(fill = "", title = "States on Track to Reach 70% Coverage Among Eligible (12+) by July 4, by Race/Ethnicity\n") +
  #caption = "*Gray colors do not have sufficient non-white population size") + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(strip.text = element_text(size = 12), text = element_text(size = 14)) +
  facet_wrap(~race_grp)
print(plot)
dev.off()

#### Acceleration to 70%
plot_data <- copy(vax_stats)
plot_data <- plot_data[Date=="2021-06-07"]
plot_data <- plot_data[, accel:=((70-vaccinated_pct_12)/as.numeric((ymd("2021-07-04")-ymd("2021-06-07"))))/100]
accel_df <- unique(plot_data[,.(accel, race_grp, state_name)])

vax_stats <- vax_stats[Date<="2021-06-07"]
vax_stats <- merge(vax_stats, accel_df, by = c("race_grp", "state_name"))

plot_data <- copy(vax_stats)
plot_data <- plot_data[Date=="2021-06-07"]
plot_data <- plot_data[, speed_factor:=accel/diff_rate]
plot_data <- plot_data[, speed_factor_pct:=(accel-diff_rate)/diff_rate]

plot_data <- unique(plot_data[,.(state_name, speed_factor, race_grp, pop_12, missing)])
plot_data <- plot_data[race_grp%in%c("White", "Asian", "Black", "Hispanic")]
plot_data <- plot_data[speed_factor<=1, col_cat:="On track to 70%"]
plot_data <- plot_data[speed_factor>1 & speed_factor<1.5, col_cat:="1-49% faster"]
plot_data <- plot_data[speed_factor>=1.5 & speed_factor<2, col_cat:="50-99% faster"]
plot_data <- plot_data[speed_factor>=2 & speed_factor<2.5, col_cat:="100-149% faster"]
plot_data <- plot_data[speed_factor>=2.5 & speed_factor<3, col_cat:="150-199% faster"]
plot_data <- plot_data[speed_factor>=3, col_cat:="200%+ faster"]
plot_data <- merge(plot_data, suppress_table, by = c("race_grp", "state_name"))
plot_data <- plot_data[nr==1 | ne==1 | missing == 1, col_cat:="Data or Estimates\nNot Reported"]

plot_data <- plot_data[, col_cat:=factor(col_cat, levels = c("On track to 70%", "1-49% faster", 
                                                             "50-99% faster", "100-149% faster", "150-199% faster", "200%+ faster",
                                                             "Data or Estimates\nNot Reported"))]

# Load this file. (Note: I stored in a folder called DATA)
spdf <- geojson_read("~/Downloads/us_states_hexgrid.geojson",  what = "sp")

# Bit of reformating
spdf@data = spdf@data %>%
  mutate(google_name = gsub(" \\(United States\\)", "", google_name))

spdf@data = spdf@data %>% mutate(google_name = gsub(" \\(United States\\)", "", google_name))
spdf_fortified <- tidy(spdf, region = "google_name")

# Calculate the centroid of each hexagon to add the label:
centers <- as.data.table(cbind.data.frame(data.frame(gCentroid(spdf, byid=TRUE), id=spdf@data$iso3166_2)))

spdf_fortified <- merge(spdf_fortified, plot_data, by.x = "id", by.y="state_name")
spdf_fortified <- as.data.table(spdf_fortified)

# Now I can plot this shape easily as described before:
pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/hex_maps_speedfactor.pdf", width = 12, height = 10)
plot <- ggplot() +
  geom_polygon(data = spdf_fortified, aes( x = long, y = lat, group = group, fill=col_cat), color="white") +
  # geom_text(data=centers, aes(x=x, y=y, label=id), color = "white") +
  geom_text(data=centers, aes(x=x, y=y, label=id), color = "white") +
  theme_void() +
  coord_map() +
  # scale_fill_manual(values = c('#00579f', '#fbe1c8', '#ffd2a0', '#ffc479', '#ffb54e', '#ffa500', '#C1C1C1')) +
  scale_fill_manual(values = c('#8fa0c8', '#f0b4a0', '#e99579', '#e07754', '#d4562f', '#c63003', '#C1C1C1')) +
  # , '#f2d2c8', '#f0b4a0', '#e99579', '#e07754', '#d4562f', '#c63003' ## '#8fa0c8'
  # '#8fa0c8', '#f2d2c8', '#f0b4a0', '#e99579', '#e07754', '#d4562f', '#c63003'
  # scale_fill_viridis_c(option = "A", na.value = "gray") +
  # scale_fill_viridis_c(option = "A", breaks = c(-17, -15, -10, -5, 0, 5, 10, 12), na.value = "gray") +
  theme(legend.position = "bottom", legend.key.width=unit(1.8,"cm")) +
  guides(fill=guide_legend(nrow=4, byrow = T)) +
  labs(fill = "", title = "Pace Required to Reach 70% Coverage Among Eligible (12+) by July 4, by Race/Ethnicity\n") +
  #caption = "*Gray colors do not have sufficient non-white population size") + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(strip.text = element_text(size = 12), text = element_text(size = 14)) +
  facet_wrap(~race_grp)
print(plot)
dev.off()

## Fill out table
plot_data <- copy(vax_stats)
plot_data <- plot_data[Date=="2021-06-07"]
plot_data <- plot_data[, accel:=((70-vaccinated_pct_12)/as.numeric((ymd("2021-07-04")-ymd("2021-06-07"))))/100]
plot_data <- plot_data[, accel_agg:=((70-agg)/as.numeric((ymd("2021-07-04")-ymd("2021-06-07"))))/100]

accel_df <- unique(plot_data[,.(accel, accel_agg, race_grp, diff_rate, diff_all_rate, state_name)])

accel_df <- accel_df[, speed_factor_race:=round(((accel-diff_rate)/diff_rate)*100, digits = 0)]
accel_df <- accel_df[, speed_factor_:=as.character(round(((accel-diff_rate)/diff_rate)*100, digits = 0))]

accel_df <- accel_df[race_grp%in%c("White", "Black", "Hispanic", "Asian")]
accel_df <- accel_df[speed_factor_race<=0, speed_factor_:="On Track"]
accel_df <- accel_df[speed_factor_race>500, speed_factor_:="500+"]
accel_df <- merge(accel_df, suppress_table, by = c("state_name", "race_grp"))
accel_df <- accel_df[ne==1, speed_factor_:="NE"]
accel_df <- accel_df[nr==1, speed_factor_:="NR"]
accel_df <- accel_df[, c("accel", "diff_rate", "speed_factor_race", "ne", "nr"):=NULL]

accel_df <- dcast(accel_df, state_name +  diff_all_rate + accel_agg ~ race_grp, value.var = c("speed_factor_"))

accel_df <- accel_df[, temp:=round(((accel_agg-diff_all_rate)/diff_all_rate)*100, digits = 0)]
accel_df <- accel_df[, All:=as.character(round(((accel_agg-diff_all_rate)/diff_all_rate)*100, digits = 0))]
accel_df <- accel_df[temp<=0, All:="On Track"]
accel_df <- accel_df[temp>500, All:="500+"]
accel_df <- accel_df[,.(state_name, All, Asian, Black, Hispanic, White)]

## National acceleration
plot_data <- copy(vax_stats_nat)
plot_data <- plot_data[is.na(project)]
plot_data <- plot_data[order(race_grp, Date)]
plot_data <- plot_data[, l2:=shift(nat_race, n = 2, type = "lag"), by = c("race_grp")]
plot_data <- plot_data[, diff:=nat_race-l2]

plot_data <- plot_data[, l2_all:=shift(nat_race_agg, n = 2, type = "lag"), by = c("race_grp")]
plot_data <- plot_data[, diff_all:=nat_race_agg-l2_all]

# Per population
plot_data <- plot_data[, diff_rate:=(diff/(14))]
plot_data <- plot_data[, diff_all_rate:=(diff_all/(14))]

plot_data <- plot_data[Date=="2021-06-07"]
plot_data <- plot_data[, accel:=((70-nat_race)/as.numeric((ymd("2021-07-04")-ymd("2021-06-07"))))]
plot_data <- plot_data[, accel_agg:=((70-nat_race_agg)/as.numeric((ymd("2021-07-04")-ymd("2021-06-07"))))]

accel_df_nat <- unique(plot_data[,.(accel, accel_agg, race_grp, diff_rate, diff_all_rate)])

accel_df_nat <- accel_df_nat[, speed_factor_race:=round(((accel-diff_rate)/diff_rate)*100, digits = 0)]
accel_df_nat <- accel_df_nat[, speed_factor_:=as.character(round(((accel-diff_rate)/diff_rate)*100, digits = 0))]

accel_df_nat <- accel_df_nat[race_grp%in%c("White", "Black", "Hispanic", "Asian")]
accel_df_nat <- accel_df_nat[speed_factor_race<=0, speed_factor_:="On Track"]
accel_df_nat <- accel_df_nat[, c("accel", "diff_rate", "speed_factor_race"):=NULL]

accel_df_nat <- dcast(accel_df_nat, diff_all_rate + accel_agg ~ race_grp, value.var = c("speed_factor_"))

accel_df_nat <- accel_df_nat[, temp:=round(((accel_agg-diff_all_rate)/diff_all_rate)*100, digits = 0)]
accel_df_nat <- accel_df_nat[, All:=as.character(round(((accel_agg-diff_all_rate)/diff_all_rate)*100, digits = 0))]
accel_df_nat <- accel_df_nat[temp<=1, All:="On Track"]
accel_df_nat <- accel_df_nat[,.(All, Asian, Black, Hispanic, White)]
accel_df_nat <- accel_df_nat[, state_name:="United States"]

accel_df <- rbind(accel_df_nat, accel_df)
accel_df <- accel_df[,.(state_name, All, Asian, Black, Hispanic, White)]

write.csv(accel_df, "~/COVID/vaccination_analysis/KFF Analysis/Final/output/acceleration_df.csv", na = "", row.names = F)


vax_stats <- copy(backup)

for (i in ymd("2021-06-07"):ymd("2021-08-31")) {
  projection <- vax_stats[Date==i & !(state_name%in%c("Nebraska"))]
  # Per population
  projection <- projection[, race_doses:=race_doses+(diff_rate*0.75*(pop_12))]
  projection <- projection[, doses_adj:=doses_adj+(diff_all_rate*0.75*(pop_12_allrace))]
  
  projection <- projection[, Date:=Date+1]
  projection <- projection[, project:=1]
  vax_stats <- rbind(vax_stats, projection, fill = T)
}

for (i in ymd("2021-05-24"):ymd("2021-08-31")) {
  projection <- vax_stats[Date==i & (state_name%in%c("Nebraska"))]
  # Per population
  projection <- projection[, race_doses:=race_doses+(diff_rate*0.75*(pop_12))]
  projection <- projection[, doses_adj:=doses_adj+(diff_all_rate*0.75*(pop_12_allrace))]
  
  projection <- projection[, Date:=Date+1]
  projection <- projection[, project:=1]
  vax_stats <- rbind(vax_stats, projection, fill = T)
}

vax_stats <- vax_stats[, vaccinated_pct_12:=(race_doses/pop_12)*100]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), race_doses:=pop_12]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), vaccinated_pct_12:=100]

vax_stats <- vax_stats[, tot_doses:=sum(race_doses, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[, agg:=(tot_doses/pop_12_allrace)*100]


pdf("~/COVID/vaccination_analysis/KFF Analysis/Final/figures/projections_state_declining_pace.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & is.na(project) & pop_12>=50000],
                 aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & is.na(project)], aes(y = agg), size = 1, color = "dark gray") +
    geom_point(size = 2) + geom_line( size = 1.1) +
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & (project == 1 | Date == "2021-06-07") & Date <="2021-07-04"], aes(y = agg), size = 1.2, color = "dark gray",
              linetype = "dashed", alpha = 0.9) +
    geom_line(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & (project == 1 | Date == "2021-06-07") & pop_12>=50000 & Date <="2021-07-04"],
              size = 1.3,
              linetype = "dashed", alpha = 0.9) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"),
                       values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage Among Eligible (12+)", color = "Race/Ethnicity", title = paste0(i),
         caption = "Gray line shows coverage among whole eligible (12+) population.\nSolid lines show historical (observed) data, dashed lines show projections.") +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 70, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 80)) +
    scale_shape_manual(values = c(19, 1))
  print(plot)
}
dev.off()

vax_stats_nat <- copy(vax_stats)
vax_stats_nat <- vax_stats_nat[, nat_doses_race:=sum(race_doses, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_race:=sum(pop_12, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_doses_agg:=sum(race_doses, na.rm=T), by = c("Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_agg:=sum(pop_12, na.rm=T), by = c("Date")]

vax_stats_nat <- vax_stats_nat[, nat_race:=(nat_doses_race/nat_pop_race)*100]
vax_stats_nat <- vax_stats_nat[, nat_race_agg:=(nat_doses_agg/nat_pop_agg)*100]

vax_stats_nat <- unique(vax_stats_nat[,.(race_grp, Date, nat_race, nat_race_agg, project)])

write.csv(vax_stats_nat[Date=="2021-07-04"], "~/COVID/vaccination_analysis/KFF Analysis/Final/output/national_declining.csv", na = "", row.names = F)


