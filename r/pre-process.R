#
# Pre-process R file
# Imports and cleans the SIREN data and generates datasets
#

here::i_am("r/pre-process.R")

# libraries
library(readstata13)
library(tidyverse)
library(janitor)
library(here)

# load data
siren_raw <- read.dta13(here("data/SIREN_Interim_20210921_v2.dta")) %>%
    clean_names() %>%
    rename(
        reinfection_pcr_date = reinfection_pc_rdate,
        primary_pcr_date = primary_pc_rdate,
        last_pcr_neg_date = last_pc_rneg_date
    )

age_group <- read.dta13(here("data/AgeGroup_LookUp.dta")) %>%
    clean_names()

# determine the start date for the analysis
start_time <- min(siren_raw$vaccine_date1, na.rm = TRUE) - 1 # this date is 7th December

# remove records with inconsistencies from the SIREN data
siren_cohort <- siren_raw %>%
    mutate(
        ar = pmax(date_enrolled, start_time)
    ) %>%
    filter(
        (vaccine_name1 <= 2 | is.na(vaccine_name1)),
        (vaccine_name2 <= 2 | is.na(vaccine_name2)),
        (vaccine_name1 == vaccine_name2 | is.na(vaccine_name1) | is.na(vaccine_name2)),
        (is.na(reinfection_pcr_date) | reinfection_pcr_date >= ar),
        !(start_date_pos_c > vaccine_date1 & !is.na(start_date_pos_c) & !is.na(vaccine_date1))
    )

# generate the time and event fields for each participant
siren <- siren_cohort %>%
    # allow people to re-enter the cohort after primary infection
    bind_rows(
        siren_cohort %>% filter(
            cohort_final == 0, !is.na(primary_pcr_date),
            (primary_pcr_date < vaccine_date1 | is.na(vaccine_date1)),
            ar < primary_pcr_date
        ) %>%
            mutate(
                cohort_final = 1,
                study_id = paste0(study_id, "_2"),
                ar = primary_pcr_date,
                start_date_pos_c = primary_pcr_date
            )
    ) %>%
    mutate(
        event = case_when(
            !is.na(primary_pcr_date) & primary_pcr_date > ar & cohort_final == 0 ~ 1,
            !is.na(reinfection_pcr_date) & reinfection_pcr_date >= ar & reinfection_pcr_date >= (start_date_pos_c) & cohort_final == 1 ~ 1,
            TRUE ~ 0
        ),
        time = case_when(
            event == 1 & cohort_final == 0 ~ primary_pcr_date,
            event == 1 & cohort_final == 1 ~ reinfection_pcr_date,
            event == 0 ~ last_pcr_neg_date
        ),
        dose_schedule = case_when(
            dose_schedule == 2 ~ 0,
            dose_schedule == 3 ~ 1
        )
    ) %>%
    filter(
        !is.na(time),
        time > ar
    )

# generate relative time variables to use in Cox proportional hazards model
siren <- siren %>%
    left_join(age_group %>% select(study_id, linked_age = agegr2), by="study_id") %>%
    mutate(
        start_date_pos_c = start_date_pos_c - start_time,
        time = time - start_time,
        ar = ar - start_time,
        vd1 = vaccine_date1 - start_time,
        vd2 = vaccine_date2 - start_time,
        follow_up_time = time,
        region = factor(region, labels = c("East Midlands","East of England","London","North East",
                                           "North West","South East","South West","West Midlands",
                                           "Yorkshire and Humber","Scotland","Northern Ireland","Wales")),
        linked_age = factor(linked_age, labels = c("<25","25-34","35-44","45-54","55-64","65+")),
        agegr2 = factor(agegr2, labels = c("<25","25-34","35-44","45-54","55-64","65+","Unknown age")),
        agegr2 = if_else(agegr2=="Unknown age", linked_age, agegr2),
        gender = factor(gender, labels = c("Male","Female","Non-binary","Unknown gender")),
        ethnic_gr = factor(ethnic_gr, labels = c("White","Mixed","Asian","Black","Other ethnicity","Unknown ethnicity")),
        work_exposure_frequency = factor(work_exposure_frequency, labels = c("Every day","Once a week","Once a month","Less than once a month","Never")),
        occ_set_cat = factor(occ_set_cat, labels = c("Office","Patient facing (non-clinical)","Outpatient",
                                                     "Maternity,Labour Ward","Ambulance,Emergency Department,Inpatient",
                                                     "Intensive care","Theatres","Other"))
    )

siren %>% count() # n = 36668 records

saveRDS(siren, "~/coviddata/siren_pre_processing.RDS")