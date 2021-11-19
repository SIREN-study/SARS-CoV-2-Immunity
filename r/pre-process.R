#
# Pre-process R file
# Imports and cleans the SIREN data and generates datasets
#

# libraries
library(readstata13)
library(tidyverse)
library(janitor)

# load data
siren_raw <- read.dta13("~/coviddata/SIREN_Interim_20210921_v2.dta") %>%
    clean_names() %>%
    rename(reinfection_pcr_date = reinfection_pc_rdate,
           primary_pcr_date = primary_pc_rdate,
           last_pcr_neg_date = last_pc_rneg_date
           )

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
        siren_cohort %>% filter(cohort_final==0, !is.na(primary_pcr_date),
                                (primary_pcr_date<vaccine_date1 | is.na(vaccine_date1)),
                                ar<primary_pcr_date) %>%
        mutate(cohort_final = 1,
               study_id = paste0(study_id,"_2"),
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
    mutate(
        start_date_pos_c = start_date_pos_c - start_time,
        time= time - start_time,
        ar = ar - start_time,

        vd1= vaccine_date1 - start_time,
        vd2 = vaccine_date2 - start_time,

        follow_up_time=time,

        region=factor(region),
        agegr2 = factor(agegr2),
        gender=factor(gender),
        ethnic_gr=factor(ethnic_gr),
        work_exposure_frequency = factor(work_exposure_frequency),
        occ_set_cat = factor(occ_set_cat)
    )

siren %>% count() # n = 36668 records

saveRDS(siren, "~/coviddata/siren_pre_processing.RDS")
