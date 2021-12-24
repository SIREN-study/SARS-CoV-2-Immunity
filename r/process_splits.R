#
# Process splits R file
# Uses the survival library and tmerge function to split data into time chunks for Cox proportional hazards models
#

# libraries
library(tidyverse)
library(survival)

# load data
siren <- readRDS("~/coviddata/siren_pre_processing.RDS")
len <- count(siren)

# Vaccine effectiveness dataset

## time splits:
## 1. pre-vaccination
## 2. at(20 27 41 55) after(vd1)
## 3. at(0 13 73 133 193) after(vd2)
## at(0 90 273 454) after(primary_inf) for regression

split_ve <- tmerge(
    data1 = siren,
    data2 = siren,
    id = study_id,
    tstop = time,
    event = event(time, event),
    eligible = tdc(ar),
    primary_inf = cumtdc(start_date_pos_c),
    primary_inf = cumtdc(start_date_pos_c + 90),
    primary_inf = cumtdc(start_date_pos_c + 273),
    primary_inf = cumtdc(start_date_pos_c + 454),
    vaccine = event(vd1, rep(1, len)),
    vaccine = event(if_else(vd1 + 20 < vd2 & !is.na(vd2), vd1 + 20, NA_real_), rep(2, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 20, NA_real_), rep(2, len)), # to deal with the vd2==NA case
    vaccine = event(if_else(vd1 + 27 < vd2 & !is.na(vd2), vd1 + 27, NA_real_), rep(3, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 27, NA_real_), rep(3, len)), # to deal with the vd2==NA case
    vaccine = event(if_else(vd1 + 41 < vd2 & !is.na(vd2), vd1 + 41, NA_real_), rep(4, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 41, NA_real_), rep(4, len)), # to deal with the vd2==NA case
    vaccine = event(if_else(vd1 + 55 < vd2 & !is.na(vd2), vd1 + 55, NA_real_), rep(5, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 55, NA_real_), rep(5, len)), # to deal with the vd2==NA case
    vaccine = event(vd2, rep(6, len)),
    vaccine = event(vd2 + 13, rep(7, len)),
    vaccine = event(vd2 + 73, rep(8, len)),
    vaccine = event(vd2 + 133, rep(9, len)),
    vaccine = event(vd2 + 193, rep(10, len))
) %>%
    group_by(study_id) %>%
    # carry over the 0 vaccine status
    mutate(
        vaccine = lag(vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(is.na(vaccine), 0, vaccine)
    ) %>%
    ungroup() %>%
    mutate(
        vaccine_cat = case_when(
            # Pfizer
            vaccine == 1 & vaccine_name1 == 1 ~ 1,
            vaccine == 2 & vaccine_name1 == 1 ~ 2,
            vaccine == 3 & vaccine_name1 == 1 ~ 3,
            vaccine == 4 & vaccine_name1 == 1 ~ 4,
            vaccine == 5 & vaccine_name1 == 1 ~ 5,
            # ChAdOx
            vaccine == 1 & vaccine_name1 == 2 ~ 6,
            vaccine == 2 & vaccine_name1 == 2 ~ 7,
            vaccine == 3 & vaccine_name1 == 2 ~ 8,
            vaccine == 4 & vaccine_name1 == 2 ~ 9,
            vaccine == 5 & vaccine_name1 == 2 ~ 10,
            # Pfizer long
            vaccine == 6 & vaccine_name2 == 1 & dose_schedule == 1 ~ 11,
            vaccine == 7 & vaccine_name2 == 1 & dose_schedule == 1 ~ 12,
            vaccine == 8 & vaccine_name2 == 1 & dose_schedule == 1 ~ 13,
            vaccine == 9 & vaccine_name2 == 1 & dose_schedule == 1 ~ 14,
            vaccine == 10 & vaccine_name2 == 1 & dose_schedule == 1 ~ 15,
            # Pfizer short
            vaccine == 6 & vaccine_name2 == 1 & dose_schedule == 0 ~ 16,
            vaccine == 7 & vaccine_name2 == 1 & dose_schedule == 0 ~ 17,
            vaccine == 8 & vaccine_name2 == 1 & dose_schedule == 0 ~ 18,
            vaccine == 9 & vaccine_name2 == 1 & dose_schedule == 0 ~ 19,
            vaccine == 10 & vaccine_name2 == 1 & dose_schedule == 0 ~ 20,
            # ChAdOx
            vaccine == 6 & vaccine_name2 == 2 ~ 21,
            vaccine == 7 & vaccine_name2 == 2 ~ 22,
            vaccine == 8 & vaccine_name2 == 2 ~ 23,
            vaccine == 9 & vaccine_name2 == 2 ~ 24,
            vaccine == 10 & vaccine_name2 == 2 ~ 24,
            # Pre-vaccination
            TRUE ~ vaccine
        ),
        vaccine_cat = factor(vaccine_cat, levels = c(0:24)),
        primary_inf = factor(primary_inf)
    ) %>%
    mutate(
        # don't consider records with an infection within 90 days
        eligible = if_else(primary_inf == 1, as.integer(0), eligible)
    )

# check categories correctly assigned
split_ve %>%
    filter(eligible == 1, primary_inf == 0) %>%
    distinct(study_id, vaccine_cat) %>% # uncomment for participant numbers
    select(vaccine_cat) %>%
    table(useNA = "ifany")


# Durability of protection dataset

## time splits:
## 1. pre-vaccination, dose 1, dose 2
## 2. at(0 90 365) after(primary_inf)

split_dura <- tmerge(
    data1 = siren,
    data2 = siren,
    id = study_id,
    tstop = time,
    event = event(time, event),
    eligible = tdc(ar),
    primary_inf = cumtdc(start_date_pos_c),
    primary_inf = cumtdc(start_date_pos_c + 90),
    primary_inf = cumtdc(start_date_pos_c + 365),
    vaccine = event(vd1, rep(1, len)),
    vaccine = event(if_else(vd1 + 20 < vd2 & !is.na(vd2), vd1 + 20, NA_real_), rep(2, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 20, NA_real_), rep(2, len)), # to deal with the vd2==NA case
    vaccine = event(if_else(vd1 + 80 < vd2 & !is.na(vd2), vd1 + 80, NA_real_), rep(3, len)),
    vaccine = event(if_else(is.na(vd2), vd1 + 80, NA_real_), rep(3, len)), # to deal with the vd2==NA case
    vaccine = event(vd2, rep(4, len)),
    vaccine = event(vd2 + 13, rep(5, len)),
    vaccine = event(vd2 + 73, rep(6, len)),
    vaccine = event(vd2 + 133, rep(7, len)),
    vaccine = event(vd2 + 193, rep(8, len))
) %>%
    group_by(study_id) %>%
    # carry over the 0 vaccine status
    mutate(
        vaccine = lag(vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(vaccine == 0, lag(vaccine), vaccine),
        vaccine = if_else(is.na(vaccine), 0, vaccine),
        # remove the third vaccine state
        vaccine = if_else(vaccine == 3, NA_real_, vaccine),
        vaccine = if_else(vaccine>3, vaccine-1, vaccine),
        primary_cat = (primary_inf * 8) + vaccine,
        primary_cat = if_else(primary_cat>15, primary_cat-8, primary_cat),
        primary_cat = factor(primary_cat, levels = c(0:23))
    ) %>%
    ungroup() %>%
    mutate(
        # don't consider records with an infection within 90 days
        eligible = if_else(primary_inf == 1, as.integer(0), eligible)
    )

# check categories correctly assigned
split_dura %>%
    filter(eligible == 1, (is.na(vaccine_name2) | vaccine_name2!=2)) %>%
    distinct(study_id, primary_cat) %>%
    select(primary_cat) %>%
    table(useNA = "ifany")

save(split_ve, split_dura, siren, file = "~/coviddata/siren_post_processing.RData")