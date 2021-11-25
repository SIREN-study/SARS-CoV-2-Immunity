* Interim analysis 2, November 2021. 
*Ferdinando Insalata, ferdinando.insalata@phe.gov.uk
*18/11/2021

clear
set more off

use ../With_Pos/SIREN_Interim_20210921_v2, clear

*Add missing values of variable AGEGR2 (grouping by age)
merge 1:m StudyId   using ../With_Pos/AgeGroup_LookUp
drop if _merge==2
replace AGEGR2 = AGEGR2_ if AGEGR2 ==99
drop _merge


*Define December 7th as start_time (day before first jab administered)
egen start_time=min(vaccine_date1-1)

*Label the variable cohort_final that assigns cohort at Dec 7th.
label define cohort_label 0 "Negative" 1 "Positive"
label values cohort_final cohort_label

**TO REMOVE IN GITHUB VERSION
*Use this to start everyone at enrollment.
*egen A = min(Date_Enrolled)
*replace start_time=A
*replace cohort_final = cohort

*Keep only participants who are vaccinated with Pfixer or ChAdOx (Pfizer =1 and ChAdOX=2, >2 is unknown or other)
drop if vaccine_name1!=. & vaccine_name1>2
drop if vaccine_name2!=. & vaccine_name2>2

*Keep only participants who received two doses of the vaccine brand
drop if vaccine_name1==1 & vaccine_name2==2
drop if vaccine_name1==2 & vaccine_name2==1
 
*Drop people that before enroll are vaccinated and then infected,
*since we only investigate immunity from infection before vaccination for now.
drop if start_date_posC >vaccine_date1   & start_date_posC !=. & vaccine_date1   !=.

*generate a variable, named "ar", containing the day a person becomes at risk (in this analysis it is the latest between Enrolment and Dec 7th)
gen ar=max(Date_Enrolled, start_time)

*Drop people with two infections before at risk, as in this analysis we investigate protection from one infection only . 
drop if ReinfectionPCRdate < ar

*Recode and relabel the values of DoseSchedule, as the values of DoseSchedule=1 have been found to be wrong
replace DoseSchedule=. if DoseSchedule==1
recode DoseSchedule 2=0 3=1
label define relabel_schedules 0 "Short Schedule" 1 "Long Schedule"
label values DoseSchedule relabel_schedules



*************************GENERATE TIME AND EVENT VARIABLES*************************

*First we focus on people who are in the negative cohort (no evidence of previous infection) on December 7th.
*Among these poeple, those with a PCR-confirmed primary infection after the time they are at risk are assigned event=1.
generate event=1 if PrimaryPCRdate !=. & PrimaryPCRdate > ar & cohort_final==0

*And the data of this primary infection is the time-to-event for them.
gen time = PrimaryPCRdate if event==1


*Now we focus on people who are in the positive cohort (evidence of previous infection, either PCR positive or Covid symptoms + subsequent antibody positive result) on December 7th.

*Among these people, those with a PCR-confirmed reinfection are assigned event=1
replace event=1 if ReinfectionPCRdate !=. & ReinfectionPCRdate >= ar & cohort_final==1

*And the date of this PCR reinfection is their time-to-event.
replace time = ReinfectionPCRdate if event==1 & cohort_final==1 & ReinfectionPCRdate >= ar

*The people whose event variable is not 1 (i.e. no primary infection or reinfection in follow up) are assigned event=0.
replace event=0 if event==.

*TO ERASE IN GITHUB VERSION
*NOTE: I tried to do the following:
*replace time= LastABneg_date if (event==0 & time==.)
*6 people are concercned, but the lastAB negative is before infection, only one is after but it is in the positive cohort so we would need LastPCRneg_date (which is missing).
*Therefore these three people will be dropped anyway because time would be < ar if we did the replacement.

*And they are censored at the time of their latest negative PCR test.
replace time = LastPCRneg_date if event==0 

*Drop those people with no test or for which the only test is before they were at risk (ar variable).
drop if time==.
drop if time < = ar 

*************************SPLITTING NEGATIVE COHORT AT FAILURE*************************

*The aim of this is allowing people infected during follow up to contribue also to the positive cohort.
*We only do this for people with a primary infection after they become at risk (ar variable) and before first dose of vaccine.
*This is because --in this analysis-- we want to investigate the protection from vaccine and/or from a single infection before vaccine  

*Selec the people concerned and label them with the variable ibv=1 (ibv stands for infected before vaccination).								
gen ibv= 1 if cohort_final ==0 & ar < PrimaryPCRdate & PrimaryPCRdate < vaccine_date1 
replace ibv=0 if ibv==.

*For the people concerned, their time-to-event is now reinfection date (if they are reinfected) or latest PCR negative (after the primary infection).
*See below to understand why this replacement is needed.
replace time = ReinfectionPCRdate if ibv ==1 & ReinfectionPCRdate !=.
replace time = LastPCRneg_date if ibv==1 & ReinfectionPCRdate ==. & LastPCRneg_date >PrimaryPCRdate & LastPCRneg_date !=.

*Declare the data as survival data. This is needed to use the stsplit command below, which allow people to move cohorts.
stset time, id(StudyId) failure(event) origin(time ar)

*If people have ibv=1 (i.e. can contribute to both negative and positive cohort), their observation is splitted into two observations on the day of their primary infection.
*the variable "after_primary" indicated which pseudo-observation, i.e. after_primary=1 for the follow-up time after primary infection.

stsplit after_primary if(ibv==1), at(0) after(PrimaryPCRdate)
replace after_primary=after_primary+1
*IMPORTANT: The splitting command also changed the time-to-event ("time" variable) of the first pseudo-observation (that with after_primary=0) to the argument of after().
*The argument of after() is the date of Primary infection. The first pseudo-observation ends on the date of the primary infection. 
*In conclusion, the "time" variable of the first pseudo-observation (referring to before primary infection, encoded by after_primary=0) is correctly the day of primary infection.

*The pseudo-observation after infection is labelled as belonging to the positive cohort.
replace cohort_final=1 if after_primary==1

*The start date ("start_date_posC") and the time at risk ("ar") of the pseudo-observation after infection becomes the day of primary infection.
replace start_date_posC = PrimaryPCRdate if after_primary ==1
replace ar=PrimaryPCRdate if after_primary==1

*The second pseudo observation has event=. after the splitting. It must be one if reinfection happens, or 0 otherwise (if the time variable is that of the latest negative PCR).
replace event =1 if event==.
replace event=0 if after_primary ==1 & time==LastPCRneg_date


*Drop people whose latest test is before they become at risk (_st==0)
*Drop the variable created by stset that we do not need. 
drop if _st==0
drop _t _t0 _st _d


*Save file
save first_pre_processing, replace


exit




