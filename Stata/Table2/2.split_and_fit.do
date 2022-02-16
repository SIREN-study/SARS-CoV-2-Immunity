clear
set more off

use first_pre_processing, clear

*Set December 7th to be time zero, start of analysis time. Dates are transformed into number of days since December 7th.
*Some variables are copied because it is useful to keep them as dates for future visualisations.
gen start_date_posC_2 = start_date_posC - start_time
gen follow_up_time = time - ar	
replace time = time - start_time
replace PrimaryPCRdate=PrimaryPCRdate - start_time +0

replace ar = ar - start_time

*The variables vd1 and vd2  are needed for the splitting below (lines 40 and 67).
gen vd1= vaccine_date1 - start_time
gen vd2 = vaccine_date2 - start_time
*hence they are also copied (as days from De 7th) because the splitting will modify the values of vd1 and vd2.
gen vax_date2=vaccine_date2 - start_time
gen vax_date1=vaccine_date1- start_time

*************************SPLITTING TIME BETWEEN JABS********************
*In order to split the time between jabs, we create a fake failure variable (fail) always equal to 1. The time of failure is the day of the second jab or,
*if the participant did not receive a second jab, it will coincide with the "time" variable (either infection or censoring time).
*This approach has "side effects": it creates spurious pseudo observations and the variable event is assigned missing (.) value when it should be zero.
*These are taken care of at the end.

					
												
*The below line creates a fake failure variable, to split time between jabs. Time of failure will be the date of second jab.
gen fail=1

*The below line adapts the strategy to people who have not received second jab or leave the cohort before first dose.
replace vd2=time if (vd2==. | time < vd1)

*For STATA to accept the data as survival data, the pseudo-observation after primary infection of people who move chort (see pre_process.do) must have a different StudyId.
*Without the following line, Stata gives an error because of a clash between start time and failure time of pseudo-observations with same StudyId.
replace StudyId =StudyId+"pos" if after_primary ==1

*stset and split, to split time between jabs.

stset vd2, id(StudyId)  fail(fail) origin(time 0) 
stsplit d1, at(0 20 27 41 55 79) after(vd1)

*Clean side effects of the previous stset and stsplit:

*The below line deletes pseudo-observations which start beyond second dose. 
*We do not consider the time since first dose when a participant has received thei second dose.
drop if _t0>=time

*The below line deals with people without first jab
replace d1=-1 if vaccine_date1 ==. & d1!=-1

*Clean the rest
sort StudyId
by StudyId: generate n1 = _n
by StudyId: generate N = _N
replace time=_t if n1!=N
replace event=0 if n1!=N
gen Enrolment= Date_Enrolled - start_time
drop fail N n1 _t0 _t _st 

*************************SPLITTING TIME AFTER SECOND JAB*************************
					
stset time, id(StudyId) failure(event) origin(time 0)
stsplit d2, at(0 13 73 133 193) after(vd2)

*Label d1=.  the pseudo-observations after second dose. This is not strictly needed, done for convenience.
replace d1=. if d2>=0

*************************TIME AFTER PRIMARY INFECTION (for positive cohort only)*************************
/*
We also want to investigate the protection conferred by primary infection. Therefore we further split into pseudo-observations
according to the time since primary infection (if any). 
In this model, the time since primary infection does not interact with vaccination status. 
The other Stata code(split_all_comb) refers to the model with combinations between vaccine status and time since primary infection.
*/
	
stsplit Primary_inf if(cohort_final==1), at(0 90 273 454) after(start_date_posC_2)

*Next line drops the first 90 days after a primary infection, since by definition of reinfection (two positive PCR at least 90 days apart)
*the participant is not at risk during that time.
drop if Primary_inf==0

*The above stsplit   creates spurious pseudo-observations, starting and ending on the same day. The next line drops them:
drop if after_primary ==1 & ibv ==1 & cohort_final ==1 & d1==-1 & ar==time & Primary_inf ==-1

*Next line labels the pseudo observation   and observations   where there is no previous primary infection with value Primary_inf=0.
*"Primary_inf" is currently -1 for the pseudo-observation before Primary infection. The value 0 is more logical.
*Also other values of Primary_inf are affected, but they are recoded below  .
replace Primary_inf = Primary_inf + 1  

*Recode and label the time chunks after primary infection
recode Primary_inf 0=0 91=1 274=2 455=3  
label define prim_inf_lab 0 "No previous inf" 1 "3-9 mth" 2 "9-15_mth" 3 "15+_mth"
label values Primary_inf prim_inf_lab

*All pseudo-observations of people in the negative cohort have Primary_inf=".", as the splitting affected only people in the positive cohort.
*We replace this missing value with the value 0, so these participant contribute follow-up time before/without primary infection.
replace Primary_inf=0 if Primary_inf==.



*************************SPLIT IN BEFORE AND AFTER ENROLMENT*************************
*This is needed because some people enrolled after Dec 7th, and we need to tell Stata to discard the time between enrolment and Dec 7th.

stsplit is_enrolled, at(0) after(Enrolment)

*Label as is_enrolled=0 the pseudo-observations before enrolment, with 1 those after enrolment (currently -1 and 0).
replace is_enrolled=is_enrolled+1

*Before this line, _st==1. always. Next line tell stata not to use time before enrolment.
replace _st=is_enrolled

drop if _st==0
*The above line does not affect the HR estimates, only affects the calculation of exposure and rates per category.
*Without dropping _st==0, exposure that should not be considered is counted, and crude rates of infection are affected too.



*************************GROUPING INTO CATEGORIES AND LABELLING***********************
/*
The variable x indicates vaccine status of the pseudo-observation, hence it is a time-varying covariate.
We encode into x number of doses (1 or 2), manufacturer, dosing interval (for Pfizer only) and time since latest dose

The value x=0 refers to unvaccinated, and labels the pseudo-observations referring to follow-up time before vaccination,
or the entire observation if the person stays unvaccinated during all their follow-up.

Then, for example, x=9 indicates that a pseudo-observation refers to follow-up interval betwene 42 and 55 days(both extremes included) after 
the participant has received the first dose of ChAdOx, provided that this is before the second dose. This means that when a participant receives 
their second dose (if they do), the time since first dose is no longer taken into account.
Other example: x=13 refers to a participant who received two doses of Pfizer, with a long (>=6 weeks) dosing interval between jabs,
between 74 and 133 days after second dose.
*/

*Unvaccinated
generate x=0 if d1==-1



*Chunks of follow-up (pseudo-observations) between first and second dose.

*Pfizer
replace x=1 if d1==0        		& vaccine_name1==1 
replace x=2 if d1==20      			& vaccine_name1==1 
replace x=3 if d1==27      			& vaccine_name1==1  
replace x=4 if d1==41      			& vaccine_name1==1   
replace x=5 if (d1==55 | d1==79)	& vaccine_name1==1  

*Chadox
replace x=6  if d1==0       		& vaccine_name1==2 
replace x=7  if d1==20        		& vaccine_name1==2 
replace x=8  if d1==27      		& vaccine_name1==2  
replace x=9 if d1==41      			& vaccine_name1==2   
replace x=10 if (d1==55 | d1==79) 	& vaccine_name1==2   



*Chunks of follow-up (pseudo-observations) after second dose:

*Pfizer two doses, long dosing interval (>=6 weeks betwene doses).
replace x= 11 if d2==0      & vaccine_name2==1   & DoseSchedule==1
replace x= 12 if d2==13     & vaccine_name2==1   & DoseSchedule==1
replace x= 13 if d2==73     & vaccine_name2==1   & DoseSchedule==1
replace x= 14 if d2==133    & vaccine_name2==1   & DoseSchedule==1
replace x= 15 if d2==193    & vaccine_name2==1   & DoseSchedule==1

*Pfizer two doses, short dosing interval (<6 weeks betwene doses).
replace x= 16 if d2==0     & vaccine_name2==1   & DoseSchedule==0
replace x= 17 if d2==13    & vaccine_name2==1   & DoseSchedule==0
replace x= 18 if d2==73    & vaccine_name2==1   & DoseSchedule==0
replace x= 19 if d2==133   & vaccine_name2==1   & DoseSchedule==0
replace x= 20 if d2==193   & vaccine_name2==1   & DoseSchedule==0

*ChAdOx two doses (regardless of dosing interval).
replace x= 21 if d2==0     & vaccine_name2==2
replace x= 22 if d2==13    & vaccine_name2==2
replace x= 23 if d2==73    & vaccine_name2==2
replace x= 24 if (d2==133 | d2==193)   & vaccine_name2==2


*Labelling the values of x
label define xlab 0 "Unvaccinated" ///
1 "d1 0-20 PF" 2 "d1 21-27 PF" 3   "d1 28-41 PF" 4 "d1 42-55 PF" 5 "d1 56+ PF" ///
6 "d1 0-20 AZ" 7 "d1 21-27 AZ" 8   "d1 28-41 AZ" 9 "d1 42-55 AZ" 10 "d1 56+ AZ" ///
11 "d2 0-13 PF_long" 12 "d2 14-73 PF_long" 13 "d2 74-133 PF_long" 14 "d2 134-193 PF_long" 15 "d2 194+ PF_long" ///
16 "d2 0-13 PF_short" 17 "d2 14-73 PF_short" 18 "d2 74-133 PF_short" 19 "d2 134-193 PF_short" 20 "d2 194+ PF_short" ///
21 "d2 0-13 AZ" 22 "d2 14-73 AZ" 23 "d2 74-133 AZ" 24 "d2 134+ AZ"

*Drop "d2 0-13 AZ", as there are not enough people for that parameter estimate to converge. This does not affect results.
drop if x==21 //Comment this line for calculation of exposures.
label values x xlab
			


*People who move cohort have now a spurious pseudo-observation crated by the stset. This drops it.
drop if after_primary ==1 & ibv ==1 & cohort_final ==1 & x==0 & Primary_inf ==0

*Exposure is used to show the information when visualising splitted observations, and generate below tables.
generate exposure=_t - _t0

*Next replace tells Stata to consider only individuals without a primary infection. 
*Remove to fit on all individual and fit model of the first submission (add controlling on i.Primary_inf).
replace _st=0 if Primary_inf !=0

stcox i.x i.GENDER i.ETHNIC_GR, vce(cluster Trust_Code) nolog strata(AGEGR2 region WORK_EXPOSURE_FREQUENCY occ_set_cat) efron



