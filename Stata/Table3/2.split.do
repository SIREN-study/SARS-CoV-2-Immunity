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
stsplit d2, at(0 13 73 133 163 193) after(vd2)

*Label d1=.  the pseudo-observations after second dose. This is not strictly needed, done for convenience.
replace d1=. if d2>=0

*************************TIME AFTER PRIMARY INFECTION (for positive cohort only)*************************
/*
We also want to investigate the protection conferred by primary infection. Therefore we further split into pseudo-observations
according to the time since primary infection (if any). 
In this model, the time since primary infection does not interact with vaccination status. 
The other Stata code(split_all_comb) refers to the model with combinations between vaccine status and time since primary infection.
*/
	
stsplit Primary_inf if(cohort_final==1), at(0 90 365) after(start_date_posC_2)

*Next line drops the first 90 days after a primary infection, since by definition of reinfection (two positive PCR at least 90 days apart)
*the participant is not at risk during that time.
drop if Primary_inf==0

*The above stsplit creates spurious pseudo-observations, starting and ending on the same day. The next line drops them:
drop if after_primary ==1 & ibv ==1 & cohort_final ==1 & d1==-1 & ar==time & Primary_inf ==-1

*Next line labels the pseudo observation  and observations where there is no previous primary infection with value Primary_inf=0.
*"Primary_inf" is currently -1 for the pseudo-observation before Primary infection. The value 0 is more logical.
*Also other values of Primary_inf are affected, but they are recoded below.
replace Primary_inf = Primary_inf + 1  

*Recode and label the time chunks after primary infection
recode Primary_inf 0=0 91=1 366=2 
label define prim_inf_lab 0 "No previous inf" 1 "<1yr" 2 ">1yr" 
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

drop if _st==0 | vaccine_name2 ==2
*The above line does not affect the HR estimates, only affects the calculation of exposure and rates per category.
*Without dropping _st==0, exposure that should not be considered is counted, and crude rates of infection are affected too.


drop if after_primary ==1 & ibv ==1 & cohort_final ==1 & d1==-1 & Primary_inf ==0

			
*************************VISUALISATION*************************

*Exposure is used to show the information when visualising splitted observations, and generate below tables.
generate exposure=_t - _t0


*Generated only for visualisation
gen LastPCRneg= LastPCRneg_date - start_time
gen vax_interval = vax_date2 - vax_date1

*order StudyId Date_Enrolled Enrolment ar  whole_time start_date_posC PrimaryPCRdate ReinfectionPCRdate  LastPCRneg vax_date1 vax_date2 _t0 _t time x d1  d2 exposure _d  Primary_inf is_enrolled _st

save after_splitting, replace
