clear
set more off

use ../With_Pos/SIREN_Interim_20210921_v2, clear


merge 1:m StudyId   using ../With_Pos/AgeGroup_LookUp
drop if _merge==2

replace AGEGR2 = AGEGR2_ if AGEGR2 ==99


egen start_time=min(vaccine_date1-1)

label define cohort_label 0 "Negative" 1 "Positive"
label values cohort_final cohort_label

*Use this to start everyone at enrollment.
*egen A = min(Date_Enrolled)
*replace start_time=A
*replace cohort_final = cohort

drop if vaccine_name1!=. & vaccine_name1>2
drop if vaccine_name2!=. & vaccine_name2>2
drop if vaccine_name1==1 & vaccine_name2==2
drop if vaccine_name1==2 & vaccine_name2==1
 
*Drop people infected after vaccine, who then become at risk
count if start_date_posC >vaccine_date1   & start_date_posC !=. & vaccine_date1   !=.
*list StudyId Date_Enrolled vaccine_date1 vaccine_date2 PrimaryPCRdate start_date_posC if start_date_posC >vaccine_date1   & start_date_posC !=. & vaccine_date1   !=.
drop if start_date_posC >vaccine_date1   & start_date_posC !=. & vaccine_date1   !=.


*generate a variable telling the day a person becomes at risk IN THIS ANALYSIS (whatever is latest between Enrolment and Dec 7th)
gen ar=max(Date_Enrolled, start_time)

*Drop people with two infections before at risk, as we want to analyse protection from one infection. 
drop if ReinfectionPCRdate < ar

*Relabel the values of DoseSchedule
replace DoseSchedule=. if DoseSchedule==1
recode DoseSchedule 2=0 3=1
label define relabel_schedules 0 "Short Schedule" 1 "Long Schedule"
label values DoseSchedule relabel_schedules



				***GENERATE TIME AND EVENT VARIABLES***

*NEGATIVE COHORT_FINAL
generate event=1 if PrimaryPCRdate !=. & PrimaryPCRdate > ar & cohort_final==0
gen time = PrimaryPCRdate if event==1


*POSITIVE COHORT_FINAL
replace event=1 if ReinfectionPCRdate !=. & ReinfectionPCRdate >= ar & cohort_final==1
replace time = ReinfectionPCRdate if event==1 & cohort_final==1 & ReinfectionPCRdate >= ar

replace event=0 if event==.
*NOTE: I tried to do the following:
*replace time= LastABneg_date if (event==0 & time==.)
*6 people are concercned, but the lastAB negative is before infection, only one is after but it is in the positive cohort so we would need LastPCRneg_date (which is missing).
*Therefore these three people will be dropped anyway because time would be < ar if we did the replacement.
replace time = LastPCRneg_date if event==0 


drop if time==.
drop if time < = ar 

								*SPLITTING NEGATIVE COHORT AT FAILURE*
gen ibv= 1 if cohort_final ==0 & ar < PrimaryPCRdate & PrimaryPCRdate < vaccine_date1 
replace ibv=0 if ibv==.


***MY VERSION
replace time = ReinfectionPCRdate if ibv ==1 & ReinfectionPCRdate !=.
replace time = LastPCRneg_date if ibv==1 & ReinfectionPCRdate ==. & LastPCRneg_date >PrimaryPCRdate & LastPCRneg_date !=.

stset time, id(StudyId) failure(event) origin(time ar)
stsplit after_primary if(ibv==1), at(0) after(PrimaryPCRdate)
*For the pre-infection pseudo-episode time is set to PrimaryPCRdate by stsplit
replace after_primary=after_primary+1

replace cohort_final=1 if after_primary==1
replace start_date_posC = PrimaryPCRdate if after_primary ==1
replace ar=PrimaryPCRdate if after_primary==1

replace event =1 if event==.
replace event=0 if after_primary ==1 & time==LastPCRneg_date



drop if _st==0
drop _t _t0 _st _d



save first_pre_processing, replace


exit




