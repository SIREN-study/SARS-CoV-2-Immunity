clear
set more off

use first_pre_processing, clear

drop if start_date_posC >vaccine_date2  & start_date_posC !=. & vaccine_date2  !=.
drop if start_date_posC >vaccine_date1  & start_date_posC !=. & vaccine_date1  !=.

*Put start time as start of the clock
gen start_date_posC_2 = start_date_posC - start_time
replace time= time - start_time
replace ar = ar - start_time
gen vd1= vaccine_date1 - start_time
gen vd2 = vaccine_date2 - start_time
gen vax_date2=vaccine_date2 - start_time
gen vax_date1=vaccine_date1- start_time
*vax_date1/2 is repeated before vd1 and vd2 will be splitted.



							*SPLITTING TIME BETWEEN JABS*
*Generate a variable to store, before splitting happens,  information on the time from start_time to event.
gen whole_time=time							
						
							
							
*The below line creates a fake failure variable, to split time between jabs.
gen fail=1

*The below line adapts the strategy to people who have not received second jab.
replace vd2=time if (vd2==. | time < vd1)

replace StudyId =StudyId+"pos" if after_primary == 1
*stset and split, to split time between jabs.
stset vd2, id(StudyId)  fail(fail) origin(time 0) 
stsplit d1, at(0) after(vd1)


*The below line deletes spurious pseudo-observations, that start beyond second dose.
drop if _t0>=time

*The below line deals with people without first jab
replace d1=-1 if vaccine_date1 ==. & d1!=-1

*Clean side effects of the previous stset and stsplit
sort StudyId
by StudyId: generate n1 = _n
by StudyId: generate N = _N
replace time=_t if n1!=N
replace event=0 if n1!=N
list StudyId vd1 _t0 _t  time vd2 d1 _st if StudyId=="R0A10032"
gen Enrolment= Date_Enrolled - start_time
drop fail N n1 _t0 _t _st 




						*SPLITTING TIME AFTER SECOND JAB
						
stset time, id(StudyId) failure(event) origin(time 0)
stsplit d2, at(0) after(vd2)



					   *SPLIT TIME AFTER PRIMARY INFECTION (for positive cohort only)*
stsplit Primary_inf if(cohort_final==1), at(0 90 273 454) after(start_date_posC_2)
drop if Primary_inf==0	

*Label the pseudo observation (first line) and observations (second line) where there is no previous primary infection with value Primary_inf=0.
replace Primary_inf = Primary_inf + 1
replace Primary_inf=0 if Primary_inf==.




						*SPLIT IN BEFORE AND AFTER ENROLMENT*
*This is needed because some people enrolled after Dec 7th, and we need to tell state that follow-up time must not be used*
stsplit is_enrolled, at(0) after(Enrolment)

*Label as 0 the pseudo-observations that need not be used, 1 those that need to be used.
replace is_enrolled=is_enrolled+1

*Before this line, _st==1. Next line tells stata not to use time before enrolment.
replace _st=is_enrolled
*IMPORTANT: together with the above stset (with origin being time 0), this only takes into accounts follow-up time since ar.

drop if _st==0

***************************************************SPLITTING DONE*************************************************




					*LABELLING
					
recode Primary_inf 0=0 91=1 274=2 455=3  
label define prim_inf_lab 0 "No previous inf" 1 "3-9 mth" 2 "9-15_mth" 3 "15+_mth"
label values Primary_inf prim_inf_lab


replace d1=. if d2>=0


*Unvaccinated
generate x=0 if d1==-1 & Primary_inf==0
replace x=1 if d1==-1 & Primary_inf==1
replace x=2 if d1==-1 & Primary_inf==2
replace x=3 if d1==-1 & Primary_inf==3


*First dose follow up chunks
replace x=4 if d1==0 & Primary_inf==0
replace x=5 if d1==0 & Primary_inf==1
replace x=6 if d1==0 & Primary_inf==2
replace x=7 if d1==0 & Primary_inf==3


*Second dose follow up chunks
replace x=8 if d2==0 & Primary_inf==0
replace x=9 if d2==0 & Primary_inf==1
replace x=10 if d2==0 & Primary_inf==2
replace x=11 if d2==0 & Primary_inf==3



label define xlab 0 "Unvaccinated, naive" 1 "Unvaccinated, 3-9"  2 "Unvaccinated, 9-15" 3  "Unvaccinated, >15" ///
4 "d1, naive" 5 "d1, 3-9"  6 "d1, 9-15" 7  "d1, >15" ///
8 "d2, naive" 9 "d2, 3-9"  10 "d2, 9-15" 11  "d2, >15" ///

label values x xlab

*Drop copies of pseudo-observations of people who got an infection (before the infection)
drop if after_primary ==1 & ibv ==1 & cohort_final ==1 & x==0 & Primary_inf ==0


			
			
			*VISUALISATION

*Expousre is used to show the information when visualising splitted observations, and generate below tables. Needed for Poisson regression too.
generate exposure=_t - _t0

preserve
collapse (sum) exposure (sum) event,   by(x)
gen rate=round((10000*event/exposure),0.01)
lis *, sep(0)
restore

preserve
collapse (sum) exposure (sum) event
gen rate=round((10000*event/exposure),0.01)
lis *, sep(0)
restore


/*
*need to drop only x=0, if x not zero and <4 they need to be counted again.
drop if after_primary ==1 & x==0
levelsof x, local(levels) 
foreach i of local levels {

	display ""
	display ""
	display ""
	display ""
	display ""
	display ""

	display `i'
	distinct StudyId if x== `i' & _st==1
	
	display ""
	display ""
	display ""
	display ""
	display ""
	display ""
	
	
	}
*/


*Generated only for visualisation
gen PrimaryPCR= PrimaryPCRdate - start_time
gen LastPCRneg= LastPCRneg_date - start_time
gen vax_interval = vax_date2 - vax_date1


order StudyId Date_Enrolled Enrolment ar  whole_time start_date_posC PrimaryPCRdate ReinfectionPCRdate  LastPCRneg vax_date1 vax_date2 _t0 _t time x d1  d2 exposure _d  Primary_inf is_enrolled _st
gen interval = vaccine_date2 - start_date_posC
save splitted_comb, replace

exit
stcox i.x, vce(cluster Trust_Code) nolog
stcox i.x i.region i.AGEGR2 i.GENDER i.ETHNIC_GR ib5.WORK_EXPOSURE_FREQUENCY i.occ_set_cat, vce(cluster Trust_Code) nolog
stcox i.x i.region i.AGEGR2 i.GENDER i.ETHNIC_GR ib5.WORK_EXPOSURE_FREQUENCY i.occ_set_cat, shared(Trust_ode) forceshared nolog
