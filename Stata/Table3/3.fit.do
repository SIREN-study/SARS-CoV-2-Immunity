use after_splitting, clear

**************************************
*VE for infected < 1 and > 1 year ago.
**************************************


*****WITHOUT PRIMARY INFECTION******
*Unvaccinated
generate x=0 if d1==-1 & Primary_inf ==0

*Chunks of follow-up (pseudo-observations) between first and second dose.
*First dose
replace x=1 if d1==0 & Primary_inf ==0
replace x=2 if d1>=20 &  d1<=79 & Primary_inf ==0

*Two doses
replace x= 3 if d2==0     & Primary_inf ==0  
replace x= 4 if d2==13    & Primary_inf ==0
replace x= 5 if d2==73    & Primary_inf ==0
replace x= 6 if (d2==133 | d2==163)   & Primary_inf ==0
replace x= 7 if d2==193   & Primary_inf ==0



*****PRIMARY INFECTION <1yr******

replace x=8 if d1==-1 & Primary_inf ==1

*Chunks of follow-up (pseudo-observations) between first and second dose.
*First dose
replace x=9 if d1==0 & Primary_inf ==1
replace x=10 if d1>=20 &  d1<=79 & Primary_inf ==1

*Two doses
replace x= 11 if d2==0     & Primary_inf ==1  
replace x= 12 if d2==13    & Primary_inf ==1
replace x= 13 if d2==73    & Primary_inf ==1
replace x= 14 if (d2==133 | d2==163)   & Primary_inf ==1

replace x= 15 if d2==193   & Primary_inf ==1


*****PRIMARY INFECTION >1yr******

replace x=16 if d1==-1 & Primary_inf ==2

*Chunks of follow-up (pseudo-observations) between first and second dose.
*First dose
replace x=17 if d1==0 & Primary_inf ==2
replace x=18 if d1>=20 &  d1<=79 & Primary_inf ==2
*Two doses
replace x= 19 if d2==0     & Primary_inf ==2  
replace x= 20 if d2==13    & Primary_inf ==2
replace x= 21 if d2==73    & Primary_inf ==2
replace x= 22 if (d2==133 | d2==163)   & Primary_inf ==2
replace x= 23 if d2==193   & Primary_inf ==2



*Labelling the values of x
label define xlab2 0 "Unvaccinated naive" ///
1 "d1 0-20 naive"  2 "d1 21+ naive" ///
3 "d2 0-13 naive" 4 "d2 14-73 naive" 5 "d2 73-133 naive" 6 "d2 133-193 naive" 7 "d2 193+ naive"  ///
8 "Unvaccinated inf<1yr" ///
9 "d1 0-20 inf<1yr" 10 "d1 21+ inf<1yr" /// 
11 "d2 0-13 inf<1yr" 12 "d2 14-73 inf<1yr" 13 "d2 73-133 inf<1yr" 14 "d2 133-193 inf<1yr" 15 "d2 193+ inf<1yr"   ///
16 "Unvaccinated inf>1yr" ///
17 "d1 0-20 inf>1yr" 18 "d1 21+ inf>1yr" ///
19 "d2 0-13 inf>1yr" 20 "d2 14-73 inf>1yr" 21 "d2 73-133 inf>1yr" 22 "d2 133-193 inf>1yr"  23 "d2 193+ inf>1yr"  
label values x xlab2

drop if _st==0
drop if vaccine_name2==2

save for_table, replace

*Exclude x=17 (because it does not converge) and AstraZeneca recipients
drop if  x==17


stcox i.x i.GENDER i.ETHNIC_GR, vce(cluster Trust_Code) nolog strata(AGEGR2 region WORK_EXPOSURE_FREQUENCY occ_set_cat) efron
