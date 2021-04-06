clear all
cd "C:\Users\jbueren\Google Drive\ABC\data"

set maxvar 20000
set more off
		/* Comments on wave use:
		 - I use wave 3 to 12
		 */

use hhid pn hhidpn rabmonth rabyear ragender raeduc racohbyr radyear using "rndhrs_p.dta"

local var_R iwstat iwendm iwendy shlt ///
 walkra dressa batha eata beda toilta /// ADLs
 mapa phonea moneya medsa shopa mealsa /// IADLs
 walksa chaira stoopa lifta clim1 dimea armsa pusha /// Other Functional Limitations
 hibpe diabe cancre lung hearte stroke psyche arthre /// chronic conditions
 imrc dlrc ser7 bwc20 scis cact pres vp mo dy yr dw /// cognitive impairment
 bmi drinkd drinkn smoken smokev /// 
 homcar nhmliv oopmd mstat govmd ///
 ipena isret igxfr

/*
- shlt: self reported health (1: Excellent, 5: Poor)

	ADLs (6):
- walkra: reports some difficulty walking across room (1: Yes, 0: No)
- dressa: reports some difficulty dressing (1: Yes, 0: No)
- batha: reports some difficulty bathing or showering (1: Yes, 0: No)
- eata: reports some difficulty eating (1: Yes, 0: No)
- beda: reports some difficulty getting in/out of bed (1: Yes, 0: No)
- toilta: reports some difficulty using toilet (1: Yes, 0: No)

	IADLS (6):
- mapa: reports some difficulty using a map (1: Yes, 0: No)
- phonea: reports some difficulty using a telephone (1: Yes, 0: No)
- moneya: reports some difficulty managing money (1: Yes, 0: No)
- medsa: reports some difficulty taking medications (1: Yes, 0: No)
- shopa: reports some difficulty (1: Yes, 0: No)
- mealsa: reports some difficulty preparing hot meals (1: Yes, 0: No)

	CHRONIC CONDITIONS (8)
-hibpe: ever had high blood pressure (I do not consider it: low mortality predictor)
-diabe: ever had high diabetes
-cancre: ever had cancer
-lung: ever had lung disease
-hearte: ever had heart problems
-stroke: ever had a stroke
-psyche: ever had psych problems
-arthre: ever had arthritis (I do not consider it: low mortality predictor)

*/

foreach var in `var_R'  {
display ("`var'")
if "`var'"=="imrc" | "`var'"=="dlrc" | "`var'"=="ser7"| "`var'"=="bwc20"| "`var'"=="scis"| "`var'"=="cact"| "`var'"=="pres"| "`var'"=="vp"| "`var'"=="mo"| "`var'"=="dy"| "`var'"=="yr"| "`var'"=="dw" {
quietly merge 1:1 hhidpn using "rndhrs_p.dta",keepusing(r3`var' r4`var' r5`var' r6`var' r7`var' r8`var' r9`var' r10`var' r11`var')
gen r12`var'=.
}
else{
quietly merge 1:1 hhidpn using "rndhrs_p.dta",keepusing(r3`var' r4`var' r5`var' r6`var' r7`var' r8`var' r9`var' r10`var' r11`var' r12`var')
}
drop _merge
}

local var_H hhres itot icap iothr
foreach var in `var_H'{
merge 1:1 hhidpn using "rndhrs_p.dta",keepusing(h3`var' h4`var' h5`var' h6`var' h7`var' h8`var' h9`var' h10`var' h11`var' h12`var')
drop _merge
}

local var_H atotb
foreach var in `var_H'{
merge 1:1 hhidpn using "rndhrs_p.dta",keepusing( h4`var' h5`var' h6`var' h7`var' h8`var' h9`var' h10`var' h11`var' h12`var')
drop _merge
}
gen h3atotb=.

local var_H hhres itot icap iothr atotb

*First wave
local f_w 3
*Last wave 
local l_w 12
								/* Age at interview if alive */
forval i=`f_w'/`l_w'{					
gen r`i'int_age=r`i'iwendy-rabyear if r`i'iwstat==1 
replace r`i'int_age=r`i'int_age+(r`i'iwendm-rabmonth)*1/12 if r`i'iwstat==1
}

drop rabmonth rabyear

								/*Reshape data*/					
foreach var in  `var_R' int_age {					
forval i=`f_w'/`l_w'{
rename r`i'`var' `var'`i'
}
}
foreach var in `var_H'{
forval i=`f_w'/`l_w'{
rename h`i'`var' `var'`i'
}
}


reshape long `var_R' int_age `var_H' ,i(hhidpn) j(wave)

drop iwendm

							/* Sample selection */
							
*Misreported indvidual (resurrection otherwise: "no estaba muerta estaba de parranda")
replace iwstat=4 if hhidpn==202147020 & wave==10

*drop after dead or dropped from sample
drop if iwstat==6 | iwstat==7 | iwstat==0

tsset hhidpn wave

*Age increases by two years between waves
by hhidpn: replace int_age=l1.int_age+2 if _n>1 
replace int_age=round(int_age)
drop if int_age<60 | int_age>98
replace int_age=int_age-1 if mod(int_age,2)==1

*Drop all first observations which are not interviews
forval w=`f_w'/`l_w'{
by hhidpn: drop if iwstat[1]!=1
}

hist int_age,bin(20)

recode walkra-arthre (.=-9)(.d=-9)(.r=-9)(.s=-9)(.x=-9)(.m=-9)(.t=-9)(.z=-9) if iwstat==1 | iwstat==4
recode walkra-arthre (.=-1) if iwstat==5
recode lung (3/5=1) //previous condition

gen generation=(int_age-60)/2+1

*drop individuals who are in the sample once but are not interviewed when alive
by hhidpn:drop if _N==1 & iwstat!=1

*Generate variable for the first and last "generation" the individual was observed
by hhidpn: gen first_age=generation[1]
by hhidpn: gen last_age=generation[_N]
order generation,after(hhidpn)

tsset hhidpn generation
tsfill,full

recode walkra-arthre (.=-1)

replace itot=itot-icap-iothr //I remove capital income and other income from household income to get a measure of permanent income
drop icap iothr

foreach var in  oopmd itot atotb ipena isret igxfr {
replace `var'=`var'*1.09986504723347 if iwendy== 1995 & `var'!=.
replace `var'=`var'*1.06955380577428 if iwendy== 1996 & `var'!=.
replace `var'=`var'*1.03887826641173 if iwendy== 1997 & `var'!=.
replace `var'=`var'*1.01557632398754 if iwendy== 1998 & `var'!=.
replace `var'=`var'*1 if iwendy== 1999 & `var'!=.
replace `var'=`var'*0.978391356542617 if iwendy== 2000 & `var'!=.
replace `var'=`var'*0.9465737514518 if iwendy== 2001 & `var'!=.
replace `var'=`var'*0.920383963862225 if iwendy== 2002 & `var'!=.
replace `var'=`var'*0.906058921623124 if iwendy== 2003 & `var'!=.
replace `var'=`var'*0.885869565217391 if iwendy== 2004 & `var'!=.
replace `var'=`var'*0.862890418210693 if iwendy== 2005 & `var'!=.
replace `var'=`var'*0.834613415258576 if iwendy== 2006 & `var'!=.
replace `var'=`var'*0.808531746031746 if iwendy== 2007 & `var'!=.
replace `var'=`var'*0.786140772250678 if iwendy== 2008 & `var'!=.
replace `var'=`var'*0.757072590720984 if iwendy== 2009 & `var'!=.
replace `var'=`var'*0.759775703025585 if iwendy== 2010 & `var'!=.
replace `var'=`var'*0.74751439997065 if iwendy== 2011 & `var'!=.
replace `var'=`var'*0.724640902644717 if iwendy== 2012 & `var'!=.
replace `var'=`var'*0.709948866259571 if iwendy== 2013 & `var'!=.
replace `var'=`var'*0.699699944624974 if iwendy== 2014 & `var'!=.
replace `var'=`var'*0.68853068396864 if iwendy== 2015 & `var'!=.
}


save data.dta,replace
*/

			/* Input data for HMC estimation */
			
clear all
use data.dta

by hhidpn: egen educ=mean(raeduc)
drop raeduc						
by hhidpn: egen gender=mean(ragender)	
drop ragender

preserve
drop shlt
keep if educ!=.
sort hhidpn generation
keep walkra-mealsa // only adls, iadls and chronic conditions
*drop arthre hibpe //low mortality predictor
export delimited using data_all.csv,replace nolab novar

restore
preserve
drop shlt
keep if educ!=.
keep hhidpn first_age last_age
collapse first_age last_age,by(hhidpn)
sort hhidpn
keep first_age last_age
export delimited using ages_all.csv,replace nolab novar

restore
preserve
drop shlt
keep if educ!=.
keep hhidpn gender
collapse gender,by(hhidpn)
sort hhidpn
keep gender
export delimited using gender_all.csv,replace nolab novar

restore
preserve
drop shlt
keep if educ!=.
keep hhidpn educ
collapse educ,by(hhidpn)
sort hhidpn
keep educ
export delimited using educ_all.csv,replace nolab novar

restore
preserve

restore
preserve
drop shlt
keep if educ!=.
sort hhidpn generation
gen loopmd=log(oopmd)
keep loopmd // only adls, iadls and chronic conditions
*drop arthre hibpe //low mortality predictor
recode loopmd (.=-1)
export delimited using loopmd.csv,replace nolab novar





