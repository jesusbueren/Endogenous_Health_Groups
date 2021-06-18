			/* Input data for HMC estimation */
			
clear all
use data_ABC.dta

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
