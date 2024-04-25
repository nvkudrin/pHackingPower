use "Data/MM data.dta" , clear
drop if mu==.
gen sign = 1 if mu>=0
replace sign = -1 if mu<0
gen ts = t*sign
drop if method!="RCT"
keep ts 
export delimited ts using "Brodeur_data.csv", replace
