version 17
set more off
cap log close
clear all
set linesize 80

cd ""

local c_date = c(current_date)
local date = subinstr("`c_date'", " ", "", .)

log using "SDOH\Logs\OS_sdoh_analysis_newADIs_`date'.log", replace


***************************************************
* OVeR-Sepsis Social Determinants of Health
* Date Created: 2023 March 27  		
* Last updated: 2024 Sept 2  
* Author: Sarah Seelye
***************************************************

*------------------------------------
* Making Analytic Dataset
*------------------------------------
/*
	* open OS VAPD
	use "OverSepsis_VAPD20172021_20221115.dta", clear
	count // 16,908,990

	* correct hospital day
	order unique_hosp_count_id
	tab hospital_day if hospital_day==1 //n=2,580,379
	bysort unique_hosp_count_id (datevalue): gen n = _n 
	tab n if n==1 //n=2,580,390

	* merge with SDOH variable data set
	merge m:1 unique_hosp_count_id using os_vapd20172021_sdoh_20220307
	drop _merge

	count //16,908,990

	* keep hospitalization-level dataset 
	keep if n==1
	count //n=2,580,390

	* drop hospital_day variable 
	drop hospital_day n datevalue

	* organize dataset 
	order patientsid new_admitdate3 new_dischargedate3 admityear, after(patienticn)
	order deathdaysafteredis-mort90_edis, after(inhosp_mort)

	count //2,580,390

	* merge in new CAN scores
	merge 1:1 unique_hosp_count_id using SDOH\Data\os_can90d_2yrto2wk_admit, keepusing(can90d_2yrto2wk_preadmit) nogen
	count //2,580,390
	order can90d_2yrto2wk_preadmit, after(can90d_upto_disch)
	
	* drop hospitalizations with missing sta6as
	count if inlist(sta6a, "*Missing*", "*Unknown at this time*")
	drop if inlist(sta6a, "*Missing*", "*Unknown at this time*") //n=887
	count //2,579,503
	
	* drop hospitals with <300 hospitalizations 
	bysort sta6a: gen hosp_count = _N
	sum hosp_count, detail

	bysort sta6a: gen num_hosps = _n
	tab num_hosps if num_hosps==1 //n=181
		
	list sta6a hosp_count if hosp_count<300 & num_hosps==1
	
	drop if hosp_count<300 //dropped 703
	count //2,578,800
	
	tab num_hosps if num_hosps==1 //n=131

	drop hosp_count num_hosps

	* only keep patients with septic shock or severe sepsis
	* use over-sepsis indicator of sepsis 
	recode os_sepsis (.=0)
	tab os_sepsis admityear, m co

	keep if os_sepsis==1 //2,433,895
	count //144,905

	* drop hospitals with <15 sepsis hospitalizations
	bysort sta6a: gen sep_hosp_count = _N
	sum sep_hosp_count, detail
	list sta6a sep_hosp_count if sep_hosp_count<15
	drop if sep_hosp_count<15 //16 dropped
	count //144,889
	
	bysort sta6a: gen num_hosps = _n
	tab num_hosps if num_hosps==1 //n=129
	
	drop num_hosps
	
	tab os_sepsis admityear, m co
	
* Create derivation and validation datasets
	* Derivation dataset: 75% random sample for years 2017-2020
	* Validation I: 25% random sample for years 2017-2020
	* Validation II: year 2021

	gen admityear_20172020 = inrange(admityear, 2017, 2020)
	tab admityear admityear_20172020

	set seed 1589756
	gen random = runiform()
	gsort -admityear_20172020 random
	sum admityear_20172020 if admityear_20172020==1
	return list
	di ceil(0.75*r(N)) //83909

	gen derivation = admityear_20172020 & _n<=83909
	gen validation1 = admityear_20172020 & _n>83909 & _n<=111878
	gen validation2 = admityear==2021

	di (111878-83909)/111878

	tab derivation validation1
	tab derivation validation2
	tab validation1 validation2_

	drop admityear_20172020 random
	
	* save new dataset 
	*save SDOH\Data\os_sdoh_20220607, replace 
*/

use SDOH\Data\os_sdoh_20220607, clear
count //144889

*--------------------------------
* Add updated DODs for cohort
*--------------------------------

* import new dods 		
import delimited "sdoh_dod_20240731.csv", clear 
count
gen dod_20240731 = date(dod, "YMD")
format dod_20240731 %td	
drop dod sdoh_cohort
sort patienticn

* check for duplicates
duplicates report patienticn
duplicates tag patienticn, gen(dup)

* keep the earliest dod on record 
bysort patienticn (dod_20240731): gen num = _n
keep if num==1
duplicates report patienticn //118624

* merge with SDOH dataset 
merge 1:m patienticn using SDOH\Data\os_sdoh_20220607

* check that the older and new dod versions are consistent 
count if dod!=dod_20240731 & !missing(dod)

* replace missing values of new dod variable with dod from older version
replace dod_20240731=dod if dod_20240731==. & dod!=.

* drop older dod variable 
drop dod 
rename dod_20240731 dod
order dod, after(sta6a)

drop dup num _merge
		
*--------------
* New ADIs
*--------------

* drop original ADIs and import new ADIs
drop adi_natrank_num adi_staternk

merge 1:1 unique_hosp_count_id using SDOH\Data\os_adi_uniq_hosp_20230717
drop if _merge==2
rename ADI_quarter adi_quarter 
rename year adi_year
rename fips_geoid fips_geoid_block
rename admit_quarter admitquarter 
order admitquarter, after(admityear)
drop admit_year quarter lat lon qtr_yr* diff_qtr diff_yr _merge
order adi_quarter, after(adi_year)
rename score pssg_match_score 

* investigate PSSG scores
sum pssg_match_score 
count if pssg_match_score<70

* per PSSG documentation, recode ADIs missing if match score<70 - deemed 
* too low to be acceptable for use 
replace adi_natrank = . if pssg_match_score<70
replace adi_staternk = . if pssg_match_score<70

rename adi_staternk adi_staterank


*-----------
* SDIs
*-----------

* create a tract-level fips code 
gen str11 fips_geoid_tract = substr(fips_geoid_block, 1, 11)
destring fips_geoid_tract, gen(censustract_fips)
format censustract_fips %14.0g

* create a year variable for merging with SDI
	* note, the last year of SDI is 2019
clonevar year = admityear
replace year = 2019 if admityear>2019
tab admityear year

* merge 2013-2017 SDI with 2017 admityear 
merge m:1 censustract_fips year using SDOH\Data\SDI\rgcsdi_2013_2017_censustract 
drop if _merge==2

drop censustract_population povertylt100_fpl_score-nonemployed_score pct_poverty_lt100-_merge

rename sdi sdi_2017
rename sdi_score sdi_score_2017

* merge 2014-2018 SDI with 2018 admityear 
merge m:1 censustract_fips year using SDOH\Data\SDI\rgcsdi_2014_2018_censustract 
drop if _merge==2

drop censustract_population povertylt100_fpl_score-nonemployed_score pct_poverty_lt100-_merge

rename sdi sdi_2018
rename sdi_score sdi_score_2018

* merge 2015-2019 SDI with 2019 admityear 
merge m:1 censustract_fips year using SDOH\Data\SDI\rgcsdi_2015_2019_censustract 
drop if _merge==2

drop censustract_population povertylt100_fpl_score-nonemployed_score pct_poverty_lt100-_merge

rename sdi sdi_2019
rename sdi_score sdi_score_2019

* create single sdi variables 
gen sdi = .
replace sdi = sdi_2017 if year==2017
replace sdi = sdi_2018 if year==2018
replace sdi = sdi_2019 if year==2019

gen sdi_score = .
replace sdi_score = sdi_score_2017 if year==2017
replace sdi_score = sdi_score_2018 if year==2018
replace sdi_score = sdi_score_2019 if year==2019

bysort admityear: tab sdi_score, m // >25% missing for 2020 & 2021; ~3% missing for 2017-2019
corr sdi_score adi_natrank 

*-----------
* SVIs
*-----------

* create year variable for merging with svi
	* note, SVIs are created in two-year periods;
	* the final SVI dataset is 2020
	* 2017-2019 merges to 2018 SVI; 2020-2021 merges to 2020 SVI 
clonevar year_svi = admityear
replace year_svi = 2018 if admityear<=2019
replace year_svi = 2020 if admityear>=2020
tab admityear year_svi

* merge 2018 svi with 2017-2019 admityear 
merge m:1 censustract_fips year_svi using SDOH\Data\svi\svi_2018_us 
drop if _merge==2

order year_svi state rpl*, after(censustract_fips)

drop st-_merge rpl_theme2-rpl_theme4

rename rpl_theme1 rpl_theme1_2018
rename rpl_themes rpl_themes_2018
rename state state_2018

* merge 2020 svi with 2020-2021 admityear 
merge m:1 censustract_fips year_svi using SDOH\Data\svi\svi_2020_us 
drop if _merge==2

order state rpl*, after(state_2018)

drop st-_merge rpl_theme2-rpl_theme4

rename rpl_theme1 rpl_theme1_2020
rename rpl_themes rpl_themes_2020
rename state state_2020

* create single svi variables 
gen svi = .
replace svi = rpl_themes_2018 if year_svi==2018
replace svi = rpl_themes_2020 if year_svi==2020
replace svi = . if svi==-999

gen svi_ses = .
replace svi_ses = rpl_theme1_2018 if year_svi==2018
replace svi_ses = rpl_theme1_2020 if year_svi==2020
replace svi_ses = . if svi_ses==-999

gen state = "."
replace state = state_2018 if year_svi==2018
replace state = state_2020 if year_svi==2020

replace state = strupper(state)

drop rpl_themes* rpl_theme1* state_*

gen svimiss = missing(svi)
version 16: table admityear, c(count unique count svi mean svi mean svimiss sum svimiss) 
	// ~ 3% missing

corr svi adi_natrank 
corr svi_ses adi_natrank


*----------------------------------------------
* Drop sensitive comorbidities from dataset 
*----------------------------------------------
	* alcohol abuse
	* drug abuse 
	* AIDS/HIV

drop etoh drug ah  


*------------------
* New Variables 
*------------------

* Covariates

foreach var in 	acute_respiratory_hosp aod_lactate aod_heme ///
				pressor_fromed_72hr aod_kidney aod_liver { 
		replace `var'=0 if `var'==. 
}

gen aod_count = acute_respiratory_hosp + aod_lactate + aod_heme + ///
				pressor_fromed_72hr + aod_kidney + aod_liver
tab aod_count

sum age
gen age_65plus = age>65

gen female = gender=="F"
gen male = gender=="M"
	
gen black = race=="BLACK OR AFRICAN AMERICAN"
gen white = race=="WHITE" | race=="WHITE NOT OF HISP ORIG"
gen other = 0
replace other=1 if black==0 & white==0

tab race white
tab race black
tab race other

tab race, m
gen race_3cat = .
replace race_3cat = 1 if race=="WHITE" | race=="WHITE NOT OF HISP ORIG"
replace race_3cat = 2 if race=="BLACK OR AFRICAN AMERICAN"
replace race_3cat = 3 if race_3cat>2
lab def race_3cat 1 "White" 2 "Black" 3 "Other"
lab val race_3cat race_3cat
tab race race_3cat

foreach var in 	cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver  	 ///
				neuro renal htn cardic_arrhym valvular_d2 pulm_cir pvd  	 ///
				paralysis pud hypothyroid  lymphoma ra coag obesity wtloss ///
				fen anemia_cbl anemia_def   psychoses depression {

		recode `var'_prior540 (.=0)
		rename `var'_prior540 `var'
}	

tab urh 
replace urh="insular islands" if urh=="I"
replace urh="highly rural" if urh=="H"
replace urh="rural" if urh=="R"
replace urh="urban" if urh=="U"
gen rural = urh=="rural" | urh=="highly rural"
tab urh rural

mkspline age_sp = age, cubic nknots(5)

gen age_cat = .
replace age_cat = 1 if age>=18 & age<45
replace age_cat = 2 if age>=45 & age<65
replace age_cat = 3 if age>=65 & age<80
replace age_cat = 4 if age>=80

sort unique_hosp_count_id

* post-discharge mortality 
gen mort_days_postdc = dod-new_dischargedate3
sum mort_days_postdc, de

gen mort90_postdc = inrange(mort_days_postdc, 2, 90)
replace mort90_postdc = . if mort_days_postdc<2

gen dc_alive = mort_days_postdc>=2

tab mort90_postdc, m
tab mort90_postdc dc_alive, mi co

gen mort360_postdc = inrange(mort_days_postdc, 2, 360)
replace mort360_postdc = . if mort_days_postdc<2

tab mort360_postdc, m
tab mort360_postdc dc_alive, mi co

gen mort720_postdc = inrange(mort_days_postdc, 2, 720)
replace mort720_postdc = . if mort_days_postdc<2

tab mort720_postdc, m
tab mort720_postdc dc_alive, mi co

*-------------------------------
* Missingness and Imputation
*-------------------------------

* investigate missingness 
misstable sum cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp liver neuro 		///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 				///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen					///
			anemia_cbl anemia_def   psychoses depression 					///	
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver					///
			aod_heme acute_respiratory_hosp sirs_temp sirs_rr sirs_pulse sirs_wbc 	///
			male age can90d_2yrto2wk_preadmit frailty_score married 				///
			housing_instab_combo_prior540 rural adi_natrank 

count if missing(can90d_2yrto2wk_preadmit) & derivation==1 	//4.3% missing		
count if missing(can90d_2yrto2wk_preadmit) & validation1==1 //4.5% missing			
count if missing(can90d_2yrto2wk_preadmit) & validation2==1 //4.3% missing		

sum can90d_2yrto2wk_preadmit if derivation==1, de
sum can90d_2yrto2wk_preadmit if validation1==1, de
sum can90d_2yrto2wk_preadmit if validation2==1, de

* use single multivariable regression imputation to impute missing values for: 
		* CAN
		* married 
		* ADI_natrank 

* CAN 
local covar male age														///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 					///
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			aod_heme acute_respiratory_hosp 								///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			housing_instab_combo_prior540 frailty_score	rural 				
regress can90d_2yrto2wk_preadmit `covar' 		
predict mkg_can_impute

gen can_impute = can90d_2yrto2wk_preadmit
replace can_impute = mkg_can_impute if can90d_2yrto2wk_preadmit==. 

sum can90d_2yrto2wk_preadmit mkg_can_impute can_impute		
drop mkg_can_impute 

sum can_impute if derivation==1, de
sum can_impute if validation1==1, de
sum can_impute if validation2==1, de


* married 
local covar male age														///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 					///
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			aod_heme acute_respiratory_hosp 								///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			housing_instab_combo_prior540 frailty_score	rural 				
logit married `covar'		
predict mkg_married_impute

gen married_impute = married 
replace married_impute = round(mkg_married_impute) if married==.

drop mkg_married_impute

* ADI 
local covar male age														///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 					///
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			aod_heme acute_respiratory_hosp 								///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			housing_instab_combo_prior540 frailty_score	rural 				
regress adi_natrank `covar'		
predict mkg_adi_impute

gen adi_impute = adi_natrank 
replace adi_impute = mkg_adi_impute if adi_natrank==.

sum adi_natrank mkg_adi_impute adi_impute
drop mkg_adi_impute
			
*---------------------------
* Supplemental Table 4
*---------------------------

* hospitalizations 
tab admityear
tab derivation
tab validation1
tab validation2

* age 
sum age, de
sum age if derivation, de
sum age if validation1, de
sum age if validation2, de

* male 
tab male
tab male if derivation 
tab male if validation1
tab male if validation2 

* race 
tab race_3cat
tab race_3cat if derivation 
tab race_3cat if validation1
tab race_3cat if validation2 

* sirs 
tab sirs_wbc
tab sirs_temp
tab sirs_pulse
tab sirs_rr

tab sirs_wbc if derivation 
tab sirs_temp if derivation 
tab sirs_pulse if derivation 
tab sirs_rr if derivation 

tab sirs_wbc if validation1 
tab sirs_temp if validation1
tab sirs_pulse if validation1 
tab sirs_rr if validation1 

tab sirs_wbc if validation2 
tab sirs_temp if validation2
tab sirs_pulse if validation2 
tab sirs_rr if validation2

* aod 
tab aod_kidney
tab aod_lactate
tab aod_heme
tab aod_liver
tab acute_respiratory_hosp
tab pressor_fromed_72hr

tab aod_kidney if derivation 
tab aod_lactate if derivation 
tab aod_heme if derivation 
tab aod_liver if derivation 
tab acute_respiratory_hosp if derivation 
tab pressor_fromed_72hr if derivation 

tab aod_kidney if validation1 
tab aod_lactate if validation1
tab aod_heme if validation1 
tab aod_liver if validation1 
tab acute_respiratory_hosp if validation1 
tab pressor_fromed_72hr if validation1 

tab aod_kidney if validation2 
tab aod_lactate if validation2
tab aod_heme if validation2 
tab aod_liver if validation2
tab acute_respiratory_hosp if validation2
tab pressor_fromed_72hr if validation2

* chronic conditions
tab renal
tab liver
tab chf
tab pulm

tab renal if derivation 
tab liver if derivation 
tab chf if derivation 
tab pulm if derivation 

tab renal if validation1 
tab liver if validation1
tab chf if validation1 
tab pulm if validation1 

tab renal if validation2 
tab liver if validation2
tab chf if validation2 
tab pulm if validation2

* social determinants of health
sum can_impute, de 
sum frailty_score, de
tab married_impute, m
tab housing_instab_combo_prior540
tab rural
sum adi_impute, de

sum can_impute if derivation, de 
sum frailty_score if derivation, de
tab married_impute if derivation, m
tab housing_instab_combo_prior540 if derivation 
tab rural if derivation 
sum adi_impute if derivation, de

sum can_impute if validation1, de 
sum frailty_score if validation1, de
tab married_impute if validation1, m
tab housing_instab_combo_prior540 if validation1 
tab rural if validation1 
sum adi_impute if validation1, de

sum can_impute if validation2, de 
sum frailty_score if validation2, de
tab married_impute if validation2, m
tab housing_instab_combo_prior540 if validation2
tab rural if validation2
sum adi_impute if validation2, de

*--------------------------------------------
* Mortality models including all variables 
*--------------------------------------------

** 30-day mortality **

* predicted 30-day mortality, full sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort30_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute || sta6a: 
predict pr_mort30_admit_all, mu	
margins, nose

predict re_mort30_all if e(sample), reffects reses(rese_mort30_all)
predict xb_mort30_all if e(sample), xb
	replace xb_mort30_all=. if derivation!=1

egen meanxb_mort30_all = mean(xb_mort30_all) if e(sample)
gen hosplogodd_mort30_all = meanxb_mort30_all + re_mort30_all
gen invlog_mort30_all = invlogit(hosplogodd_mort30_all)

sum pr_mort30_admit_all invlog_mort30_all, de


* predicted 30-day mortality, derivation sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort30_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if derivation==1 || sta6a: 
predict pr_mort30_admit_deriv if e(sample), mu	
predict re_mort30_deriv if e(sample), reffects reses(rese_mort30_deriv)
predict xb_mort30_deriv if e(sample), xb
	replace xb_mort30_deriv=. if derivation!=1

egen meanxb_mort30_deriv = mean(xb_mort30_deriv) if e(sample)
gen hosplogodd_mort30_deriv = meanxb_mort30_deriv + re_mort30_deriv
gen invlog_mort30_deriv = invlogit(hosplogodd_mort30_deriv)

sum pr_mort30_admit_deriv invlog_mort30_deriv, de


* predicted 30-day mortality, validation 1 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort30_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation1==1 || sta6a: 
predict pr_mort30_admit_valid1 if e(sample), mu	
predict re_mort30_valid1 if e(sample), reffects reses(rese_mort30_valid1)
predict xb_mort30_valid1 if e(sample), xb
	replace xb_mort30_valid1=. if validation1!=1

egen meanxb_mort30_valid1 = mean(xb_mort30_valid1) if e(sample)
gen hosplogodd_mort30_valid1 = meanxb_mort30_valid1 + re_mort30_valid1
gen invlog_mort30_valid1 = invlogit(hosplogodd_mort30_valid1)

sum pr_mort30_admit_valid1 invlog_mort30_valid1, de


* predicted 30-day mortality, validation 2 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort30_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation2==1 || sta6a: 
predict pr_mort30_admit_valid2 if e(sample), mu	
predict re_mort30_valid2 if e(sample), reffects reses(rese_mort30_valid2)
predict xb_mort30_valid2 if e(sample), xb
	replace xb_mort30_valid2=. if validation2!=1

egen meanxb_mort30_valid2 = mean(xb_mort30_valid2) if e(sample)
gen hosplogodd_mort30_valid2 = meanxb_mort30_valid2 + re_mort30_valid2
gen invlog_mort30_valid2 = invlogit(hosplogodd_mort30_valid2)

sum pr_mort30_admit_valid2 invlog_mort30_valid2, de


** 90-day mortality **

* predicted 90-day mortality, full sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute || sta6a: 
predict pr_mort90_admit_all, mu	
predict re_mort90_all if e(sample), reffects reses(rese_mort90_all)
predict xb_mort90_all if e(sample), xb

egen meanxb_mort90_all = mean(xb_mort90_all) if e(sample)
gen hosplogodd_mort90_all = meanxb_mort90_all + re_mort90_all
gen invlog_mort90_all = invlogit(hosplogodd_mort90_all)

sum pr_mort90_admit_all invlog_mort90_all, de


* predicted 90-day mortality, derivation sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if derivation==1 || sta6a: 
predict pr_mort90_admit_deriv if e(sample), mu	
predict re_mort90_deriv if e(sample), reffects reses(rese_mort90_deriv)
predict xb_mort90_deriv if e(sample), xb
	replace xb_mort90_deriv=. if derivation!=1

egen meanxb_mort90_deriv = mean(xb_mort90_deriv) if e(sample)
gen hosplogodd_mort90_deriv = meanxb_mort90_deriv + re_mort90_deriv
gen invlog_mort90_deriv = invlogit(hosplogodd_mort90_deriv)

sum pr_mort90_admit_deriv invlog_mort90_deriv, de


* predicted 90-day mortality, validation 1 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation1==1 || sta6a: 
predict pr_mort90_admit_valid1 if e(sample), mu	
predict re_mort90_valid1 if e(sample), reffects reses(rese_mort90_valid1)
predict xb_mort90_valid1 if e(sample), xb
	replace xb_mort90_valid1=. if validation1!=1

egen meanxb_mort90_valid1 = mean(xb_mort90_valid1) if e(sample)
gen hosplogodd_mort90_valid1 = meanxb_mort90_valid1 + re_mort90_valid1
gen invlog_mort90_valid1 = invlogit(hosplogodd_mort90_valid1)

sum pr_mort90_admit_valid1 invlog_mort90_valid1, de


* predicted 90-day mortality, validation 2 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_admit male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation2==1 || sta6a: 
predict pr_mort90_admit_valid2 if e(sample), mu	
predict re_mort90_valid2 if e(sample), reffects reses(rese_mort90_valid2)
predict xb_mort90_valid2 if e(sample), xb
	replace xb_mort90_valid2=. if validation2!=1

egen meanxb_mort90_valid2 = mean(xb_mort90_valid2) if e(sample)
gen hosplogodd_mort90_valid2 = meanxb_mort90_valid2 + re_mort90_valid2
gen invlog_mort90_valid2 = invlogit(hosplogodd_mort90_valid2)

sum pr_mort90_admit_valid2 invlog_mort90_valid2, de


** inhospital mortality **

* predicted inhospital mortality, full sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit inhosp_mort male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute || sta6a: 
predict pr_mortinhosp_admit_all, mu	
predict re_mortinhosp_all if e(sample), reffects reses(rese_mortinhosp_all)
predict xb_mortinhosp_all if e(sample), xb

egen meanxb_mortinhosp_all = mean(xb_mortinhosp_all) if e(sample)
gen hosplogodd_mortinhosp_all = meanxb_mortinhosp_all + re_mortinhosp_all
gen invlog_mortinhosp_all = invlogit(hosplogodd_mortinhosp_all)

sum pr_mortinhosp_admit_all invlog_mortinhosp_all, de


* predicted inhospital mortality, derivation sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit inhosp_mort male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if derivation==1 || sta6a: 
predict pr_mortinhosp_admit_deriv if e(sample), mu	
predict re_mortinhosp_deriv if e(sample), reffects reses(rese_mortinhosp_deriv)
predict xb_mortinhosp_deriv if e(sample), xb
	replace xb_mortinhosp_deriv=. if derivation!=1

egen meanxb_mortinhosp_deriv = mean(xb_mortinhosp_deriv) if e(sample)
gen hosplogodd_mortinhosp_deriv = meanxb_mortinhosp_deriv + re_mortinhosp_deriv
gen invlog_mortinhosp_deriv = invlogit(hosplogodd_mortinhosp_deriv)

sum pr_mortinhosp_admit_deriv invlog_mortinhosp_deriv, de


* predicted inhospital mortality, validation 1 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit inhosp_mort male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation1==1 || sta6a: 
predict pr_mortinhosp_admit_valid1 if e(sample), mu	
predict re_mortinhosp_valid1 if e(sample), reffects reses(rese_mortinhosp_valid1)
predict xb_mortinhosp_valid1 if e(sample), xb
	replace xb_mortinhosp_valid1=. if validation1!=1

egen meanxb_mortinhosp_valid1 = mean(xb_mortinhosp_valid1) if e(sample)
gen hosplogodd_mortinhosp_valid1 = meanxb_mortinhosp_valid1 + re_mortinhosp_valid1
gen invlog_mortinhosp_valid1 = invlogit(hosplogodd_mortinhosp_valid1)

sum pr_mortinhosp_admit_valid1 invlog_mortinhosp_valid1, de


* predicted inhospital mortality, validation 2 sample 
local comorbid 															///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 	///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver				///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit inhosp_mort male age_sp* i.race_3cat `sirs' `aod' `comorbid' 	///
		can_impute frailty_score married_impute 			///
		housing_instab_combo_prior540 rural adi_impute if validation2==1 || sta6a: 
predict pr_mortinhosp_admit_valid2 if e(sample), mu	
predict re_mortinhosp_valid2 if e(sample), reffects reses(rese_mortinhosp_valid2)
predict xb_mortinhosp_valid2 if e(sample), xb
	replace xb_mortinhosp_valid2=. if validation2!=1

egen meanxb_mortinhosp_valid2 = mean(xb_mortinhosp_valid2) if e(sample)
gen hosplogodd_mortinhosp_valid2 = meanxb_mortinhosp_valid2 + re_mortinhosp_valid2
gen invlog_mortinhosp_valid2 = invlogit(hosplogodd_mortinhosp_valid2)

sum pr_mortinhosp_admit_valid2 invlog_mortinhosp_valid2, de


*--------------------------------------------------------------
* Case-mix variation among hospitalizations at VA hospitals 
*--------------------------------------------------------------

* categorize hospitals into quartiles using each hospital's median ADI score
bysort sta6a: egen adi_impute_median=median(adi_impute)
bysort sta6a: gen mkg_adi_hospitalquart = adi_impute_median if _n==1

sum mkg_adi_hospitalquart, de	

xtile adi_impute_quart = mkg_adi_hospitalquart, nq(4) //129 hospitals in quartiles
bysort adi_impute_quart: sum mkg_adi_hospitalquart, det

bysort sta6a: egen adi_hospitalquart = max(adi_impute_quart)
tab adi_hospitalquart

sort mkg_adi_hospitalquart
list sta6a mkg_adi_hospitalquart adi_hospitalquart if mkg_adi_hospitalquart!=.

	* double-checking quartiles
	preserve 
		collapse adi_impute_median, by(sta6a)
		xtile adi_impute_quart2 = adi_impute_median, nq(4)
		tab adi_impute_quart2
		bysort adi_impute_quart2: tab sta6a
		bysort adi_impute_quart2: sum adi_impute_median, det
		sort adi_impute_median
		list sta6a adi_impute_median
	restore

* count number of hospitals in each quartile
preserve 
	collapse (count) counthospsbyquart=unique_hosp_count_id  , by(sta6a adi_hospitalquart)
	sort adi_hospitalquart
	bysort adi_hospitalquart: tab sta6a
	list sta6a adi_hospitalquart counthospsbyquart
restore 

* Table *
preserve 

	collapse	(count) counthospsbyquart=unique_hosp_count_id  				///
				(p50) age_quartp50=age 											///
				(p25) age_quartp25=age 											///
				(p75) age_quartp75=age											///
				(sum) white (sum) black (sum) other								///
				(sum) male (sum) sirs_wbc (sum) sirs_temp (sum) sirs_pulse   	///
				(sum) sirs_rr (sum) aod_kidney (sum) aod_lactate 				///
				(sum) aod_heme (sum) aod_liver (sum) acute_respiratory_hosp		///
				(sum) pressor_fromed_72hr (sum) renal (sum) liver 				///
				(sum) chf (sum) pulm 											///
				(p50) can_quartp50=can_impute 					///
				(p25) can_quartp25=can_impute 					///
				(p75) can_quartp75=can_impute 					///
				(p50) frailty_quartp50=frailty_score							///
				(p25) frailty_quartp25=frailty_score							///
				(p75) frailty_quartp75=frailty_score							///
				(p50) adi_impute_quartp50=adi_impute								///
				(p25) adi_impute_quartp25=adi_impute								///
				(p75) adi_impute_quartp75=adi_impute								///
				(sum) married_impute (sum) housing_instab_combo_prior540 (sum) rural 	///
				(p50) mort30_quart50 = pr_mort30_admit_all						///
				(p25) mort30_quart25 = pr_mort30_admit_all						///
				(p75) mort30_quart75 = pr_mort30_admit_all						///
				(p50) mort90_quart50 = pr_mort90_admit_all						///
				(p25) mort90_quart25 = pr_mort90_admit_all						///
				(p75) mort90_quart75 = pr_mort90_admit_all						///
				(p50) mortinhosp_quart50 = pr_mortinhosp_admit_all				///
				(p25) mortinhosp_quart25 = pr_mortinhosp_admit_all				///
				(p75) mortinhosp_quart75 = pr_mortinhosp_admit_all				///
				(p50) mort30_quart50invlog = invlog_mort30_all						///
				(p25) mort30_quart25invlog = invlog_mort30_all						///
				(p75) mort30_quart75invlog = invlog_mort30_all						///
				(p50) mort90_quart50invlog = invlog_mort90_all						///
				(p25) mort90_quart25invlog = invlog_mort90_all						///
				(p75) mort90_quart75invlog = invlog_mort90_all						///
				(p50) mortinhosp_quart50invlog = invlog_mortinhosp_all				///
				(p25) mortinhosp_quart25invlog = invlog_mortinhosp_all			///
				(p75) mortinhosp_quart75invlog = invlog_mortinhosp_all				///
				, by(adi_hospitalquart)			
						
		gen pwhite = white/counthosps
		gen pblack = black/counthosps
		gen pother = other/counthosps
		gen pmale = male/counthosps
		gen psirs_wbc = sirs_wbc/counthosps 
		gen psirs_temp = sirs_temp/counthosps 
		gen psirs_pulse = sirs_pulse/counthosps  
		gen psirs_rr = sirs_rr/counthosps 
		gen paod_kidney = aod_kidney/counthosps 
		gen paod_lactate = aod_lactate/counthosps 	
		gen paod_heme = aod_heme/counthosps 
		gen paod_liver = aod_liver/counthosps 
		gen pacute_respiratory_hosp = acute_respiratory_hosp/counthosps	
		gen ppressor_fromed_72hr = pressor_fromed_72hr/counthosps 
		gen prenal = renal/counthosps 
		gen pliver = liver/counthosps
		gen pchf = chf/counthosps 
		gen ppulm = pulm/counthosps 
		gen pmarried_impute = married_impute/counthosps 
		gen phousing_instab = housing_instab_combo_prior540/counthosps
		gen prural = rural/counthosps 
		
		drop white black other male sirs_wbc sirs_temp sirs_pulse sirs_rr aod_kidney aod_lactate ///
			 aod_heme aod_liver acute_respiratory_hosp pressor_fromed_72hr renal liver 	///
			 chf pulm married_impute housing_instab_combo_prior540	rural 
			 
		list adi_hospitalquart age_quartp50 age_quartp25 age_quartp75
		
		list adi_hospitalquart pmale 
		list adi_hospitalquart pwhite
		list adi_hospitalquart pblack
		list adi_hospitalquart pother

		list adi_hospitalquart psirs_wbc 
		list adi_hospitalquart psirs_temp
		list adi_hospitalquart psirs_pulse
		list adi_hospitalquart psirs_rr

		list adi_hospitalquart paod_kidney
		list adi_hospitalquart paod_lactate
		list adi_hospitalquart paod_heme
		list adi_hospitalquart paod_liver
		list adi_hospitalquart pacute_respiratory_hosp
		list adi_hospitalquart ppressor_fromed_72hr

		list adi_hospitalquart prenal
		list adi_hospitalquart pliver
		list adi_hospitalquart pchf
		list adi_hospitalquart ppulm

		list adi_hospitalquart can_quartp50 can_quartp25 can_quartp75
		list adi_hospitalquart frailty_quartp50 frailty_quartp25 frailty_quartp75

		list adi_hospitalquart pmarried_impute
		list adi_hospitalquart phousing_instab
		list adi_hospitalquart prural
		list adi_hospitalquart adi_impute_quartp50 adi_impute_quartp25 adi_impute_quartp75

		list adi_hospitalquart mort30_quart50 mort30_quart25 mort30_quart75
		list adi_hospitalquart mort90_quart50 mort90_quart25 mort90_quart75
		list adi_hospitalquart mortinhosp_quart50 mortinhosp_quart25 mortinhosp_quart75
		
		list adi_hospitalquart mort30_quart50invlog mort30_quart25invlog mort30_quart75invlog
		list adi_hospitalquart mort90_quart50invlog mort90_quart25invlog mort90_quart75invlog
		list adi_hospitalquart mortinhosp_quart50invlog mortinhosp_quart25invlog mortinhosp_quart75invlog
		
restore 


* ADI hospital quartile test of trends
* asymptotic p-values
nptrend age, group(adi_hospitalquart) jterpstra 
nptrend male, group(adi_hospitalquart) carmitage
nptrend sirs_wbc, group(adi_hospitalquart) carmitage
nptrend sirs_temp, group(adi_hospitalquart) carmitage
nptrend sirs_pulse, group(adi_hospitalquart) carmitage
nptrend sirs_rr, group(adi_hospitalquart) carmitage
nptrend aod_kidney, group(adi_hospitalquart) carmitage
nptrend aod_lactate, group(adi_hospitalquart) carmitage
nptrend aod_heme, group(adi_hospitalquart) carmitage
nptrend aod_liver, group(adi_hospitalquart) carmitage
nptrend acute_respiratory_hosp, group(adi_hospitalquart) carmitage
nptrend pressor_fromed_72hr, group(adi_hospitalquart) carmitage
nptrend renal, group(adi_hospitalquart) carmitage
nptrend liver, group(adi_hospitalquart) carmitage
nptrend chf, group(adi_hospitalquart) carmitage
nptrend pulm, group(adi_hospitalquart) carmitage
nptrend can_impute, group(adi_hospitalquart) jterpstra 
nptrend frailty_score, group(adi_hospitalquart) jterpstra 
nptrend married_impute, group(adi_hospitalquart) carmitage
nptrend housing_instab_combo_prior540, group(adi_hospitalquart) carmitage
nptrend rural, group(adi_hospitalquart) carmitage
nptrend pr_mort30_admit_all, group(adi_hospitalquart) jterpstra 
nptrend pr_mort90_admit_all, group(adi_hospitalquart) jterpstra 
nptrend pr_mortinhosp_admit_all, group(adi_hospitalquart) jterpstra 

* exact p-values
nptrend age, group(adi_hospitalquart) jterpstra exact 
nptrend male, group(adi_hospitalquart) carmitage exact
nptrend black, group(adi_hospitalquart) carmitage exact
nptrend white, group(adi_hospitalquart) carmitage exact
nptrend other, group(adi_hospitalquart) carmitage exact
nptrend sirs_wbc, group(adi_hospitalquart) carmitage exact
nptrend sirs_temp, group(adi_hospitalquart) carmitage exact
nptrend sirs_pulse, group(adi_hospitalquart) carmitage exact
nptrend sirs_rr, group(adi_hospitalquart) carmitage exact
nptrend aod_kidney, group(adi_hospitalquart) carmitage exact
nptrend aod_lactate, group(adi_hospitalquart) carmitage exact
nptrend aod_heme, group(adi_hospitalquart) carmitage exact
nptrend aod_liver, group(adi_hospitalquart) carmitage exact
nptrend acute_respiratory_hosp, group(adi_hospitalquart) carmitage exact
nptrend pressor_fromed_72hr, group(adi_hospitalquart) carmitage exact
nptrend renal, group(adi_hospitalquart) carmitage exact
nptrend liver, group(adi_hospitalquart) carmitage exact
nptrend chf, group(adi_hospitalquart) carmitage exact
nptrend pulm, group(adi_hospitalquart) carmitage exact
nptrend can_impute, group(adi_hospitalquart) jterpstra exact
nptrend frailty_score, group(adi_hospitalquart) jterpstra exact
nptrend married_impute, group(adi_hospitalquart) carmitage exact
nptrend housing_instab_combo_prior540, group(adi_hospitalquart) carmitage exact
nptrend rural, group(adi_hospitalquart) carmitage exact
nptrend adi_impute, group(adi_hospitalquart) jterpstra exact
nptrend pr_mort30_admit_all, group(adi_hospitalquart) jterpstra exact
nptrend pr_mort90_admit_all, group(adi_hospitalquart) jterpstra exact
nptrend pr_mortinhosp_admit_all, group(adi_hospitalquart) jterpstra exact						

	
* drop all previous predicted mortality variables 
drop pr_mort*
				
********************************************************************************
* Table 2: Comparison of model performance in derivation and validation 1 & 2
********************************************************************************

*----------------------
* M0: Community Only
*----------------------

* DERIVATION *
melogit mort90_edis adi_impute if derivation==1 || sta6a: 
estimates store derivation_m0

predict pr_derivation_m0 if e(sample)
predict xb_derivation_m0 if e(sample), xb
brier mort90_edis pr_derivation_m0 

* standardized mortality ratio
total mort90_edis if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m0 if derivation==1 
	matrix list e(b)
	local prderivationdeathm0=e(b)[1,1]
	di `prderivationdeathm0'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm0'	
	
* excess deaths 
di `derivationdeath' - `prderivationdeathm0'	
di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m0

predict xb_validation1_m0 if validation1==1, xb 
	replace xb_validation1_m0=. if validation1!=1
predict pr_validation1_m0 if validation1==1	
predict re_validation1_m0 if validation1==1, reffects reses(rese_validation1_m0)
predict stdp_validation1_m0 if validation1==1, stdp 
	replace stdp_validation1_m0=. if validation1!=1

egen meanxb_validation1_m0 = mean(xb_validation1_m0) if validation1==1
gen hosplogodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0
gen invlog_validation1_m0 = invlogit(hosplogodd_validation1_m0)

gen hi_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 + 1.4*stdp_validation1_m0
gen lo_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 - 1.4*stdp_validation1_m0

gen hi_invlog_validation1_m0 = invlogit(hi_logodd_validation1_m0)
gen lo_invlog_validation1_m0 = invlogit(lo_logodd_validation1_m0)

brier mort90_edis pr_validation1_m0 


* standardized mortality ratio
total mort90_edis if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m0 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm0=e(b)[1,1]
	di `prvalidation1deathm0'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm0'	
	
* excess deaths 
di `validation1death' - `prvalidation1deathm0'	
di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000


* VALIDATION *
estimates restore derivation_m0

predict xb_validation2_m0 if validation2==1, xb 
	replace xb_validation2_m0=. if validation2!=1
predict pr_validation2_m0 if validation2==1	
brier mort90_edis pr_validation2_m0 

* standardized mortality ratio
total mort90_edis if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m0 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm0=e(b)[1,1]
	di `prvalidation2deathm0'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm0'	
	
* excess deaths 
di `validation2death' - `prvalidation2deathm0'	
di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' if derivation==1 || sta6a:  
estimates store derivation_m1
predict pr_derivation_m1 if e(sample)	
predict re_derivation_m1 if e(sample), reffects reses(rese_derivation_m1)
predict xb_derivation_m1 if e(sample), xb
	replace xb_derivation_m1=. if derivation!=1

egen meanxb_derivation_m1 = mean(xb_derivation_m1) if e(sample)
gen hosplogodd_derivation_m1 = meanxb_derivation_m1 + re_derivation_m1
gen invlog_derivation_m1 = invlogit(hosplogodd_derivation_m1)

brier mort90_edis pr_derivation_m1 

// standardized mortality ratio
total mort90_edis if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m1 if derivation==1 
	matrix list e(b)
	local prderivationdeathm1=e(b)[1,1]
	di `prderivationdeathm1'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm1'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm1'	
di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000


preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m1 = pr_derivation_m1, nq(`ngroups')
	collapse (sum) mort90_edis pr_derivation_m1 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m1)
	list numhosps mort90_edis pr_derivation_m1, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m1
predict pr_validation1_m1 if validation1==1	
predict re_validation1_m1 if validation1==1, reffects reses(rese_validation1_m1)
predict xb_validation1_m1 if validation1==1, xb
	replace xb_validation1_m1=. if validation1!=1
predict stdp_validation1_m1 if validation1==1, stdp 
	replace stdp_validation1_m1=. if validation1!=1

egen meanxb_validation1_m1 = mean(xb_validation1_m1) if validation1==1
gen hosplogodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1
gen invlog_validation1_m1 = invlogit(hosplogodd_validation1_m1)

gen hi_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 + 1.4*stdp_validation1_m1
gen lo_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 - 1.4*stdp_validation1_m1

gen hi_invlog_validation1_m1 = invlogit(hi_logodd_validation1_m1)
gen lo_invlog_validation1_m1 = invlogit(lo_logodd_validation1_m1)

brier mort90_edis pr_validation1_m1


// standardized mortality ratio
total mort90_edis if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m1 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm1=e(b)[1,1]
	di `prvalidation1deathm1'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm1'	

// excess deaths 
di `validation1death' - `prvalidation1deathm1'	
di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m1
predict pr_validation2_m1 if validation2==1	
predict re_validation2_m1 if validation2==1, reffects reses(rese_validation2_m1)
predict xb_validation2_m1 if validation2==1, xb
	replace xb_validation2_m1=. if validation2!=1
predict stdp_validation2_m1 if validation2==1, stdp 
	replace stdp_validation2_m1=. if validation2!=1

egen meanxb_validation2_m1 = mean(xb_validation2_m1) if validation2==1
gen hosplogodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1
gen invlog_validation2_m1 = invlogit(hosplogodd_validation2_m1)

gen hi_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 + 1.4*stdp_validation2_m1
gen lo_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 - 1.4*stdp_validation2_m1

gen hi_invlog_validation2_m1 = invlogit(hi_logodd_validation2_m1)
gen lo_invlog_validation2_m1 = invlogit(lo_logodd_validation2_m1)

brier mort90_edis pr_validation2_m1

// standardized mortality ratio
total mort90_edis if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m1 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm1=e(b)[1,1]
	di `prvalidation2deathm1'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm1'	

// excess deaths 
di `validation2death' - `prvalidation2deathm1'	
di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*-------------------------------------------------------------------------------
* M2: Base + Individual (race, marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual'		///
		 if derivation==1 || sta6a: 
estimates store derivation_m2
predict pr_derivation_m2 if e(sample)	
predict re_derivation_m2 if e(sample), reffects reses(rese_derivation_m2)
predict xb_derivation_m2 if e(sample), xb
	replace xb_derivation_m2=. if derivation!=1

egen meanxb_derivation_m2 = mean(xb_derivation_m2) if e(sample)
gen hosplogodd_derivation_m2 = meanxb_derivation_m2 + re_derivation_m2
gen invlog_derivation_m2 = invlogit(hosplogodd_derivation_m2)

brier mort90_edis pr_derivation_m2 

// standardized mortality ratio
total mort90_edis if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m2 if derivation==1 
	matrix list e(b)
	local prderivationdeathm2=e(b)[1,1]
	di `prderivationdeathm2'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm2'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm2'	
di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m2
predict pr_validation1_m2 if validation1==1	
predict re_validation1_m2 if validation1==1, reffects reses(rese_validation1_m2)
predict xb_validation1_m2 if validation1==1, xb
	replace xb_validation1_m2=. if validation1!=1
predict stdp_validation1_m2 if validation1==1, stdp 
	replace stdp_validation1_m2=. if validation1!=1

egen meanxb_validation1_m2 = mean(xb_validation1_m2) if validation1==1
gen hosplogodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2
gen invlog_validation1_m2 = invlogit(hosplogodd_validation1_m2)

gen hi_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 + 1.4*stdp_validation1_m2
gen lo_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 - 1.4*stdp_validation1_m2

gen hi_invlog_validation1_m2 = invlogit(hi_logodd_validation1_m2)
gen lo_invlog_validation1_m2 = invlogit(lo_logodd_validation1_m2)

brier mort90_edis pr_validation1_m2


// standardized mortality ratio
total mort90_edis if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m2 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm2=e(b)[1,1]
	di `prvalidation1deathm2'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm2'	

// excess deaths 
di `validation1death' - `prvalidation1deathm2'	
di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m2
predict pr_validation2_m2 if validation2==1	
predict re_validation2_m2 if validation2==1, reffects reses(rese_validation2_m2)
predict xb_validation2_m2 if validation2==1, xb
	replace xb_validation2_m2=. if validation2!=1
predict stdp_validation2_m2 if validation2==1, stdp 
	replace stdp_validation2_m2=. if validation2!=1

egen meanxb_validation2_m2 = mean(xb_validation2_m2) if validation2==1
gen hosplogodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2
gen invlog_validation2_m2 = invlogit(hosplogodd_validation2_m2)

gen hi_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 + 1.4*stdp_validation2_m2
gen lo_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 - 1.4*stdp_validation2_m2

gen hi_invlog_validation2_m2 = invlogit(hi_logodd_validation2_m2)
gen lo_invlog_validation2_m2 = invlogit(lo_logodd_validation2_m2)

brier mort90_edis pr_validation2_m2

// standardized mortality ratio
total mort90_edis if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m2 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm2=e(b)[1,1]
	di `prvalidation2deathm2'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm2'	

// excess deaths 
di `validation2death' - `prvalidation2deathm2'	
di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000
	
	
*----------------------------
* M3: Base + Community 
*----------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' `community'		///	
		if derivation==1 || sta6a: 
estimates store derivation_m3
predict pr_derivation_m3 if e(sample)	
predict re_derivation_m3 if e(sample), reffects reses(rese_derivation_m3)
predict xb_derivation_m3 if e(sample), xb
	replace xb_derivation_m3=. if derivation!=1

egen meanxb_derivation_m3 = mean(xb_derivation_m3) if e(sample)
gen hosplogodd_derivation_m3 = meanxb_derivation_m3 + re_derivation_m3
gen invlog_derivation_m3 = invlogit(hosplogodd_derivation_m3)

brier mort90_edis pr_derivation_m3 

// standardized mortality ratio
total mort90_edis if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m3 if derivation==1 
	matrix list e(b)
	local prderivationdeathm3=e(b)[1,1]
	di `prderivationdeathm3'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm3'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm3'	
di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m3
predict pr_validation1_m3 if validation1==1	
predict re_validation1_m3 if validation1==1, reffects reses(rese_validation1_m3)
predict xb_validation1_m3 if validation1==1, xb
	replace xb_validation1_m3=. if validation1!=1
predict stdp_validation1_m3 if validation1==1, stdp 
	replace stdp_validation1_m3=. if validation1!=1

egen meanxb_validation1_m3 = mean(xb_validation1_m3) if validation1==1
gen hosplogodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3
gen invlog_validation1_m3 = invlogit(hosplogodd_validation1_m3)

gen hi_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 + 1.4*stdp_validation1_m3
gen lo_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 - 1.4*stdp_validation1_m3

gen hi_invlog_validation1_m3 = invlogit(hi_logodd_validation1_m3)
gen lo_invlog_validation1_m3 = invlogit(lo_logodd_validation1_m3)

brier mort90_edis pr_validation1_m3

// standardized mortality ratio
total mort90_edis if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m3 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm3=e(b)[1,1]
	di `prvalidation1deathm3'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm3'	

// excess deaths 
di `validation1death' - `prvalidation1deathm3'	
di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m3
predict pr_validation2_m3 if validation2==1	
predict re_validation2_m3 if validation2==1, reffects reses(rese_validation2_m3)
predict xb_validation2_m3 if validation2==1, xb
	replace xb_validation2_m3=. if validation2!=1
predict stdp_validation2_m3 if validation2==1, stdp 
	replace stdp_validation2_m3=. if validation2!=1

egen meanxb_validation2_m3 = mean(xb_validation2_m3) if validation2==1
gen hosplogodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3
gen invlog_validation2_m3 = invlogit(hosplogodd_validation2_m3)

gen hi_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 + 1.4*stdp_validation2_m3
gen lo_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 - 1.4*stdp_validation2_m3

gen hi_invlog_validation2_m3 = invlogit(hi_logodd_validation2_m3)
gen lo_invlog_validation2_m3 = invlogit(lo_logodd_validation2_m3)

brier mort90_edis pr_validation2_m3

// standardized mortality ratio
total mort90_edis if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m3 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm3=e(b)[1,1]
	di `prvalidation2deathm3'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm3'	

// excess deaths 
di `validation2death' - `prvalidation2deathm3'	
di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' ///
		if derivation==1 || sta6a: 
estimates store derivation_m4
predict pr_derivation_m4 if e(sample)	
predict re_derivation_m4 if e(sample), reffects reses(rese_derivation_m4)
predict xb_derivation_m4 if e(sample), xb
	replace xb_derivation_m4=. if derivation!=1

egen meanxb_derivation_m4 = mean(xb_derivation_m4) if e(sample)
gen hosplogodd_derivation_m4 = meanxb_derivation_m4 + re_derivation_m4
gen invlog_derivation_m4 = invlogit(hosplogodd_derivation_m4)

brier mort90_edis pr_derivation_m4 

// standardized mortality ratio
total mort90_edis if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m4 if derivation==1 
	matrix list e(b)
	local prderivationdeathm4=e(b)[1,1]
	di `prderivationdeathm4'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm4'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm4'	
di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m4 = pr_derivation_m4, nq(`ngroups')
	bysort quantile_derivation_m4: sum pr_derivation_m4, de
	collapse (sum) mort90_edis pr_derivation_m4 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m4)
	list numhosps mort90_edis pr_derivation_m4, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m4
predict pr_validation1_m4 if validation1==1	
predict re_validation1_m4 if validation1==1, reffects reses(rese_validation1_m4)
predict xb_validation1_m4 if validation1==1, xb
	replace xb_validation1_m4=. if validation1!=1
predict stdp_validation1_m4 if validation1==1, stdp 
	replace stdp_validation1_m4=. if validation1!=1

egen meanxb_validation1_m4 = mean(xb_validation1_m4) if validation1==1
gen hosplogodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4
gen invlog_validation1_m4 = invlogit(hosplogodd_validation1_m4)

gen hi_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 + 1.4*stdp_validation1_m4
gen lo_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 - 1.4*stdp_validation1_m4

gen hi_invlog_validation1_m4 = invlogit(hi_logodd_validation1_m4)
gen lo_invlog_validation1_m4 = invlogit(lo_logodd_validation1_m4)

brier mort90_edis pr_validation1_m4


// standardized mortality ratio
total mort90_edis if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m4 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm4=e(b)[1,1]
	di `prvalidation1deathm4'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm4'	

// excess deaths 
di `validation1death' - `prvalidation1deathm4'	
di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000

preserve 
	keep if validation1==1 
	local ngroups 10
	xtile quantile_validation1_m4 = pr_validation1_m4, nq(`ngroups')
	bysort quantile_validation1_m4: sum pr_validation1_m4, de
	collapse (sum) mort90_edis pr_validation1_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation1_m4)
	list numhosps mort90_edis pr_validation1_m4, abb(20)
restore


* VALIDATION 2021 * 
estimates restore derivation_m4
predict pr_validation2_m4 if validation2==1	
predict re_validation2_m4 if validation2==1, reffects reses(rese_validation2_m4)
predict xb_validation2_m4 if validation2==1, xb
	replace xb_validation2_m4=. if validation2!=1
predict stdp_validation2_m4 if validation2==1, stdp 
	replace stdp_validation2_m4=. if validation2!=1

egen meanxb_validation2_m4 = mean(xb_validation2_m4) if validation2==1
gen hosplogodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4
gen invlog_validation2_m4 = invlogit(hosplogodd_validation2_m4)

gen hi_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 + 1.4*stdp_validation2_m4
gen lo_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 - 1.4*stdp_validation2_m4

gen hi_invlog_validation2_m4 = invlogit(hi_logodd_validation2_m4)
gen lo_invlog_validation2_m4 = invlogit(lo_logodd_validation2_m4)

brier mort90_edis pr_validation2_m4

// standardized mortality ratio
total mort90_edis if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m4 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm4=e(b)[1,1]
	di `prvalidation2deathm4'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm4'	

// excess deaths 
di `validation2death' - `prvalidation2deathm4'	
di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


preserve 
	keep if validation2==1 
	local ngroups 10
	xtile quantile_validation2_m4 = pr_validation2_m4, nq(`ngroups')
	bysort quantile_validation2_m4: sum pr_validation2_m4, de
	collapse (sum) mort90_edis pr_validation2_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation2_m4)
	list numhosps mort90_edis pr_validation2_m4, abb(20)
restore


*-----------------------------------------
* Comparing models 0 & 1 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_edis xb_derivation_m0 xb_derivation_m1 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_edis xb_validation1_m0 xb_validation1_m1 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_edis xb_validation2_m0 xb_validation2_m1 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 2 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_edis xb_derivation_m1 xb_derivation_m2 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_edis xb_validation1_m1 xb_validation1_m2 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_edis xb_validation2_m1 xb_validation2_m2 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 3 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_edis xb_derivation_m1 xb_derivation_m3 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_edis xb_validation1_m1 xb_validation1_m3 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_edis xb_validation2_m1 xb_validation2_m3 	///
					if validation2==1
					
*-----------------------------------------
* Comparing models 1 & 4 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_edis xb_derivation_m1 xb_derivation_m4 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_edis xb_validation1_m1 xb_validation1_m4 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_edis xb_validation2_m1 xb_validation2_m4 	///
					if validation2==1
					

*-------------------------------------------------------------------------------
* NEW Figure 1: Comparison of base model and models 0, 2, 3, 4 in validation 1 
*-------------------------------------------------------------------------------			
					
preserve
	keep if validation1==1
	collapse re_validation1_m0 re_validation1_m1							///
			 re_validation1_m2 re_validation1_m3 re_validation1_m4			///
			 (mean) mean_validation1_m0=invlog_validation1_m0 				///
			 (mean) mean_validation1_m1=invlog_validation1_m1 				///
			 (mean) mean_validation1_m2=invlog_validation1_m2 				///
			 (mean) mean_validation1_m3=invlog_validation1_m3 				///
			 (mean) mean_validation1_m4=invlog_validation1_m4				///
			 (count) numhosps_m0=pr_validation1_m0							///
			 (count) numhosps_m1=pr_validation1_m1							///
			 (count) numhosps_m2=pr_validation1_m2							///
			 (count) numhosps_m3=pr_validation1_m3							///
			 (count) numhosps_m4=pr_validation1_m4							///
 			 (sum) sum_validation1_m0=invlog_validation1_m0 				///
			 (sum) sum_validation1_m1=invlog_validation1_m1 				///
			 (sum) sum_validation1_m2=invlog_validation1_m2 				///
			 (sum) sum_validation1_m3=invlog_validation1_m3 				///
			 (sum) sum_validation1_m4=invlog_validation1_m4					///
			  , by(sta6a)
	// random effects - M1 v M4
	twoway (scatter re_validation1_m4 re_validation1_m1, msymbol(Oh) msize(small))	///
		(lfit re_validation1_m4 re_validation1_m1),									///
		legend(off)	aspectratio(1)																		///
		ylab(0(.10).30, labsize(small) nogrid angle(0)) graphregion(color(white)) ymtick(##2)				///	
		xlab(0(.10).30, labsize(small))	xmtick(##2)											///
		ytitle("Hospital random effects in model 4", size(small) margin(medsmall))	///
		xtitle("Hospital random effects in model 1", size(small) margin(medsmall))	///
		name(scatter_validation1m4_re, replace)
	
	// predicted mortality from mixed model - M1 v M0
	twoway (scatter mean_validation1_m0 mean_validation1_m1, msymbol(Oh) msize(small))	///
		(lfit mean_validation1_m0 mean_validation1_m1),									///
		legend(off)	aspectratio(1)																		///
		ylab(0(.10).30, labsize(small) nogrid angle(0)) graphregion(color(white)) ymtick(##2)			///	
		xlab(0(.10).30, labsize(small))	xmtick(##2)											///
		ytitle("90-Day Adjusted Mortality Rate (ADI Only)", size(small) margin(medsmall))	///
		xtitle("90-Day Adjusted Mortality Rate (Base)", size(small) margin(medsmall))	///
		name(scatter_validation1m0_pr, replace)
	
	// predicted mortality from mixed model - M1 v M2
	twoway (scatter mean_validation1_m2 mean_validation1_m1, msymbol(Oh) msize(small))	///
		(lfit mean_validation1_m2 mean_validation1_m1),									///
		legend(off)	aspectratio(1)																		///
		ylab(0(.10).30, labsize(small) nogrid angle(0)) graphregion(color(white)) ymtick(##2)			///	
		xlab(0(.10).30, labsize(small))	xmtick(##2)											///
		ytitle("90-Day Adjusted Mortality Rate (Expanded)", size(small) margin(medsmall))	///
		xtitle("90-Day Adjusted Mortality Rate (Base)", size(small) margin(medsmall))	///
		name(scatter_validation1m2_pr, replace)
	
	// predicted mortality from mixed model - M1 v M3
	twoway (scatter mean_validation1_m3 mean_validation1_m1, msymbol(Oh) msize(small))	///
		(lfit mean_validation1_m3 mean_validation1_m1),									///
		legend(off)	aspectratio(1)																		///
		ylab(0(.10).30, labsize(small) nogrid angle(0)) graphregion(color(white)) ymtick(##2)			///	
		xlab(0(.10).30, labsize(small))	xmtick(##2)											///
		ytitle("90-Day Adjusted Mortality Rate (Base + Community SDOH)", size(small) margin(medsmall))	///
		xtitle("90-Day Adjusted Mortality Rate (Base)", size(small) margin(medsmall))	///
		name(scatter_validation1m3_pr, replace)
	
	// predicted mortality from mixed model - M1 v M4
	twoway (scatter mean_validation1_m4 mean_validation1_m1, msymbol(Oh) msize(small))	///
		(lfit mean_validation1_m4 mean_validation1_m1),									///
		legend(off)	aspectratio(1)																		///
		ylab(0(.10).30, labsize(small) nogrid angle(0)) graphregion(color(white)) ymtick(##2)			///	
		xlab(0(.10).30, labsize(small))	xmtick(##2)											///
		ytitle("90-Day Adjusted Mortality Rate (All)", size(small) margin(medsmall))	///
		xtitle("90-Day Adjusted Mortality Rate (Base)", size(small) margin(medsmall))	///
		name(scatter_validation1m4_pr, replace)

	
	graph save "scatter_validation1m0_pr" "SDOH\Figures\scatter_validation1m0_pr_20231221", replace 
	graph save "scatter_validation1m2_pr" "SDOH\Figures\scatter_validation1m2_pr_20231221", replace 
	graph save "scatter_validation1m3_pr" "SDOH\Figures\scatter_validation1m3_pr_20231221", replace 
	graph save "scatter_validation1m4_pr" "SDOH\Figures\scatter_validation1m4_pr_20231221", replace 
	
restore


*-------------------------------------------------------------------------------
* Figure: Comparison of hospital rankings with and without social determinants 
*-------------------------------------------------------------------------------

* Internal Validation Sample
preserve 
	collapse invlog_validation1_m1 invlog_validation1_m4 			///
			 hi_invlog_validation1_m1 lo_invlog_validation1_m1		///
			 hi_invlog_validation1_m4 lo_invlog_validation1_m4, by(sta6a)
	sort invlog_validation1_m1
	egen hosprank_m1 = rank(invlog_validation1_m1)
	sort invlog_validation1_m4 
	egen hosprank_m4 = rank(invlog_validation1_m4)
	br sta6a hosprank_m1 hosprank_m4
	xtile validation_m1_quintile = invlog_validation1_m1, nq(5) 
	xtile validation_m4_quintile = invlog_validation1_m4, nq(5) 
	sort validation_m1_quintile validation_m4_quintile
	br sta6a validation_m1_quintile validation_m4_quintile
	
	* create 3-category variable for color-coding quintiles in model 1
	gen model_1 = 0
	replace model_1 = 1 if validation_m1_quintile==1
	replace model_1 = 2 if validation_m1_quintile==5
	
	tab model_1 validation_m1_quintile
	
	* median line
	sum invlog_validation1_m4, de
	return list
	local m4median = r(p50)
	display `m4median'
	
	* caterpillar plot 
	sort invlog_validation1_m4 
		
	twoway rspike hi_invlog_validation1_m4 lo_invlog_validation1_m4 hosprank_m4 if model_1==0, 	///
				lstyle(ci) lcolor(gs8) ||		///
		   rspike hi_invlog_validation1_m4 lo_invlog_validation1_m4 hosprank_m4 if model_1==1 , 	///
				lstyle(ci) lcolor(lime) ||		///
		   rspike hi_invlog_validation1_m4 lo_invlog_validation1_m4 hosprank_m4 if model_1==2 , 	///
				lstyle(ci) lcolor(orange_red) ||		///		
		   scatter invlog_validation1_m4 hosprank_m4 if model_1==0,		///
			   yline(`m4median', lstyle(foreground))	///
			   mstyle(p1) msize(tiny) mcolor(gs8)	///
			   legend(off)	plotregion(lstyle(none))	///
			   ylab(0(.10).30, labsize(small) nogrid  angle(0)) graphregion(color(white)) ymtick(##2) ///			   
			   xlab(0 "1" 27 "2" 53 "3" 78 "4" 104 "5", grid labsize(vsmall)) 			///
			   ytitle("90-Day Adjusted Mortality Rate", size(small) margin(medsmall))  ///
			   xtitle("Hospital Rank by Quintile", size(small) margin(medsmall)) 	||  ///
		   scatter invlog_validation1_m4 hosprank_m4 if model_1==1,		///
			   yline(`m4median', lstyle(foreground))	///
			   mstyle(p1) msize(tiny) mcolor(lime)	///
			   legend(off)	plotregion(lstyle(none))	///
			   ylab(0(.10).30, labsize(small) nogrid  angle(0)) graphregion(color(white)) ymtick(##2)	///
			   xlab(0 "1" 27 "2" 53 "3" 78 "4" 104 "5", grid labsize(vsmall)) 			///
			   ytitle("90-Day Adjusted Mortality Rate", size(small) margin(medsmall)) ///
			   xtitle("Hospital Rank by Quintile", size(small) margin(medsmall))  ||  ///
		   scatter invlog_validation1_m4 hosprank_m4 if model_1==2,		///
			   yline(`m4median', lstyle(foreground))	///
			   mstyle(p1) msize(tiny) mcolor(orange_red)	///
			   legend(off)	plotregion(lstyle(none))	///
			   ylab(0(.10).30, labsize(small) nogrid  angle(0)) graphregion(color(white)) ymtick(##2)	///
			   xlab(0 "1" 27 "2" 53 "3" 78 "4" 104 "5", grid labsize(vsmall)) 			///
			   ytitle("90-Day Adjusted Mortality Rate", size(small) margin(medsmall)) ///
			   xtitle("Hospital Rank by Quintile", size(small) margin(medsmall))    ///	   
			   name(caterpillar, replace)
	
	graph save "caterpillar" "SDOH\Figures\caterpillar_20231221", replace 
	
restore		

		
********************************************************************************
* 90-day mortality to obtain ORs for covariates	using overall data set 	
********************************************************************************

* M0: ADI only

melogit mort90_edis adi_impute || sta6a: , or

* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)

* model using spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' || sta6a:  , or

* model using age_cat 
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male ib3.age_cat `sirs' `aod' `comorbid' || sta6a:  , or


* M2: Base + Individual (race, marital, housing instability, frailty, CAN) 

* model using spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' || sta6a: , or

* model using age cat 
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_edis male ib3.age_cat i.race_3cat `sirs' `aod' `comorbid' `individual' || sta6a: , or

* M3: Base + Community 

* model using age spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' `community' || sta6a: , or

* model using age cat 
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort90_edis male ib3.age_cat `sirs' `aod' `comorbid' `community' || sta6a: , or

* M4: Base + Expanded & Community 

* model using age spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

* model using age cat
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
		
melogit mort90_edis male ib3.age_cat i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

	
********************************************************************************
* 90-day post-discharge mortality to obtain ORs for covariates using overall 
* data set 	
********************************************************************************

* M4: Base + Expanded & Community 

* model using age spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort90_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

* model using age cat
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort90_postdc male ib3.age_cat i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

	
********************************************************************************
* 360-day post-discharge mortality to obtain ORs for covariates using overall 
* data set 	
********************************************************************************

* M4: Base + Expanded & Community 

* model using age spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort360_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

* model using age cat
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort360_postdc male ib3.age_cat i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

	
********************************************************************************
* 720-day post-discharge mortality to obtain ORs for covariates using overall 
* data set 	
********************************************************************************

* M4: Base + Expanded & Community 

* model using age spline
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
	
melogit mort720_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or

* model using age cat
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort720_postdc male ib3.age_cat i.race_3cat `sirs' `aod' `comorbid' `individual' `community' || sta6a: , or


********************************************************************************
* Subgroup analyses for derivation, internal validation, validation samples
********************************************************************************

*------------------
* AGE >=65 years 
*------------------

*** M0: Community Only ***

* DERIVATION *
melogit mort90_edis adi_impute if derivation==1 & age_65plus==1 || sta6a: 
estimates store deriv65plus_m0

predict pr_deriv65plus_m0 if e(sample)	
brier mort90_edis pr_deriv65plus_m0 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & age_65plus==1 
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_deriv65plus_m0 if derivation==1 & age_65plus==1 
		matrix list e(b)
		local prderivationdeathm0=e(b)[1,1]
		di `prderivationdeathm0'	
	total derivation if derivation==1 & age_65plus==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm0'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm0'	
	di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore deriv65plus_m0

predict pr_valid165plus_m0 if validation1==1 & age_65plus==1
brier mort90_edis pr_valid165plus_m0 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & age_65plus==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid165plus_m0 if validation1==1 & age_65plus==1
		matrix list e(b)
		local prvalidation1deathm0=e(b)[1,1]
		di `prvalidation1deathm0'	
	total validation1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm0'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm0'	
	di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000

* VALIDATION *
estimates restore deriv65plus_m0

predict pr_valid265plus_m0 if validation2==1 & age_65plus==1
brier mort90_edis pr_valid265plus_m0 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & age_65plus==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid265plus_m0 if validation2==1 & age_65plus==1
		matrix list e(b)
		local prvalidation2deathm0=e(b)[1,1]
		di `prvalidation2deathm0'	
	total validation2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm0'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm0'	
	di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000

	
*** M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria) ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male `sirs' `aod' `comorbid' if derivation==1 & age_65plus==1 || sta6a: 
estimates store deriv65plus_m1

predict pr_deriv65plus_m1 if e(sample)	
brier mort90_edis pr_deriv65plus_m1 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & age_65plus==1 
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_deriv65plus_m1 if derivation==1 & age_65plus==1
		matrix list e(b)
		local prderivationdeathm1=e(b)[1,1]
		di `prderivationdeathm1'	
	total derivation if derivation==1 & age_65plus==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm1'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm1'	
	di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore deriv65plus_m1

predict pr_valid165plus_m1 if validation1==1 & age_65plus==1
brier mort90_edis pr_valid165plus_m1 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & age_65plus==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid165plus_m1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local prvalidation1deathm1=e(b)[1,1]
		di `prvalidation1deathm1'	
	total validation1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm1'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm1'	
	di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000


* VALIDATION *
estimates restore deriv65plus_m1

predict pr_valid265plus_m1 if validation2==1 & age_65plus==1
brier mort90_edis pr_valid265plus_m1 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & age_65plus==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid265plus_m1 if validation2==1 & age_65plus==1
		matrix list e(b)
		local prvalidation2deathm1=e(b)[1,1]
		di `prvalidation2deathm1'	
	total validation2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm1'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm1'	
	di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000
	

*** M2: Base + Individual (race, marital, housing instability, frailty, CAN)  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
melogit mort90_edis male i.race_3cat `sirs' `aod' `comorbid' `individual' if derivation==1 & age_65plus==1 || sta6a: 
estimates store deriv65plus_m2

predict pr_deriv65plus_m2 if e(sample)	
brier mort90_edis pr_deriv65plus_m2 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & age_65plus==1 
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_deriv65plus_m2 if derivation==1 & age_65plus==1
		matrix list e(b)
		local prderivationdeathm2=e(b)[1,1]
		di `prderivationdeathm2'	
	total derivation if derivation==1 & age_65plus==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm2'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm2'	
	di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000
	
	
* INTERNAL VALIDATION *
estimates restore deriv65plus_m2

predict pr_valid165plus_m2 if validation1==1 & age_65plus==1
brier mort90_edis pr_valid165plus_m2 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & age_65plus==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid165plus_m2 if validation1==1 & age_65plus==1
		matrix list e(b)
		local prvalidation1deathm2=e(b)[1,1]
		di `prvalidation1deathm2'	
	total validation1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm2'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm2'	
	di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000


* VALIDATION *
estimates restore deriv65plus_m2

predict pr_valid265plus_m2 if validation2==1 & age_65plus==1
brier mort90_edis pr_valid265plus_m2 
	
	* standardized mortality ratio
	total mort90_edis if validation2==1 & age_65plus==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid265plus_m2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local prvalidation2deathm2=e(b)[1,1]
		di `prvalidation2deathm2'	
	total validation2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm2'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm2'	
	di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000
	

*** M3: Base + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
		
local community 															///	
			rural adi_impute
			
melogit mort90_edis male `sirs' `aod' `comorbid' `community' 	///
			if derivation==1 & age_65plus==1 || sta6a: 
estimates store deriv65plus_m3

predict pr_deriv65plus_m3 if e(sample)	
brier mort90_edis pr_deriv65plus_m3 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & age_65plus==1 
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_deriv65plus_m3 if derivation==1 & age_65plus==1
		matrix list e(b)
		local prderivationdeathm3=e(b)[1,1]
		di `prderivationdeathm3'	
	total derivation if derivation==1 & age_65plus==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm3'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm3'	
	di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000
	
* INTERNAL VALIDATION *
estimates restore deriv65plus_m3

predict pr_valid165plus_m3 if validation1==1 & age_65plus==1
brier mort90_edis pr_valid165plus_m3 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & age_65plus==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid165plus_m3 if validation1==1 & age_65plus==1
		matrix list e(b)
		local prvalidation1deathm3=e(b)[1,1]
		di `prvalidation1deathm3'	
	total validation1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm3'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm3'	
	di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000

* VALIDATION *
estimates restore deriv65plus_m3

predict pr_valid265plus_m3 if validation2==1 & age_65plus==1
brier mort90_edis pr_valid265plus_m3 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & age_65plus==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid265plus_m3 if validation2==1 & age_65plus==1
		matrix list e(b)
		local prvalidation2deathm3=e(b)[1,1]
		di `prvalidation2deathm3'	
	total validation2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm3'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm3'	
	di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000
	
*** M4: Base + Individual + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis male i.race_3cat `sirs' `aod' `comorbid' `individual' `community' 	///
			if derivation==1 & age_65plus==1 || sta6a: 
estimates store deriv65plus_m4

predict pr_deriv65plus_m4 if e(sample)	
brier mort90_edis pr_deriv65plus_m4 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & age_65plus==1 
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_deriv65plus_m4 if derivation==1 & age_65plus==1
		matrix list e(b)
		local prderivationdeathm4=e(b)[1,1]
		di `prderivationdeathm4'	
	total derivation if derivation==1 & age_65plus==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm4'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm4'	
	di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore deriv65plus_m4

predict pr_valid165plus_m4 if validation1==1 & age_65plus==1
brier mort90_edis pr_valid165plus_m4 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & age_65plus==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid165plus_m4 if validation1==1 & age_65plus==1
		matrix list e(b)
		local prvalidation1deathm4=e(b)[1,1]
		di `prvalidation1deathm4'	
	total validation1 if validation1==1 & age_65plus==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm4'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm4'	
	di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000

* VALIDATION *
estimates restore deriv65plus_m4

predict pr_valid265plus_m4 if validation2==1 & age_65plus==1
brier mort90_edis pr_valid265plus_m4 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & age_65plus==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid265plus_m4 if validation2==1 & age_65plus==1
		matrix list e(b)
		local prvalidation2deathm4=e(b)[1,1]
		di `prvalidation2deathm4'	
	total validation2 if validation2==1 & age_65plus==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm4'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm4'	
	di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


*------------------
* Women 
*------------------

*** M0: Community Only ***

* DERIVATION *
melogit mort90_edis adi_impute if derivation==1 & female==1 || sta6a: 
estimates store derivfemale_m0

predict pr_derivfemale_m0 if e(sample)	
brier mort90_edis pr_derivfemale_m0 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & female==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivfemale_m0 if derivation==1 & female==1
		matrix list e(b)
		local prderivationdeathm0=e(b)[1,1]
		di `prderivationdeathm0'	
	total derivation if derivation==1 & female==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm0'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm0'	
	di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivfemale_m0

predict pr_valid1female_m0 if validation1==1 & female==1
brier mort90_edis pr_valid1female_m0 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & female==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1female_m0 if validation1==1 & female==1
		matrix list e(b)
		local prvalidation1deathm0=e(b)[1,1]
		di `prvalidation1deathm0'	
	total validation1 if validation1==1 & female==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm0'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm0'	
	di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000

* VALIDATION *
estimates restore derivfemale_m0

predict pr_valid2female_m0 if validation2==1 & female==1
brier mort90_edis pr_valid2female_m0 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & female==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2female_m0 if validation2==1 & female==1
		matrix list e(b)
		local prvalidation2deathm0=e(b)[1,1]
		di `prvalidation2deathm0'	
	total validation2 if validation2==1 & female==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm0'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm0'	
	di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000

	
*** M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria) ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis age_sp* `sirs' `aod' `comorbid' if derivation==1 & female==1 || sta6a: 
estimates store derivfemale_m1

predict pr_derivfemale_m1 if e(sample)	
brier mort90_edis pr_derivfemale_m1 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & female==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivfemale_m1 if derivation==1 & female==1
		matrix list e(b)
		local prderivationdeathm1=e(b)[1,1]
		di `prderivationdeathm1'	
	total derivation if derivation==1 & female==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm1'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm1'	
	di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivfemale_m1

predict pr_valid1female_m1 if validation1==1 & female==1
brier mort90_edis pr_valid1female_m1 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & female==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1female_m1 if validation1==1 & female==1
		matrix list e(b)
		local prvalidation1deathm1=e(b)[1,1]
		di `prvalidation1deathm1'	
	total validation1 if validation1==1 & female==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm1'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm1'	
	di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000

* VALIDATION *
estimates restore derivfemale_m1

predict pr_valid2female_m1 if validation2==1 & female==1
brier mort90_edis pr_valid2female_m1 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & female==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2female_m1 if validation2==1 & female==1
		matrix list e(b)
		local prvalidation2deathm1=e(b)[1,1]
		di `prvalidation2deathm1'	
	total validation2 if validation2==1 & female==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm1'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm1'	
	di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000

	
*** M2: Base + Individual (marital, housing instability, frailty, CAN)  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
melogit mort90_edis age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' if derivation==1 & female==1 || sta6a: 
estimates store derivfemale_m2

predict pr_derivfemale_m2 if e(sample)	
brier mort90_edis pr_derivfemale_m2 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & female==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivfemale_m2 if derivation==1 & female==1
		matrix list e(b)
		local prderivationdeathm2=e(b)[1,1]
		di `prderivationdeathm2'	
	total derivation if derivation==1 & female==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm2'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm2'	
	di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivfemale_m2

predict pr_valid1female_m2 if validation1==1 & female==1
brier mort90_edis pr_valid1female_m2 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & female==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1female_m2 if validation1==1 & female==1
		matrix list e(b)
		local prvalidation1deathm2=e(b)[1,1]
		di `prvalidation1deathm2'	
	total validation1 if validation1==1 & female==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm2'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm2'	
	di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000

* VALIDATION *
estimates restore derivfemale_m2

predict pr_valid2female_m2 if validation2==1 & female==1
brier mort90_edis pr_valid2female_m2 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & female==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2female_m2 if validation2==1 & female==1
		matrix list e(b)
		local prvalidation2deathm2=e(b)[1,1]
		di `prvalidation2deathm2'	
	total validation2 if validation2==1 & female==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm2'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm2'	
	di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000


*** M3: Base + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
		
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* `sirs' `aod' `comorbid' `community' 	///
			if derivation==1 & female==1 || sta6a: 
estimates store derivfemale_m3

predict pr_derivfemale_m3 if e(sample)	
brier mort90_edis pr_derivfemale_m3 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & female==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivfemale_m3 if derivation==1 & female==1
		matrix list e(b)
		local prderivationdeathm3=e(b)[1,1]
		di `prderivationdeathm3'	
	total derivation if derivation==1 & female==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm3'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm3'	
	di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivfemale_m3

predict pr_valid1female_m3 if validation1==1 & female==1
brier mort90_edis pr_valid1female_m3 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & female==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1female_m3 if validation1==1 & female==1
		matrix list e(b)
		local prvalidation1deathm3=e(b)[1,1]
		di `prvalidation1deathm3'	
	total validation1 if validation1==1 & female==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm3'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm3'	
	di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000

* VALIDATION *
estimates restore derivfemale_m3

predict pr_valid2female_m3 if validation2==1 & female==1
brier mort90_edis pr_valid2female_m3 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & female==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2female_m3 if validation2==1 & female==1
		matrix list e(b)
		local prvalidation2deathm3=e(b)[1,1]
		di `prvalidation2deathm3'	
	total validation2 if validation2==1 & female==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm3'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm3'	
	di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000


*** M4: Base + Individual + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' 	///
			if derivation==1 & female==1 || sta6a: 
estimates store derivfemale_m4

predict pr_derivfemale_m4 if e(sample)	
brier mort90_edis pr_derivfemale_m4 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & female==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivfemale_m4 if derivation==1 & female==1
		matrix list e(b)
		local prderivationdeathm4=e(b)[1,1]
		di `prderivationdeathm4'	
	total derivation if derivation==1 & female==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm4'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm4'	
	di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivfemale_m4

predict pr_valid1female_m4 if validation1==1 & female==1
brier mort90_edis pr_valid1female_m4 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & female==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1female_m4 if validation1==1 & female==1
		matrix list e(b)
		local prvalidation1deathm4=e(b)[1,1]
		di `prvalidation1deathm4'	
	total validation1 if validation1==1 & female==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm4'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm4'	
	di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000


* VALIDATION *
estimates restore derivfemale_m4

predict pr_valid2female_m4 if validation2==1 & female==1
brier mort90_edis pr_valid2female_m4 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & female==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2female_m4 if validation2==1 & female==1
		matrix list e(b)
		local prvalidation2deathm4=e(b)[1,1]
		di `prvalidation2deathm4'	
	total validation2 if validation2==1 & female==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm4'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm4'	
	di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


*------------------
* Black race
*------------------

*** M0: Community Only ***

* DERIVATION *
melogit mort90_edis adi_impute if derivation==1 & black==1 || sta6a: 
estimates store derivblack_m0

predict pr_derivblack_m0 if e(sample)	
brier mort90_edis pr_derivblack_m0 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & black==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivblack_m0 if derivation==1 & black==1
		matrix list e(b)
		local prderivationdeathm0=e(b)[1,1]
		di `prderivationdeathm0'	
	total derivation if derivation==1 & black==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm0'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm0'	
	di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000

						
* INTERNAL VALIDATION *
estimates restore derivblack_m0

predict pr_valid1black_m0 if validation1==1 & black==1
brier mort90_edis pr_valid1black_m0 
	
	* standardized mortality ratio
	total mort90_edis if validation1==1 & black==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1black_m0 if validation1==1 & black==1
		matrix list e(b)
		local prvalidation1deathm0=e(b)[1,1]
		di `prvalidation1deathm0'	
	total validation1 if validation1==1 & black==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm0'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm0'	
	di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000
			
* VALIDATION *
estimates restore derivblack_m0

predict pr_valid2black_m0 if validation2==1 & black==1
brier mort90_edis pr_valid2black_m0 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & black==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2black_m0 if validation2==1 & black==1
		matrix list e(b)
		local prvalidation2deathm0=e(b)[1,1]
		di `prvalidation2deathm0'	
	total validation2 if validation2==1 & black==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm0'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm0'	
	di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*** M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria) ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' if derivation==1 & black==1 || sta6a: 
estimates store derivblack_m1

predict pr_derivblack_m1 if e(sample)	
brier mort90_edis pr_derivblack_m1 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & black==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivblack_m1 if derivation==1 & black==1
		matrix list e(b)
		local prderivationdeathm1=e(b)[1,1]
		di `prderivationdeathm1'	
	total derivation if derivation==1 & black==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm1'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm1'	
	di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivblack_m1

predict pr_valid1black_m1 if validation1==1 & black==1
brier mort90_edis pr_valid1black_m1 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & black==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1black_m1 if validation1==1 & black==1
		matrix list e(b)
		local prvalidation1deathm1=e(b)[1,1]
		di `prvalidation1deathm1'	
	total validation1 if validation1==1 & black==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm1'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm1'	
	di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000
	
* VALIDATION *
estimates restore derivblack_m1

predict pr_valid2black_m1 if validation2==1 & black==1
brier mort90_edis pr_valid2black_m1 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & black==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2black_m1 if validation2==1 & black==1
		matrix list e(b)
		local prvalidation2deathm1=e(b)[1,1]
		di `prvalidation2deathm1'	
	total validation2 if validation2==1 & black==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm1'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm1'	
	di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*** M2: Base + Individual (marital, housing instability, frailty, CAN)  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' `individual' if derivation==1 & black==1 || sta6a: 
estimates store derivblack_m2

predict pr_derivblack_m2 if e(sample)	
brier mort90_edis pr_derivblack_m2 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & black==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivblack_m2 if derivation==1 & black==1
		matrix list e(b)
		local prderivationdeathm2=e(b)[1,1]
		di `prderivationdeathm2'	
	total derivation if derivation==1 & black==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm2'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm2'	
	di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivblack_m2

predict pr_valid1black_m2 if validation1==1 & black==1
brier mort90_edis pr_valid1black_m2 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & black==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1black_m2 if validation1==1 & black==1
		matrix list e(b)
		local prvalidation1deathm2=e(b)[1,1]
		di `prvalidation1deathm2'	
	total validation1 if validation1==1 & black==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm2'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm2'	
	di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000
	
* VALIDATION *
estimates restore derivblack_m2

predict pr_valid2black_m2 if validation2==1 & black==1
brier mort90_edis pr_valid2black_m2 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & black==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2black_m2 if validation2==1 & black==1
		matrix list e(b)
		local prvalidation2deathm2=e(b)[1,1]
		di `prvalidation2deathm2'	
	total validation2 if validation2==1 & black==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm2'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm2'	
	di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000


*** M3: Base + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
		
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' `community' 	///
			if derivation==1 & black==1 || sta6a: 
estimates store derivblack_m3

predict pr_derivblack_m3 if e(sample)	
brier mort90_edis pr_derivblack_m3 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & black==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivblack_m3 if derivation==1 & black==1
		matrix list e(b)
		local prderivationdeathm3=e(b)[1,1]
		di `prderivationdeathm3'	
	total derivation if derivation==1 & black==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm3'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm3'	
	di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivblack_m3

predict pr_valid1black_m3 if validation1==1 & black==1
brier mort90_edis pr_valid1black_m3 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & black==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1black_m3 if validation1==1 & black==1
		matrix list e(b)
		local prvalidation1deathm3=e(b)[1,1]
		di `prvalidation1deathm3'	
	total validation1 if validation1==1 & black==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm3'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm3'	
	di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000
	
* VALIDATION *
estimates restore derivblack_m3

predict pr_valid2black_m3 if validation2==1 & black==1
brier mort90_edis pr_valid2black_m3 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & black==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2black_m3 if validation2==1 & black==1
		matrix list e(b)
		local prvalidation2deathm3=e(b)[1,1]
		di `prvalidation2deathm3'	
	total validation2 if validation2==1 & black==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm3'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm3'	
	di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000


*** M4: Base + Individual + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' `individual' `community' 	///
			if derivation==1 & black==1 || sta6a: 
estimates store derivblack_m4

predict pr_derivblack_m4 if e(sample)	
brier mort90_edis pr_derivblack_m4 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & black==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivblack_m4 if derivation==1 & black==1
		matrix list e(b)
		local prderivationdeathm4=e(b)[1,1]
		di `prderivationdeathm4'	
	total derivation if derivation==1 & black==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm4'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm4'	
	di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

* INTERNAL VALIDATION *
estimates restore derivblack_m4

predict pr_valid1black_m4 if validation1==1 & black==1
brier mort90_edis pr_valid1black_m4 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & black==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1black_m4 if validation1==1 & black==1
		matrix list e(b)
		local prvalidation1deathm4=e(b)[1,1]
		di `prvalidation1deathm4'	
	total validation1 if validation1==1 & black==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm4'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm4'	
	di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000
	
* VALIDATION *
estimates restore derivblack_m4

predict pr_valid2black_m4 if validation2==1 & black==1
brier mort90_edis pr_valid2black_m4 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & black==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2black_m4 if validation2==1 & black==1
		matrix list e(b)
		local prvalidation2deathm4=e(b)[1,1]
		di `prvalidation2deathm4'	
	total validation2 if validation2==1 & black==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm4'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm4'	
	di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


*------------------
* COVID+
*------------------

*** M0: Community Only ***

* DERIVATION *
melogit mort90_edis adi_impute if derivation==1 & covid_diag_hosp==1 || sta6a: 
estimates store derivcovid_m0

predict pr_derivcovid_m0 if e(sample)	
brier mort90_edis pr_derivcovid_m0 
	
	* standardized mortality ratio
	total mort90_edis if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivcovid_m0 if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local prderivationdeathm0=e(b)[1,1]
		di `prderivationdeathm0'	
	total derivation if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm0'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm0'	
	di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000
					
* INTERNAL VALIDATION *
estimates restore derivcovid_m0

predict pr_valid1covid_m0 if validation1==1 & covid_diag_hosp==1
brier mort90_edis pr_valid1covid_m0 
	
	* standardized mortality ratio
	total mort90_edis if validation1==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1covid_m0 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation1deathm0=e(b)[1,1]
		di `prvalidation1deathm0'	
	total validation1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm0'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm0'	
	di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000
				
* VALIDATION *
estimates restore derivcovid_m0

predict pr_valid2covid_m0 if validation2==1 & covid_diag_hosp==1
brier mort90_edis pr_valid2covid_m0 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2covid_m0 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation2deathm0=e(b)[1,1]
		di `prvalidation2deathm0'	
	total validation2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm0'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm0'	
	di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*** M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria) ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' if derivation==1 & covid_diag_hosp==1 || sta6a: 
estimates store derivcovid_m1

predict pr_derivcovid_m1 if e(sample)	
brier mort90_edis pr_derivcovid_m1 
	
	* standardized mortality ratio
	total mort90_edis if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivcovid_m1 if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local prderivationdeathm1=e(b)[1,1]
		di `prderivationdeathm1'	
	total derivation if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm1'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm1'	
	di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000
		
* INTERNAL VALIDATION *
estimates restore derivcovid_m1

predict pr_valid1covid_m1 if validation1==1 & covid_diag_hosp==1
brier mort90_edis pr_valid1covid_m1 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1covid_m1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation1deathm1=e(b)[1,1]
		di `prvalidation1deathm1'	
	total validation1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm1'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm1'	
	di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000
			
* VALIDATION *
estimates restore derivcovid_m1

predict pr_valid2covid_m1 if validation2==1 & covid_diag_hosp==1
brier mort90_edis pr_valid2covid_m1 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2covid_m1 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation2deathm1=e(b)[1,1]
		di `prvalidation2deathm1'	
	total validation2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm1'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm1'	
	di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*** M2: Base + Individual (marital, housing instability, frailty, CAN)  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
melogit mort90_edis age_sp* male i.race_3cat `sirs' `aod' `comorbid' `individual' if derivation==1 & covid_diag_hosp==1 || sta6a: 
estimates store derivcovid_m2

predict pr_derivcovid_m2 if e(sample)	
brier mort90_edis pr_derivcovid_m2 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivcovid_m2 if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local prderivationdeathm2=e(b)[1,1]
		di `prderivationdeathm2'	
	total derivation if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm2'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm2'	
	di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000
	
* INTERNAL VALIDATION *
estimates restore derivcovid_m2

predict pr_valid1covid_m2 if validation1==1 & covid_diag_hosp==1
brier mort90_edis pr_valid1covid_m2 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1covid_m2 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation1deathm2=e(b)[1,1]
		di `prvalidation1deathm2'	
	total validation1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm2'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm2'	
	di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000
		
* VALIDATION *
estimates restore derivcovid_m2

predict pr_valid2covid_m2 if validation2==1 & covid_diag_hosp==1
brier mort90_edis pr_valid2covid_m2 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2covid_m2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation2deathm2=e(b)[1,1]
		di `prvalidation2deathm2'	
	total validation2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm2'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm2'	
	di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000


*** M3: Base + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
		
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* male `sirs' `aod' `comorbid' `community' 	///
			if derivation==1 & covid_diag_hosp==1 || sta6a: 
estimates store derivcovid_m3

predict pr_derivcovid_m3 if e(sample)	
brier mort90_edis pr_derivcovid_m3 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivcovid_m3 if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local prderivationdeathm3=e(b)[1,1]
		di `prderivationdeathm3'	
	total derivation if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm3'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm3'	
	di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000
	
* INTERNAL VALIDATION *
estimates restore derivcovid_m3

predict pr_valid1covid_m3 if validation1==1 & covid_diag_hosp==1
brier mort90_edis pr_valid1covid_m3 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1covid_m3 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation1deathm3=e(b)[1,1]
		di `prvalidation1deathm3'	
	total validation1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm3'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm3'	
	di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000
	
* VALIDATION *
estimates restore derivcovid_m3

predict pr_valid2covid_m3 if validation2==1 & covid_diag_hosp==1
brier mort90_edis pr_valid2covid_m3 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2covid_m3 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation2deathm3=e(b)[1,1]
		di `prvalidation2deathm3'	
	total validation2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm3'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm3'	
	di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000


*** M4: Base + Individual + Community  ***

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute	
			
local community 															///	
			rural adi_impute
			
melogit mort90_edis age_sp* male i.race_3cat `sirs' `aod' `comorbid' `individual' `community' 	///
			if derivation==1 & covid_diag_hosp==1 || sta6a: 
estimates store derivcovid_m4

predict pr_derivcovid_m4 if e(sample)	
brier mort90_edis pr_derivcovid_m4 

	* standardized mortality ratio
	total mort90_edis if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivationdeath=e(b)[1,1]
		di `derivationdeath'
	total pr_derivcovid_m4 if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local prderivationdeathm4=e(b)[1,1]
		di `prderivationdeathm4'	
	total derivation if derivation==1 & covid_diag_hosp==1
		matrix list e(b)
		local derivation=e(b)[1,1]
		di `derivation'

	di `derivationdeath'/`prderivationdeathm4'	
		
	* excess deaths 
	di `derivationdeath' - `prderivationdeathm4'	
	di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000
	
* INTERNAL VALIDATION *
estimates restore derivcovid_m4

predict pr_valid1covid_m4 if validation1==1 & covid_diag_hosp==1
brier mort90_edis pr_valid1covid_m4 

	* standardized mortality ratio
	total mort90_edis if validation1==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation1death=e(b)[1,1]
		di `validation1death'
	total pr_valid1covid_m4 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation1deathm4=e(b)[1,1]
		di `prvalidation1deathm4'	
	total validation1 if validation1==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation1=e(b)[1,1]
		di `validation1'

	di `validation1death'/`prvalidation1deathm4'	
		
	* excess deaths 
	di `validation1death' - `prvalidation1deathm4'	
	di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000
	
* VALIDATION *
estimates restore derivcovid_m4

predict pr_valid2covid_m4 if validation2==1 & covid_diag_hosp==1
brier mort90_edis pr_valid2covid_m4 

	* standardized mortality ratio
	total mort90_edis if validation2==1 & covid_diag_hosp==1 
		matrix list e(b)
		local validation2death=e(b)[1,1]
		di `validation2death'
	total pr_valid2covid_m4 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local prvalidation2deathm4=e(b)[1,1]
		di `prvalidation2deathm4'	
	total validation2 if validation2==1 & covid_diag_hosp==1
		matrix list e(b)
		local validation2=e(b)[1,1]
		di `validation2'

	di `validation2death'/`prvalidation2deathm4'	
		
	* excess deaths 
	di `validation2death' - `prvalidation2deathm4'	
	di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


*******************************
* Karandeep Calibration Plots *	 
*******************************

foreach var in derivation validation1 validation2 {
	preserve 
		keep if `var'==1
		keep unique_hosp_count_id mort90_edis pr_`var'_m4 
		count
		save SDOH\Data\karandeep_`var'_20240102, replace
	restore
}
	

************************
* Sensitivity Analysis *
************************

*---------------------------
*  Missingness of SDI/SVI 
*---------------------------

misstable sum sdi_score svi 
	 * sdi % missing: 13.16%
	 * svi % missing: 3.46%

gen sdimiss = missing(sdi_score)	 
	 
version 16: table admityear, c(count unique count sdi_score mean sdi_score mean sdimiss sum sdimiss) 
version 16: table admityear, c(count unique count svi mean svi mean svimiss sum svimiss) 

* use single multivariable regression imputation to impute missing values 

* SVI
local covar male age i.race_3cat											///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 					///
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			aod_heme acute_respiratory_hosp 								///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			housing_instab_combo_prior540 frailty	rural can_impute				
regress svi `covar' 		
predict mkg_svi_impute

gen svi_impute = svi
replace svi_impute = mkg_svi_impute if svi==. 

sum svi mkg_svi_impute svi_impute		
drop mkg_svi_impute 

sum svi svi_impute if derivation==1, de
sum svi svi_impute if validation1==1, de
sum svi svi_impute if validation2==1, de

* SDI 
local covar male age i.race_3cat											///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 					///
			aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			aod_heme acute_respiratory_hosp 								///
			sirs_temp sirs_rr sirs_pulse sirs_wbc							///
			housing_instab_combo_prior540 frailty_score	rural can_impute				
regress sdi_score `covar' 		
predict mkg_sdi_impute

gen sdi_impute = sdi_score
replace sdi_impute = mkg_sdi_impute if sdi_score==. 

sum sdi_score mkg_sdi_impute sdi_impute		
drop mkg_sdi_impute 

sum sdi_score sdi_impute if derivation==1, de
sum sdi_score sdi_impute if validation1==1, de
sum sdi_score sdi_impute if validation2==1, de


*****************************
* Model Performance
*****************************

** SVI **

*----------------------
* M0: Community Only
*----------------------

melogit mort90_edis svi_impute || sta6a: 
estimates store svi_m0

predict pr_svi_m0 if e(sample)	
brier mort90_edis pr_svi_m0 


* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local svideath=e(b)[1,1]
	di `svideath'
total pr_svi_m0
	matrix list e(b)
	local prsvideathm0=e(b)[1,1]
	di `prsvideathm0'	

di `svideath'/`prsvideathm0'	
	
* excess deaths 
di `svideath' - `prsvideathm0'	
di ((`svideath' - `prsvideathm0')/144889)*10000

*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid'  || sta6a:  
estimates store svi_m1

predict pr_svi_m1 if e(sample)	
brier mort90_edis pr_svi_m1 


* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local svideath=e(b)[1,1]
	di `svideath'
total pr_svi_m1
	matrix list e(b)
	local prsvideathm1=e(b)[1,1]
	di `prsvideathm1'	

di `svideath'/`prsvideathm1'	
	
* excess deaths 
di `svideath' - `prsvideathm1'	
di ((`svideath' - `prsvideathm1')/144889)*10000

*-------------------------------------------------------------------------------
* M2: Base + Individual (marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' || sta6a: 
estimates store svi_m2

predict pr_svi_m2 if e(sample)	
brier mort90_edis pr_svi_m2 
	
* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local svideath=e(b)[1,1]
	di `svideath'
total pr_svi_m2
	matrix list e(b)
	local prsvideathm2=e(b)[1,1]
	di `prsvideathm2'	

di `svideath'/`prsvideathm2'	
	
* excess deaths 
di `svideath' - `prsvideathm2'	
di ((`svideath' - `prsvideathm2')/144889)*10000

*----------------------------
* M3: Base + Community 
*----------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural svi_impute
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' `community' || sta6a: 
estimates store svi_m3

predict pr_svi_m3 if e(sample)	
brier mort90_edis pr_svi_m3 

* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local svideath=e(b)[1,1]
	di `svideath'
total pr_svi_m3
	matrix list e(b)
	local prsvideathm3=e(b)[1,1]
	di `prsvideathm3'	

di `svideath'/`prsvideathm3'	
	
* excess deaths 
di `svideath' - `prsvideathm3'	
di ((`svideath' - `prsvideathm3')/144889)*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural svi_impute
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community'  || sta6a: 
estimates store svi_m4

predict pr_svi_m4 if e(sample)	
brier mort90_edis pr_svi_m4 

* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local svideath=e(b)[1,1]
	di `svideath'
total pr_svi_m4
	matrix list e(b)
	local prsvideathm4=e(b)[1,1]
	di `prsvideathm4'	

di `svideath'/`prsvideathm4'	
	
* excess deaths 
di `svideath' - `prsvideathm4'	
di ((`svideath' - `prsvideathm4')/144889)*10000


** SDI **

*----------------------
* M0: SDI Only
*----------------------

melogit mort90_edis sdi_impute || sta6a: 
estimates store sdi_m0

predict pr_sdi_m0 if e(sample)	
brier mort90_edis pr_sdi_m0 


* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local sdideath=e(b)[1,1]
	di `sdideath'
total pr_sdi_m0
	matrix list e(b)
	local prsdideathm0=e(b)[1,1]
	di `prsdideathm0'	

di `sdideath'/`prsdideathm0'	
	
* excess deaths 
di `sdideath' - `prsdideathm0'	
di ((`sdideath' - `prsdideathm0')/144889)*10000

*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid'  || sta6a:  
estimates store sdi_m1

predict pr_sdi_m1 if e(sample)	
brier mort90_edis pr_sdi_m1 

* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local sdideath=e(b)[1,1]
	di `sdideath'
total pr_sdi_m1
	matrix list e(b)
	local prsdideathm1=e(b)[1,1]
	di `prsdideathm1'	

di `sdideath'/`prsdideathm1'	
	
* excess deaths 
di `sdideath' - `prsdideathm1'	
di ((`sdideath' - `prsdideathm1')/144889)*10000

*-------------------------------------------------------------------------------
* M2: Base + Individual (marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' || sta6a: 
estimates store sdi_m2

predict pr_sdi_m2 if e(sample)	
brier mort90_edis pr_sdi_m2 
	
* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local sdideath=e(b)[1,1]
	di `sdideath'
total pr_sdi_m2
	matrix list e(b)
	local prsdideathm2=e(b)[1,1]
	di `prsdideathm2'	

di `sdideath'/`prsdideathm2'	
	
* excess deaths 
di `sdideath' - `prsdideathm2'	
di ((`sdideath' - `prsdideathm2')/144889)*10000

*----------------------------
* M3: Base + Community 
*----------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural sdi_impute
			
melogit mort90_edis male age_sp* `sirs' `aod' `comorbid' `community' || sta6a: 
estimates store sdi_m3

predict pr_sdi_m3 if e(sample)	
brier mort90_edis pr_sdi_m3 

* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local sdideath=e(b)[1,1]
	di `sdideath'
total pr_sdi_m3
	matrix list e(b)
	local prsdideathm3=e(b)[1,1]
	di `prsdideathm3'	

di `sdideath'/`prsdideathm3'	
	
* excess deaths 
di `sdideath' - `prsdideathm3'	
di ((`sdideath' - `prsdideathm3')/144889)*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural sdi_impute
			
melogit mort90_edis male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community'  || sta6a: 
estimates store sdi_m4

predict pr_sdi_m4 if e(sample)	
brier mort90_edis pr_sdi_m4 

* standardized mortality ratio
total mort90_edis 
	matrix list e(b)
	local sdideath=e(b)[1,1]
	di `sdideath'
total pr_sdi_m4
	matrix list e(b)
	local prsdideathm4=e(b)[1,1]
	di `prsdideathm4'	

di `sdideath'/`prsdideathm4'	
	
* excess deaths 
di `sdideath' - `prsdideathm4'	
di ((`sdideath' - `prsdideathm4')/144889)*10000


*******************************************************************************
* Sensitivity Analysis: Comparison of Model Performance - Post-Discharge Mort *
*******************************************************************************
				
*-------------------------------------------------------------------------------
* Supplemental Table: Comparison of model performance - 90-day Post-Discharge
*-------------------------------------------------------------------------------

drop pr_derivation_m0-lo_invlog_validation2_m4

*----------------------
* M0: Community Only
*----------------------

* DERIVATION *
melogit mort90_postdc adi_impute if derivation==1 || sta6a: 
estimates store derivation_m0

predict pr_derivation_m0 if e(sample)
predict xb_derivation_m0 if e(sample), xb
brier mort90_postdc pr_derivation_m0 

* standardized mortality ratio
total mort90_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m0 if derivation==1 
	matrix list e(b)
	local prderivationdeathm0=e(b)[1,1]
	di `prderivationdeathm0'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm0'	
	
* excess deaths 
di `derivationdeath' - `prderivationdeathm0'	
di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m0

predict xb_validation1_m0 if validation1==1, xb 
	replace xb_validation1_m0=. if validation1!=1
predict pr_validation1_m0 if validation1==1	
predict re_validation1_m0 if validation1==1, reffects reses(rese_validation1_m0)
predict stdp_validation1_m0 if validation1==1, stdp 
	replace stdp_validation1_m0=. if validation1!=1

egen meanxb_validation1_m0 = mean(xb_validation1_m0) if validation1==1
gen hosplogodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0
gen invlog_validation1_m0 = invlogit(hosplogodd_validation1_m0)

gen hi_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 + 1.4*stdp_validation1_m0
gen lo_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 - 1.4*stdp_validation1_m0

gen hi_invlog_validation1_m0 = invlogit(hi_logodd_validation1_m0)
gen lo_invlog_validation1_m0 = invlogit(lo_logodd_validation1_m0)

brier mort90_postdc pr_validation1_m0 


* standardized mortality ratio
total mort90_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m0 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm0=e(b)[1,1]
	di `prvalidation1deathm0'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm0'	
	
* excess deaths 
di `validation1death' - `prvalidation1deathm0'	
di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000


* VALIDATION *
estimates restore derivation_m0

predict xb_validation2_m0 if validation2==1, xb 
	replace xb_validation2_m0=. if validation2!=1
predict pr_validation2_m0 if validation2==1	
brier mort90_postdc pr_validation2_m0 

* standardized mortality ratio
total mort90_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m0 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm0=e(b)[1,1]
	di `prvalidation2deathm0'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm0'	
	
* excess deaths 
di `validation2death' - `prvalidation2deathm0'	
di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort90_postdc male age_sp* `sirs' `aod' `comorbid' if derivation==1 || sta6a:  
estimates store derivation_m1
predict pr_derivation_m1 if e(sample)	
predict re_derivation_m1 if e(sample), reffects reses(rese_derivation_m1)
predict xb_derivation_m1 if e(sample), xb
	replace xb_derivation_m1=. if derivation!=1

egen meanxb_derivation_m1 = mean(xb_derivation_m1) if e(sample)
gen hosplogodd_derivation_m1 = meanxb_derivation_m1 + re_derivation_m1
gen invlog_derivation_m1 = invlogit(hosplogodd_derivation_m1)

brier mort90_postdc pr_derivation_m1 

// standardized mortality ratio
total mort90_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m1 if derivation==1 
	matrix list e(b)
	local prderivationdeathm1=e(b)[1,1]
	di `prderivationdeathm1'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm1'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm1'	
di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000


preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m1 = pr_derivation_m1, nq(`ngroups')
	collapse (sum) mort90_postdc pr_derivation_m1 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m1)
	list numhosps mort90_postdc pr_derivation_m1, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m1
predict pr_validation1_m1 if validation1==1	
predict re_validation1_m1 if validation1==1, reffects reses(rese_validation1_m1)
predict xb_validation1_m1 if validation1==1, xb
	replace xb_validation1_m1=. if validation1!=1
predict stdp_validation1_m1 if validation1==1, stdp 
	replace stdp_validation1_m1=. if validation1!=1

egen meanxb_validation1_m1 = mean(xb_validation1_m1) if validation1==1
gen hosplogodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1
gen invlog_validation1_m1 = invlogit(hosplogodd_validation1_m1)

gen hi_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 + 1.4*stdp_validation1_m1
gen lo_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 - 1.4*stdp_validation1_m1

gen hi_invlog_validation1_m1 = invlogit(hi_logodd_validation1_m1)
gen lo_invlog_validation1_m1 = invlogit(lo_logodd_validation1_m1)

brier mort90_postdc pr_validation1_m1 


// standardized mortality ratio
total mort90_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m1 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm1=e(b)[1,1]
	di `prvalidation1deathm1'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm1'	

// excess deaths 
di `validation1death' - `prvalidation1deathm1'	
di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m1
predict pr_validation2_m1 if validation2==1	
predict re_validation2_m1 if validation2==1, reffects reses(rese_validation2_m1)
predict xb_validation2_m1 if validation2==1, xb
	replace xb_validation2_m1=. if validation2!=1
predict stdp_validation2_m1 if validation2==1, stdp 
	replace stdp_validation2_m1=. if validation2!=1

egen meanxb_validation2_m1 = mean(xb_validation2_m1) if validation2==1
gen hosplogodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1
gen invlog_validation2_m1 = invlogit(hosplogodd_validation2_m1)

gen hi_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 + 1.4*stdp_validation2_m1
gen lo_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 - 1.4*stdp_validation2_m1

gen hi_invlog_validation2_m1 = invlogit(hi_logodd_validation2_m1)
gen lo_invlog_validation2_m1 = invlogit(lo_logodd_validation2_m1)

brier mort90_postdc pr_validation2_m1

// standardized mortality ratio
total mort90_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m1 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm1=e(b)[1,1]
	di `prvalidation2deathm1'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm1'	

// excess deaths 
di `validation2death' - `prvalidation2deathm1'	
di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*-------------------------------------------------------------------------------
* M2: Base + Individual (race, marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort90_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual'		///
		 if derivation==1 || sta6a: 
estimates store derivation_m2
predict pr_derivation_m2 if e(sample)	
predict re_derivation_m2 if e(sample), reffects reses(rese_derivation_m2)
predict xb_derivation_m2 if e(sample), xb
	replace xb_derivation_m2=. if derivation!=1

egen meanxb_derivation_m2 = mean(xb_derivation_m2) if e(sample)
gen hosplogodd_derivation_m2 = meanxb_derivation_m2 + re_derivation_m2
gen invlog_derivation_m2 = invlogit(hosplogodd_derivation_m2)

brier mort90_postdc pr_derivation_m2 

// standardized mortality ratio
total mort90_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m2 if derivation==1 
	matrix list e(b)
	local prderivationdeathm2=e(b)[1,1]
	di `prderivationdeathm2'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm2'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm2'	
di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m2
predict pr_validation1_m2 if validation1==1	
predict re_validation1_m2 if validation1==1, reffects reses(rese_validation1_m2)
predict xb_validation1_m2 if validation1==1, xb
	replace xb_validation1_m2=. if validation1!=1
predict stdp_validation1_m2 if validation1==1, stdp 
	replace stdp_validation1_m2=. if validation1!=1

egen meanxb_validation1_m2 = mean(xb_validation1_m2) if validation1==1
gen hosplogodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2
gen invlog_validation1_m2 = invlogit(hosplogodd_validation1_m2)

gen hi_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 + 1.4*stdp_validation1_m2
gen lo_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 - 1.4*stdp_validation1_m2

gen hi_invlog_validation1_m2 = invlogit(hi_logodd_validation1_m2)
gen lo_invlog_validation1_m2 = invlogit(lo_logodd_validation1_m2)

brier mort90_postdc pr_validation1_m2


// standardized mortality ratio
total mort90_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m2 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm2=e(b)[1,1]
	di `prvalidation1deathm2'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm2'	

// excess deaths 
di `validation1death' - `prvalidation1deathm2'	
di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m2
predict pr_validation2_m2 if validation2==1	
predict re_validation2_m2 if validation2==1, reffects reses(rese_validation2_m2)
predict xb_validation2_m2 if validation2==1, xb
	replace xb_validation2_m2=. if validation2!=1
predict stdp_validation2_m2 if validation2==1, stdp 
	replace stdp_validation2_m2=. if validation2!=1

egen meanxb_validation2_m2 = mean(xb_validation2_m2) if validation2==1
gen hosplogodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2
gen invlog_validation2_m2 = invlogit(hosplogodd_validation2_m2)

gen hi_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 + 1.4*stdp_validation2_m2
gen lo_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 - 1.4*stdp_validation2_m2

gen hi_invlog_validation2_m2 = invlogit(hi_logodd_validation2_m2)
gen lo_invlog_validation2_m2 = invlogit(lo_logodd_validation2_m2)

brier mort90_postdc pr_validation2_m2

// standardized mortality ratio
total mort90_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m2 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm2=e(b)[1,1]
	di `prvalidation2deathm2'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm2'	

// excess deaths 
di `validation2death' - `prvalidation2deathm2'	
di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000
	
	
*----------------------------
* M3: Base + Community 
*----------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort90_postdc male age_sp* `sirs' `aod' `comorbid' `community'		///	
		if derivation==1 || sta6a: 
estimates store derivation_m3
predict pr_derivation_m3 if e(sample)	
predict re_derivation_m3 if e(sample), reffects reses(rese_derivation_m3)
predict xb_derivation_m3 if e(sample), xb
	replace xb_derivation_m3=. if derivation!=1

egen meanxb_derivation_m3 = mean(xb_derivation_m3) if e(sample)
gen hosplogodd_derivation_m3 = meanxb_derivation_m3 + re_derivation_m3
gen invlog_derivation_m3 = invlogit(hosplogodd_derivation_m3)

brier mort90_postdc pr_derivation_m3 

// standardized mortality ratio
total mort90_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m3 if derivation==1 
	matrix list e(b)
	local prderivationdeathm3=e(b)[1,1]
	di `prderivationdeathm3'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm3'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm3'	
di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m3
predict pr_validation1_m3 if validation1==1	
predict re_validation1_m3 if validation1==1, reffects reses(rese_validation1_m3)
predict xb_validation1_m3 if validation1==1, xb
	replace xb_validation1_m3=. if validation1!=1
predict stdp_validation1_m3 if validation1==1, stdp 
	replace stdp_validation1_m3=. if validation1!=1

egen meanxb_validation1_m3 = mean(xb_validation1_m3) if validation1==1
gen hosplogodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3
gen invlog_validation1_m3 = invlogit(hosplogodd_validation1_m3)

gen hi_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 + 1.4*stdp_validation1_m3
gen lo_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 - 1.4*stdp_validation1_m3

gen hi_invlog_validation1_m3 = invlogit(hi_logodd_validation1_m3)
gen lo_invlog_validation1_m3 = invlogit(lo_logodd_validation1_m3)

brier mort90_postdc pr_validation1_m3


// standardized mortality ratio
total mort90_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m3 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm3=e(b)[1,1]
	di `prvalidation1deathm3'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm3'	

// excess deaths 
di `validation1death' - `prvalidation1deathm3'	
di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m3
predict pr_validation2_m3 if validation2==1	
predict re_validation2_m3 if validation2==1, reffects reses(rese_validation2_m3)
predict xb_validation2_m3 if validation2==1, xb
	replace xb_validation2_m3=. if validation2!=1
predict stdp_validation2_m3 if validation2==1, stdp 
	replace stdp_validation2_m3=. if validation2!=1

egen meanxb_validation2_m3 = mean(xb_validation2_m3) if validation2==1
gen hosplogodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3
gen invlog_validation2_m3 = invlogit(hosplogodd_validation2_m3)

gen hi_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 + 1.4*stdp_validation2_m3
gen lo_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 - 1.4*stdp_validation2_m3

gen hi_invlog_validation2_m3 = invlogit(hi_logodd_validation2_m3)
gen lo_invlog_validation2_m3 = invlogit(lo_logodd_validation2_m3)

brier mort90_postdc pr_validation2_m3

// standardized mortality ratio
total mort90_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m3 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm3=e(b)[1,1]
	di `prvalidation2deathm3'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm3'	

// excess deaths 
di `validation2death' - `prvalidation2deathm3'	
di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort90_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' ///
		if derivation==1 || sta6a: 
estimates store derivation_m4
predict pr_derivation_m4 if e(sample)	
predict re_derivation_m4 if e(sample), reffects reses(rese_derivation_m4)
predict xb_derivation_m4 if e(sample), xb
	replace xb_derivation_m4=. if derivation!=1

egen meanxb_derivation_m4 = mean(xb_derivation_m4) if e(sample)
gen hosplogodd_derivation_m4 = meanxb_derivation_m4 + re_derivation_m4
gen invlog_derivation_m4 = invlogit(hosplogodd_derivation_m4)

brier mort90_postdc pr_derivation_m4 

// standardized mortality ratio
total mort90_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m4 if derivation==1 
	matrix list e(b)
	local prderivationdeathm4=e(b)[1,1]
	di `prderivationdeathm4'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm4'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm4'	
di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m4 = pr_derivation_m4, nq(`ngroups')
	bysort quantile_derivation_m4: sum pr_derivation_m4, de
	collapse (sum) mort90_postdc pr_derivation_m4 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m4)
	list numhosps mort90_postdc pr_derivation_m4, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m4
predict pr_validation1_m4 if validation1==1	
predict re_validation1_m4 if validation1==1, reffects reses(rese_validation1_m4)
predict xb_validation1_m4 if validation1==1, xb
	replace xb_validation1_m4=. if validation1!=1
predict stdp_validation1_m4 if validation1==1, stdp 
	replace stdp_validation1_m4=. if validation1!=1

egen meanxb_validation1_m4 = mean(xb_validation1_m4) if validation1==1
gen hosplogodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4
gen invlog_validation1_m4 = invlogit(hosplogodd_validation1_m4)

gen hi_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 + 1.4*stdp_validation1_m4
gen lo_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 - 1.4*stdp_validation1_m4

gen hi_invlog_validation1_m4 = invlogit(hi_logodd_validation1_m4)
gen lo_invlog_validation1_m4 = invlogit(lo_logodd_validation1_m4)

brier mort90_postdc pr_validation1_m4


// standardized mortality ratio
total mort90_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m4 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm4=e(b)[1,1]
	di `prvalidation1deathm4'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm4'	

// excess deaths 
di `validation1death' - `prvalidation1deathm4'	
di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000

preserve 
	keep if validation1==1 
	local ngroups 10
	xtile quantile_validation1_m4 = pr_validation1_m4, nq(`ngroups')
	bysort quantile_validation1_m4: sum pr_validation1_m4, de
	collapse (sum) mort90_postdc pr_validation1_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation1_m4)
	list numhosps mort90_postdc pr_validation1_m4, abb(20)
restore


* VALIDATION 2021 * 
estimates restore derivation_m4
predict pr_validation2_m4 if validation2==1	
predict re_validation2_m4 if validation2==1, reffects reses(rese_validation2_m4)
predict xb_validation2_m4 if validation2==1, xb
	replace xb_validation2_m4=. if validation2!=1
predict stdp_validation2_m4 if validation2==1, stdp 
	replace stdp_validation2_m4=. if validation2!=1

egen meanxb_validation2_m4 = mean(xb_validation2_m4) if validation2==1
gen hosplogodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4
gen invlog_validation2_m4 = invlogit(hosplogodd_validation2_m4)

gen hi_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 + 1.4*stdp_validation2_m4
gen lo_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 - 1.4*stdp_validation2_m4

gen hi_invlog_validation2_m4 = invlogit(hi_logodd_validation2_m4)
gen lo_invlog_validation2_m4 = invlogit(lo_logodd_validation2_m4)

brier mort90_postdc pr_validation2_m4

// standardized mortality ratio
total mort90_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m4 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm4=e(b)[1,1]
	di `prvalidation2deathm4'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm4'	

// excess deaths 
di `validation2death' - `prvalidation2deathm4'	
di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


preserve 
	keep if validation2==1 
	local ngroups 10
	xtile quantile_validation2_m4 = pr_validation2_m4, nq(`ngroups')
	bysort quantile_validation2_m4: sum pr_validation2_m4, de
	collapse (sum) mort90_postdc pr_validation2_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation2_m4)
	list numhosps mort90_postdc pr_validation2_m4, abb(20)
restore


*-----------------------------------------
* Comparing models 0 & 1 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_postdc xb_derivation_m0 xb_derivation_m1 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_postdc xb_validation1_m0 xb_validation1_m1 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_postdc xb_validation2_m0 xb_validation2_m1 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 2 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_postdc xb_derivation_m1 xb_derivation_m2 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_postdc xb_validation1_m1 xb_validation1_m2 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_postdc xb_validation2_m1 xb_validation2_m2 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 3 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_postdc xb_derivation_m1 xb_derivation_m3 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_postdc xb_validation1_m1 xb_validation1_m3 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_postdc xb_validation2_m1 xb_validation2_m3 	///
					if validation2==1
					
*-----------------------------------------
* Comparing models 1 & 4 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort90_postdc xb_derivation_m1 xb_derivation_m4 	///
					if derivation==1

* Internal validation data set 
roccomp mort90_postdc xb_validation1_m1 xb_validation1_m4 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort90_postdc xb_validation2_m1 xb_validation2_m4 	///
					if validation2==1
			

*-------------------------------------------------------------------------------
* Supplemental Table: Comparison of model performance - 360-day Post-Discharge
*-------------------------------------------------------------------------------

drop pr_derivation_m0-lo_invlog_validation2_m4

*----------------------
* M0: Community Only
*----------------------

* DERIVATION *
melogit mort360_postdc adi_impute if derivation==1 || sta6a: 
estimates store derivation_m0

predict pr_derivation_m0 if e(sample)
predict xb_derivation_m0 if e(sample), xb
brier mort360_postdc pr_derivation_m0 

* standardized mortality ratio
total mort360_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m0 if derivation==1 
	matrix list e(b)
	local prderivationdeathm0=e(b)[1,1]
	di `prderivationdeathm0'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm0'	
	
* excess deaths 
di `derivationdeath' - `prderivationdeathm0'	
di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m0

predict xb_validation1_m0 if validation1==1, xb 
	replace xb_validation1_m0=. if validation1!=1
predict pr_validation1_m0 if validation1==1	
predict re_validation1_m0 if validation1==1, reffects reses(rese_validation1_m0)
predict stdp_validation1_m0 if validation1==1, stdp 
	replace stdp_validation1_m0=. if validation1!=1

egen meanxb_validation1_m0 = mean(xb_validation1_m0) if validation1==1
gen hosplogodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0
gen invlog_validation1_m0 = invlogit(hosplogodd_validation1_m0)

gen hi_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 + 1.4*stdp_validation1_m0
gen lo_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 - 1.4*stdp_validation1_m0

gen hi_invlog_validation1_m0 = invlogit(hi_logodd_validation1_m0)
gen lo_invlog_validation1_m0 = invlogit(lo_logodd_validation1_m0)

brier mort360_postdc pr_validation1_m0 


* standardized mortality ratio
total mort360_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m0 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm0=e(b)[1,1]
	di `prvalidation1deathm0'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm0'	
	
* excess deaths 
di `validation1death' - `prvalidation1deathm0'	
di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000


* VALIDATION *
estimates restore derivation_m0

predict xb_validation2_m0 if validation2==1, xb 
	replace xb_validation2_m0=. if validation2!=1
predict pr_validation2_m0 if validation2==1	
brier mort360_postdc pr_validation2_m0 

* standardized mortality ratio
total mort360_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m0 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm0=e(b)[1,1]
	di `prvalidation2deathm0'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm0'	
	
* excess deaths 
di `validation2death' - `prvalidation2deathm0'	
di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort360_postdc male age_sp* `sirs' `aod' `comorbid' if derivation==1 || sta6a:  
estimates store derivation_m1
predict pr_derivation_m1 if e(sample)	
predict re_derivation_m1 if e(sample), reffects reses(rese_derivation_m1)
predict xb_derivation_m1 if e(sample), xb
	replace xb_derivation_m1=. if derivation!=1

egen meanxb_derivation_m1 = mean(xb_derivation_m1) if e(sample)
gen hosplogodd_derivation_m1 = meanxb_derivation_m1 + re_derivation_m1
gen invlog_derivation_m1 = invlogit(hosplogodd_derivation_m1)

brier mort360_postdc pr_derivation_m1 

// standardized mortality ratio
total mort360_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m1 if derivation==1 
	matrix list e(b)
	local prderivationdeathm1=e(b)[1,1]
	di `prderivationdeathm1'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm1'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm1'	
di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000


preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m1 = pr_derivation_m1, nq(`ngroups')
	collapse (sum) mort360_postdc pr_derivation_m1 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m1)
	list numhosps mort360_postdc pr_derivation_m1, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m1
predict pr_validation1_m1 if validation1==1	
predict re_validation1_m1 if validation1==1, reffects reses(rese_validation1_m1)
predict xb_validation1_m1 if validation1==1, xb
	replace xb_validation1_m1=. if validation1!=1
predict stdp_validation1_m1 if validation1==1, stdp 
	replace stdp_validation1_m1=. if validation1!=1

egen meanxb_validation1_m1 = mean(xb_validation1_m1) if validation1==1
gen hosplogodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1
gen invlog_validation1_m1 = invlogit(hosplogodd_validation1_m1)

gen hi_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 + 1.4*stdp_validation1_m1
gen lo_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 - 1.4*stdp_validation1_m1

gen hi_invlog_validation1_m1 = invlogit(hi_logodd_validation1_m1)
gen lo_invlog_validation1_m1 = invlogit(lo_logodd_validation1_m1)

brier mort360_postdc pr_validation1_m1


// standardized mortality ratio
total mort360_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m1 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm1=e(b)[1,1]
	di `prvalidation1deathm1'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm1'	

// excess deaths 
di `validation1death' - `prvalidation1deathm1'	
di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m1
predict pr_validation2_m1 if validation2==1	
predict re_validation2_m1 if validation2==1, reffects reses(rese_validation2_m1)
predict xb_validation2_m1 if validation2==1, xb
	replace xb_validation2_m1=. if validation2!=1
predict stdp_validation2_m1 if validation2==1, stdp 
	replace stdp_validation2_m1=. if validation2!=1

egen meanxb_validation2_m1 = mean(xb_validation2_m1) if validation2==1
gen hosplogodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1
gen invlog_validation2_m1 = invlogit(hosplogodd_validation2_m1)

gen hi_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 + 1.4*stdp_validation2_m1
gen lo_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 - 1.4*stdp_validation2_m1

gen hi_invlog_validation2_m1 = invlogit(hi_logodd_validation2_m1)
gen lo_invlog_validation2_m1 = invlogit(lo_logodd_validation2_m1)

brier mort360_postdc pr_validation2_m1

// standardized mortality ratio
total mort360_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m1 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm1=e(b)[1,1]
	di `prvalidation2deathm1'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm1'	

// excess deaths 
di `validation2death' - `prvalidation2deathm1'	
di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*-------------------------------------------------------------------------------
* M2: Base + Individual (race, marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort360_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual'		///
		 if derivation==1 || sta6a: 
estimates store derivation_m2
predict pr_derivation_m2 if e(sample)	
predict re_derivation_m2 if e(sample), reffects reses(rese_derivation_m2)
predict xb_derivation_m2 if e(sample), xb
	replace xb_derivation_m2=. if derivation!=1

egen meanxb_derivation_m2 = mean(xb_derivation_m2) if e(sample)
gen hosplogodd_derivation_m2 = meanxb_derivation_m2 + re_derivation_m2
gen invlog_derivation_m2 = invlogit(hosplogodd_derivation_m2)

brier mort360_postdc pr_derivation_m2 

// standardized mortality ratio
total mort360_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m2 if derivation==1 
	matrix list e(b)
	local prderivationdeathm2=e(b)[1,1]
	di `prderivationdeathm2'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm2'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm2'	
di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m2
predict pr_validation1_m2 if validation1==1	
predict re_validation1_m2 if validation1==1, reffects reses(rese_validation1_m2)
predict xb_validation1_m2 if validation1==1, xb
	replace xb_validation1_m2=. if validation1!=1
predict stdp_validation1_m2 if validation1==1, stdp 
	replace stdp_validation1_m2=. if validation1!=1

egen meanxb_validation1_m2 = mean(xb_validation1_m2) if validation1==1
gen hosplogodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2
gen invlog_validation1_m2 = invlogit(hosplogodd_validation1_m2)

gen hi_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 + 1.4*stdp_validation1_m2
gen lo_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 - 1.4*stdp_validation1_m2

gen hi_invlog_validation1_m2 = invlogit(hi_logodd_validation1_m2)
gen lo_invlog_validation1_m2 = invlogit(lo_logodd_validation1_m2)

brier mort360_postdc pr_validation1_m2


// standardized mortality ratio
total mort360_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m2 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm2=e(b)[1,1]
	di `prvalidation1deathm2'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm2'	

// excess deaths 
di `validation1death' - `prvalidation1deathm2'	
di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m2
predict pr_validation2_m2 if validation2==1	
predict re_validation2_m2 if validation2==1, reffects reses(rese_validation2_m2)
predict xb_validation2_m2 if validation2==1, xb
	replace xb_validation2_m2=. if validation2!=1
predict stdp_validation2_m2 if validation2==1, stdp 
	replace stdp_validation2_m2=. if validation2!=1

egen meanxb_validation2_m2 = mean(xb_validation2_m2) if validation2==1
gen hosplogodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2
gen invlog_validation2_m2 = invlogit(hosplogodd_validation2_m2)

gen hi_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 + 1.4*stdp_validation2_m2
gen lo_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 - 1.4*stdp_validation2_m2

gen hi_invlog_validation2_m2 = invlogit(hi_logodd_validation2_m2)
gen lo_invlog_validation2_m2 = invlogit(lo_logodd_validation2_m2)

brier mort360_postdc pr_validation2_m2

// standardized mortality ratio
total mort360_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m2 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm2=e(b)[1,1]
	di `prvalidation2deathm2'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm2'	

// excess deaths 
di `validation2death' - `prvalidation2deathm2'	
di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000
	
	
*----------------------------
* M3: Base + Community 
*----------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort360_postdc male age_sp* `sirs' `aod' `comorbid' `community'		///	
		if derivation==1 || sta6a: 
estimates store derivation_m3
predict pr_derivation_m3 if e(sample)	
predict re_derivation_m3 if e(sample), reffects reses(rese_derivation_m3)
predict xb_derivation_m3 if e(sample), xb
	replace xb_derivation_m3=. if derivation!=1

egen meanxb_derivation_m3 = mean(xb_derivation_m3) if e(sample)
gen hosplogodd_derivation_m3 = meanxb_derivation_m3 + re_derivation_m3
gen invlog_derivation_m3 = invlogit(hosplogodd_derivation_m3)

brier mort360_postdc pr_derivation_m3 

// standardized mortality ratio
total mort360_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m3 if derivation==1 
	matrix list e(b)
	local prderivationdeathm3=e(b)[1,1]
	di `prderivationdeathm3'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm3'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm3'	
di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m3
predict pr_validation1_m3 if validation1==1	
predict re_validation1_m3 if validation1==1, reffects reses(rese_validation1_m3)
predict xb_validation1_m3 if validation1==1, xb
	replace xb_validation1_m3=. if validation1!=1
predict stdp_validation1_m3 if validation1==1, stdp 
	replace stdp_validation1_m3=. if validation1!=1

egen meanxb_validation1_m3 = mean(xb_validation1_m3) if validation1==1
gen hosplogodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3
gen invlog_validation1_m3 = invlogit(hosplogodd_validation1_m3)

gen hi_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 + 1.4*stdp_validation1_m3
gen lo_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 - 1.4*stdp_validation1_m3

gen hi_invlog_validation1_m3 = invlogit(hi_logodd_validation1_m3)
gen lo_invlog_validation1_m3 = invlogit(lo_logodd_validation1_m3)

brier mort360_postdc pr_validation1_m3


// standardized mortality ratio
total mort360_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m3 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm3=e(b)[1,1]
	di `prvalidation1deathm3'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm3'	

// excess deaths 
di `validation1death' - `prvalidation1deathm3'	
di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m3
predict pr_validation2_m3 if validation2==1	
predict re_validation2_m3 if validation2==1, reffects reses(rese_validation2_m3)
predict xb_validation2_m3 if validation2==1, xb
	replace xb_validation2_m3=. if validation2!=1
predict stdp_validation2_m3 if validation2==1, stdp 
	replace stdp_validation2_m3=. if validation2!=1

egen meanxb_validation2_m3 = mean(xb_validation2_m3) if validation2==1
gen hosplogodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3
gen invlog_validation2_m3 = invlogit(hosplogodd_validation2_m3)

gen hi_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 + 1.4*stdp_validation2_m3
gen lo_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 - 1.4*stdp_validation2_m3

gen hi_invlog_validation2_m3 = invlogit(hi_logodd_validation2_m3)
gen lo_invlog_validation2_m3 = invlogit(lo_logodd_validation2_m3)

brier mort360_postdc pr_validation2_m3

// standardized mortality ratio
total mort360_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m3 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm3=e(b)[1,1]
	di `prvalidation2deathm3'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm3'	

// excess deaths 
di `validation2death' - `prvalidation2deathm3'	
di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort360_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' ///
		if derivation==1 || sta6a: 
estimates store derivation_m4
predict pr_derivation_m4 if e(sample)	
predict re_derivation_m4 if e(sample), reffects reses(rese_derivation_m4)
predict xb_derivation_m4 if e(sample), xb
	replace xb_derivation_m4=. if derivation!=1

egen meanxb_derivation_m4 = mean(xb_derivation_m4) if e(sample)
gen hosplogodd_derivation_m4 = meanxb_derivation_m4 + re_derivation_m4
gen invlog_derivation_m4 = invlogit(hosplogodd_derivation_m4)

brier mort360_postdc pr_derivation_m4 

// standardized mortality ratio
total mort360_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m4 if derivation==1 
	matrix list e(b)
	local prderivationdeathm4=e(b)[1,1]
	di `prderivationdeathm4'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm4'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm4'	
di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m4 = pr_derivation_m4, nq(`ngroups')
	bysort quantile_derivation_m4: sum pr_derivation_m4, de
	collapse (sum) mort360_postdc pr_derivation_m4 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m4)
	list numhosps mort360_postdc pr_derivation_m4, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m4
predict pr_validation1_m4 if validation1==1	
predict re_validation1_m4 if validation1==1, reffects reses(rese_validation1_m4)
predict xb_validation1_m4 if validation1==1, xb
	replace xb_validation1_m4=. if validation1!=1
predict stdp_validation1_m4 if validation1==1, stdp 
	replace stdp_validation1_m4=. if validation1!=1

egen meanxb_validation1_m4 = mean(xb_validation1_m4) if validation1==1
gen hosplogodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4
gen invlog_validation1_m4 = invlogit(hosplogodd_validation1_m4)

gen hi_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 + 1.4*stdp_validation1_m4
gen lo_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 - 1.4*stdp_validation1_m4

gen hi_invlog_validation1_m4 = invlogit(hi_logodd_validation1_m4)
gen lo_invlog_validation1_m4 = invlogit(lo_logodd_validation1_m4)

brier mort360_postdc pr_validation1_m4


// standardized mortality ratio
total mort360_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m4 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm4=e(b)[1,1]
	di `prvalidation1deathm4'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm4'	

// excess deaths 
di `validation1death' - `prvalidation1deathm4'	
di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000

preserve 
	keep if validation1==1 
	local ngroups 10
	xtile quantile_validation1_m4 = pr_validation1_m4, nq(`ngroups')
	bysort quantile_validation1_m4: sum pr_validation1_m4, de
	collapse (sum) mort360_postdc pr_validation1_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation1_m4)
	list numhosps mort360_postdc pr_validation1_m4, abb(20)
restore


* VALIDATION 2021 * 
estimates restore derivation_m4
predict pr_validation2_m4 if validation2==1	
predict re_validation2_m4 if validation2==1, reffects reses(rese_validation2_m4)
predict xb_validation2_m4 if validation2==1, xb
	replace xb_validation2_m4=. if validation2!=1
predict stdp_validation2_m4 if validation2==1, stdp 
	replace stdp_validation2_m4=. if validation2!=1

egen meanxb_validation2_m4 = mean(xb_validation2_m4) if validation2==1
gen hosplogodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4
gen invlog_validation2_m4 = invlogit(hosplogodd_validation2_m4)

gen hi_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 + 1.4*stdp_validation2_m4
gen lo_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 - 1.4*stdp_validation2_m4

gen hi_invlog_validation2_m4 = invlogit(hi_logodd_validation2_m4)
gen lo_invlog_validation2_m4 = invlogit(lo_logodd_validation2_m4)

brier mort360_postdc pr_validation2_m4

// standardized mortality ratio
total mort360_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m4 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm4=e(b)[1,1]
	di `prvalidation2deathm4'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm4'	

// excess deaths 
di `validation2death' - `prvalidation2deathm4'	
di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


preserve 
	keep if validation2==1 
	local ngroups 10
	xtile quantile_validation2_m4 = pr_validation2_m4, nq(`ngroups')
	bysort quantile_validation2_m4: sum pr_validation2_m4, de
	collapse (sum) mort360_postdc pr_validation2_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation2_m4)
	list numhosps mort360_postdc pr_validation2_m4, abb(20)
restore


*-----------------------------------------
* Comparing models 0 & 1 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort360_postdc xb_derivation_m0 xb_derivation_m1 	///
					if derivation==1

* Internal validation data set 
roccomp mort360_postdc xb_validation1_m0 xb_validation1_m1 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort360_postdc xb_validation2_m0 xb_validation2_m1 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 2 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort360_postdc xb_derivation_m1 xb_derivation_m2 	///
					if derivation==1

* Internal validation data set 
roccomp mort360_postdc xb_validation1_m1 xb_validation1_m2 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort360_postdc xb_validation2_m1 xb_validation2_m2 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 3 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort360_postdc xb_derivation_m1 xb_derivation_m3 	///
					if derivation==1

* Internal validation data set 
roccomp mort360_postdc xb_validation1_m1 xb_validation1_m3 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort360_postdc xb_validation2_m1 xb_validation2_m3 	///
					if validation2==1
					
*-----------------------------------------
* Comparing models 1 & 4 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort360_postdc xb_derivation_m1 xb_derivation_m4 	///
					if derivation==1

* Internal validation data set 
roccomp mort360_postdc xb_validation1_m1 xb_validation1_m4 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort360_postdc xb_validation2_m1 xb_validation2_m4 	///
					if validation2==1
			

*-------------------------------------------------------------------------------
* Supplemental Table: Comparison of model performance - 360-day Post-Discharge
*-------------------------------------------------------------------------------

drop pr_derivation_m0-lo_invlog_validation2_m4

*----------------------
* M0: Community Only
*----------------------

* DERIVATION *
melogit mort720_postdc adi_impute if derivation==1 || sta6a: 
estimates store derivation_m0

predict pr_derivation_m0 if e(sample)
predict xb_derivation_m0 if e(sample), xb
brier mort720_postdc pr_derivation_m0 

* standardized mortality ratio
total mort720_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m0 if derivation==1 
	matrix list e(b)
	local prderivationdeathm0=e(b)[1,1]
	di `prderivationdeathm0'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm0'	
	
* excess deaths 
di `derivationdeath' - `prderivationdeathm0'	
di ((`derivationdeath' - `prderivationdeathm0')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m0

predict xb_validation1_m0 if validation1==1, xb 
	replace xb_validation1_m0=. if validation1!=1
predict pr_validation1_m0 if validation1==1	
predict re_validation1_m0 if validation1==1, reffects reses(rese_validation1_m0)
predict stdp_validation1_m0 if validation1==1, stdp 
	replace stdp_validation1_m0=. if validation1!=1

egen meanxb_validation1_m0 = mean(xb_validation1_m0) if validation1==1
gen hosplogodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0
gen invlog_validation1_m0 = invlogit(hosplogodd_validation1_m0)

gen hi_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 + 1.4*stdp_validation1_m0
gen lo_logodd_validation1_m0 = meanxb_validation1_m0 + re_validation1_m0 - 1.4*stdp_validation1_m0

gen hi_invlog_validation1_m0 = invlogit(hi_logodd_validation1_m0)
gen lo_invlog_validation1_m0 = invlogit(lo_logodd_validation1_m0)

brier mort720_postdc pr_validation1_m0 


* standardized mortality ratio
total mort720_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m0 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm0=e(b)[1,1]
	di `prvalidation1deathm0'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm0'	
	
* excess deaths 
di `validation1death' - `prvalidation1deathm0'	
di ((`validation1death' - `prvalidation1deathm0')/`validation1')*10000


* VALIDATION *
estimates restore derivation_m0

predict xb_validation2_m0 if validation2==1, xb 
	replace xb_validation2_m0=. if validation2!=1
predict pr_validation2_m0 if validation2==1	
brier mort720_postdc pr_validation2_m0 

* standardized mortality ratio
total mort720_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m0 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm0=e(b)[1,1]
	di `prvalidation2deathm0'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm0'	
	
* excess deaths 
di `validation2death' - `prvalidation2deathm0'	
di ((`validation2death' - `prvalidation2deathm0')/`validation2')*10000


*-------------------------------------------------------------------------------
* M1: Base (age, gender, comorbidities, acute organ dysfunctions, SIRS criteria)
*-------------------------------------------------------------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc
			
melogit mort720_postdc male age_sp* `sirs' `aod' `comorbid' if derivation==1 || sta6a:  
estimates store derivation_m1
predict pr_derivation_m1 if e(sample)	
predict re_derivation_m1 if e(sample), reffects reses(rese_derivation_m1)
predict xb_derivation_m1 if e(sample), xb
	replace xb_derivation_m1=. if derivation!=1

egen meanxb_derivation_m1 = mean(xb_derivation_m1) if e(sample)
gen hosplogodd_derivation_m1 = meanxb_derivation_m1 + re_derivation_m1
gen invlog_derivation_m1 = invlogit(hosplogodd_derivation_m1)

brier mort720_postdc pr_derivation_m1 

// standardized mortality ratio
total mort720_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m1 if derivation==1 
	matrix list e(b)
	local prderivationdeathm1=e(b)[1,1]
	di `prderivationdeathm1'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm1'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm1'	
di ((`derivationdeath' - `prderivationdeathm1')/`derivation')*10000


preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m1 = pr_derivation_m1, nq(`ngroups')
	collapse (sum) mort720_postdc pr_derivation_m1 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m1)
	list numhosps mort720_postdc pr_derivation_m1, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m1
predict pr_validation1_m1 if validation1==1	
predict re_validation1_m1 if validation1==1, reffects reses(rese_validation1_m1)
predict xb_validation1_m1 if validation1==1, xb
	replace xb_validation1_m1=. if validation1!=1
predict stdp_validation1_m1 if validation1==1, stdp 
	replace stdp_validation1_m1=. if validation1!=1

egen meanxb_validation1_m1 = mean(xb_validation1_m1) if validation1==1
gen hosplogodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1
gen invlog_validation1_m1 = invlogit(hosplogodd_validation1_m1)

gen hi_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 + 1.4*stdp_validation1_m1
gen lo_logodd_validation1_m1 = meanxb_validation1_m1 + re_validation1_m1 - 1.4*stdp_validation1_m1

gen hi_invlog_validation1_m1 = invlogit(hi_logodd_validation1_m1)
gen lo_invlog_validation1_m1 = invlogit(lo_logodd_validation1_m1)

brier mort720_postdc pr_validation1_m1


// standardized mortality ratio
total mort720_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m1 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm1=e(b)[1,1]
	di `prvalidation1deathm1'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm1'	

// excess deaths 
di `validation1death' - `prvalidation1deathm1'	
di ((`validation1death' - `prvalidation1deathm1')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m1
predict pr_validation2_m1 if validation2==1	
predict re_validation2_m1 if validation2==1, reffects reses(rese_validation2_m1)
predict xb_validation2_m1 if validation2==1, xb
	replace xb_validation2_m1=. if validation2!=1
predict stdp_validation2_m1 if validation2==1, stdp 
	replace stdp_validation2_m1=. if validation2!=1

egen meanxb_validation2_m1 = mean(xb_validation2_m1) if validation2==1
gen hosplogodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1
gen invlog_validation2_m1 = invlogit(hosplogodd_validation2_m1)

gen hi_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 + 1.4*stdp_validation2_m1
gen lo_logodd_validation2_m1 = meanxb_validation2_m1 + re_validation2_m1 - 1.4*stdp_validation2_m1

gen hi_invlog_validation2_m1 = invlogit(hi_logodd_validation2_m1)
gen lo_invlog_validation2_m1 = invlogit(lo_logodd_validation2_m1)

brier mort720_postdc pr_validation2_m1

// standardized mortality ratio
total mort720_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m1 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm1=e(b)[1,1]
	di `prvalidation2deathm1'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm1'	

// excess deaths 
di `validation2death' - `prvalidation2deathm1'	
di ((`validation2death' - `prvalidation2deathm1')/`validation2')*10000


*-------------------------------------------------------------------------------
* M2: Base + Individual (race, marital, housing instability, frailty, CAN) 
*-------------------------------------------------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute			
			
melogit mort720_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual'		///
		 if derivation==1 || sta6a: 
estimates store derivation_m2
predict pr_derivation_m2 if e(sample)	
predict re_derivation_m2 if e(sample), reffects reses(rese_derivation_m2)
predict xb_derivation_m2 if e(sample), xb
	replace xb_derivation_m2=. if derivation!=1

egen meanxb_derivation_m2 = mean(xb_derivation_m2) if e(sample)
gen hosplogodd_derivation_m2 = meanxb_derivation_m2 + re_derivation_m2
gen invlog_derivation_m2 = invlogit(hosplogodd_derivation_m2)

brier mort720_postdc pr_derivation_m2 

// standardized mortality ratio
total mort720_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m2 if derivation==1 
	matrix list e(b)
	local prderivationdeathm2=e(b)[1,1]
	di `prderivationdeathm2'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm2'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm2'	
di ((`derivationdeath' - `prderivationdeathm2')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m2
predict pr_validation1_m2 if validation1==1	
predict re_validation1_m2 if validation1==1, reffects reses(rese_validation1_m2)
predict xb_validation1_m2 if validation1==1, xb
	replace xb_validation1_m2=. if validation1!=1
predict stdp_validation1_m2 if validation1==1, stdp 
	replace stdp_validation1_m2=. if validation1!=1

egen meanxb_validation1_m2 = mean(xb_validation1_m2) if validation1==1
gen hosplogodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2
gen invlog_validation1_m2 = invlogit(hosplogodd_validation1_m2)

gen hi_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 + 1.4*stdp_validation1_m2
gen lo_logodd_validation1_m2 = meanxb_validation1_m2 + re_validation1_m2 - 1.4*stdp_validation1_m2

gen hi_invlog_validation1_m2 = invlogit(hi_logodd_validation1_m2)
gen lo_invlog_validation1_m2 = invlogit(lo_logodd_validation1_m2)

brier mort720_postdc pr_validation1_m2


// standardized mortality ratio
total mort720_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m2 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm2=e(b)[1,1]
	di `prvalidation1deathm2'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm2'	

// excess deaths 
di `validation1death' - `prvalidation1deathm2'	
di ((`validation1death' - `prvalidation1deathm2')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m2
predict pr_validation2_m2 if validation2==1	
predict re_validation2_m2 if validation2==1, reffects reses(rese_validation2_m2)
predict xb_validation2_m2 if validation2==1, xb
	replace xb_validation2_m2=. if validation2!=1
predict stdp_validation2_m2 if validation2==1, stdp 
	replace stdp_validation2_m2=. if validation2!=1

egen meanxb_validation2_m2 = mean(xb_validation2_m2) if validation2==1
gen hosplogodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2
gen invlog_validation2_m2 = invlogit(hosplogodd_validation2_m2)

gen hi_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 + 1.4*stdp_validation2_m2
gen lo_logodd_validation2_m2 = meanxb_validation2_m2 + re_validation2_m2 - 1.4*stdp_validation2_m2

gen hi_invlog_validation2_m2 = invlogit(hi_logodd_validation2_m2)
gen lo_invlog_validation2_m2 = invlogit(lo_logodd_validation2_m2)

brier mort720_postdc pr_validation2_m2

// standardized mortality ratio
total mort720_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m2 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm2=e(b)[1,1]
	di `prvalidation2deathm2'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm2'	

// excess deaths 
di `validation2death' - `prvalidation2deathm2'	
di ((`validation2death' - `prvalidation2deathm2')/`validation2')*10000
	
	
*----------------------------
* M3: Base + Community 
*----------------------------

* DERIVATION *
local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local community 															///	
			rural adi_impute
			
melogit mort720_postdc male age_sp* `sirs' `aod' `comorbid' `community'		///	
		if derivation==1 || sta6a: 
estimates store derivation_m3
predict pr_derivation_m3 if e(sample)	
predict re_derivation_m3 if e(sample), reffects reses(rese_derivation_m3)
predict xb_derivation_m3 if e(sample), xb
	replace xb_derivation_m3=. if derivation!=1

egen meanxb_derivation_m3 = mean(xb_derivation_m3) if e(sample)
gen hosplogodd_derivation_m3 = meanxb_derivation_m3 + re_derivation_m3
gen invlog_derivation_m3 = invlogit(hosplogodd_derivation_m3)

brier mort720_postdc pr_derivation_m3 

// standardized mortality ratio
total mort720_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m3 if derivation==1 
	matrix list e(b)
	local prderivationdeathm3=e(b)[1,1]
	di `prderivationdeathm3'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm3'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm3'	
di ((`derivationdeath' - `prderivationdeathm3')/`derivation')*10000


* INTERNAL VALIDATION *
estimates restore derivation_m3
predict pr_validation1_m3 if validation1==1	
predict re_validation1_m3 if validation1==1, reffects reses(rese_validation1_m3)
predict xb_validation1_m3 if validation1==1, xb
	replace xb_validation1_m3=. if validation1!=1
predict stdp_validation1_m3 if validation1==1, stdp 
	replace stdp_validation1_m3=. if validation1!=1

egen meanxb_validation1_m3 = mean(xb_validation1_m3) if validation1==1
gen hosplogodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3
gen invlog_validation1_m3 = invlogit(hosplogodd_validation1_m3)

gen hi_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 + 1.4*stdp_validation1_m3
gen lo_logodd_validation1_m3 = meanxb_validation1_m3 + re_validation1_m3 - 1.4*stdp_validation1_m3

gen hi_invlog_validation1_m3 = invlogit(hi_logodd_validation1_m3)
gen lo_invlog_validation1_m3 = invlogit(lo_logodd_validation1_m3)

brier mort720_postdc pr_validation1_m3


// standardized mortality ratio
total mort720_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m3 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm3=e(b)[1,1]
	di `prvalidation1deathm3'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm3'	

// excess deaths 
di `validation1death' - `prvalidation1deathm3'	
di ((`validation1death' - `prvalidation1deathm3')/`validation1')*10000


* VALIDATION 2021 * 
estimates restore derivation_m3
predict pr_validation2_m3 if validation2==1	
predict re_validation2_m3 if validation2==1, reffects reses(rese_validation2_m3)
predict xb_validation2_m3 if validation2==1, xb
	replace xb_validation2_m3=. if validation2!=1
predict stdp_validation2_m3 if validation2==1, stdp 
	replace stdp_validation2_m3=. if validation2!=1

egen meanxb_validation2_m3 = mean(xb_validation2_m3) if validation2==1
gen hosplogodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3
gen invlog_validation2_m3 = invlogit(hosplogodd_validation2_m3)

gen hi_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 + 1.4*stdp_validation2_m3
gen lo_logodd_validation2_m3 = meanxb_validation2_m3 + re_validation2_m3 - 1.4*stdp_validation2_m3

gen hi_invlog_validation2_m3 = invlogit(hi_logodd_validation2_m3)
gen lo_invlog_validation2_m3 = invlogit(lo_logodd_validation2_m3)

brier mort720_postdc pr_validation2_m3

// standardized mortality ratio
total mort720_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m3 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm3=e(b)[1,1]
	di `prvalidation2deathm3'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm3'	

// excess deaths 
di `validation2death' - `prvalidation2deathm3'	
di ((`validation2death' - `prvalidation2deathm3')/`validation2')*10000

	
*----------------------------------------
* M4: Base + Individual & Community 
*----------------------------------------

* DERIVATION *

local comorbid 																///
			cancer_nonmet cancer_met pulm chf dm_uncomp dm_comp	liver neuro ///
			renal htn cardic_arrhym valvular_d2 pulm_cir pvd paralysis 		///
			pud hypothyroid  lymphoma ra coag obesity wtloss fen			///
			anemia_cbl anemia_def   psychoses depression 

local aod 																	///
			 aod_lactate aod_kidney pressor_fromed_72hr aod_liver			///
			 aod_heme acute_respiratory_hosp 			
			
local sirs																	///
			sirs_temp sirs_rr sirs_pulse sirs_wbc

local individual															///
			married_impute housing_instab_combo_prior540 frailty_score				///
			can_impute
			
local community 															///	
			rural adi_impute
			
melogit mort720_postdc male age_sp* i.race_3cat `sirs' `aod' `comorbid' `individual' `community' ///
		if derivation==1 || sta6a: 
estimates store derivation_m4
predict pr_derivation_m4 if e(sample)	
predict re_derivation_m4 if e(sample), reffects reses(rese_derivation_m4)
predict xb_derivation_m4 if e(sample), xb
	replace xb_derivation_m4=. if derivation!=1

egen meanxb_derivation_m4 = mean(xb_derivation_m4) if e(sample)
gen hosplogodd_derivation_m4 = meanxb_derivation_m4 + re_derivation_m4
gen invlog_derivation_m4 = invlogit(hosplogodd_derivation_m4)

brier mort720_postdc pr_derivation_m4 

// standardized mortality ratio
total mort720_postdc if derivation==1 
	matrix list e(b)
	local derivationdeath=e(b)[1,1]
	di `derivationdeath'
total pr_derivation_m4 if derivation==1 
	matrix list e(b)
	local prderivationdeathm4=e(b)[1,1]
	di `prderivationdeathm4'	
total derivation if derivation==1
	matrix list e(b)
	local derivation=e(b)[1,1]
	di `derivation'

di `derivationdeath'/`prderivationdeathm4'	

// excess deaths 
di `derivationdeath' - `prderivationdeathm4'	
di ((`derivationdeath' - `prderivationdeathm4')/`derivation')*10000

preserve 
	keep if derivation==1 
	local ngroups 10
	xtile quantile_derivation_m4 = pr_derivation_m4, nq(`ngroups')
	bysort quantile_derivation_m4: sum pr_derivation_m4, de
	collapse (sum) mort720_postdc pr_derivation_m4 (count) numhosps=unique_hosp_count_id, by(quantile_derivation_m4)
	list numhosps mort720_postdc pr_derivation_m4, abb(20)
restore


* INTERNAL VALIDATION *
estimates restore derivation_m4
predict pr_validation1_m4 if validation1==1	
predict re_validation1_m4 if validation1==1, reffects reses(rese_validation1_m4)
predict xb_validation1_m4 if validation1==1, xb
	replace xb_validation1_m4=. if validation1!=1
predict stdp_validation1_m4 if validation1==1, stdp 
	replace stdp_validation1_m4=. if validation1!=1

egen meanxb_validation1_m4 = mean(xb_validation1_m4) if validation1==1
gen hosplogodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4
gen invlog_validation1_m4 = invlogit(hosplogodd_validation1_m4)

gen hi_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 + 1.4*stdp_validation1_m4
gen lo_logodd_validation1_m4 = meanxb_validation1_m4 + re_validation1_m4 - 1.4*stdp_validation1_m4

gen hi_invlog_validation1_m4 = invlogit(hi_logodd_validation1_m4)
gen lo_invlog_validation1_m4 = invlogit(lo_logodd_validation1_m4)

brier mort720_postdc pr_validation1_m4


// standardized mortality ratio
total mort720_postdc if validation1==1 
	matrix list e(b)
	local validation1death=e(b)[1,1]
	di `validation1death'
total pr_validation1_m4 if validation1==1 
	matrix list e(b)
	local prvalidation1deathm4=e(b)[1,1]
	di `prvalidation1deathm4'	
total validation1 if validation1==1
	matrix list e(b)
	local validation1=e(b)[1,1]
	di `validation1'

di `validation1death'/`prvalidation1deathm4'	

// excess deaths 
di `validation1death' - `prvalidation1deathm4'	
di ((`validation1death' - `prvalidation1deathm4')/`validation1')*10000

preserve 
	keep if validation1==1 
	local ngroups 10
	xtile quantile_validation1_m4 = pr_validation1_m4, nq(`ngroups')
	bysort quantile_validation1_m4: sum pr_validation1_m4, de
	collapse (sum) mort720_postdc pr_validation1_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation1_m4)
	list numhosps mort720_postdc pr_validation1_m4, abb(20)
restore


* VALIDATION 2021 * 
estimates restore derivation_m4
predict pr_validation2_m4 if validation2==1	
predict re_validation2_m4 if validation2==1, reffects reses(rese_validation2_m4)
predict xb_validation2_m4 if validation2==1, xb
	replace xb_validation2_m4=. if validation2!=1
predict stdp_validation2_m4 if validation2==1, stdp 
	replace stdp_validation2_m4=. if validation2!=1

egen meanxb_validation2_m4 = mean(xb_validation2_m4) if validation2==1
gen hosplogodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4
gen invlog_validation2_m4 = invlogit(hosplogodd_validation2_m4)

gen hi_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 + 1.4*stdp_validation2_m4
gen lo_logodd_validation2_m4 = meanxb_validation2_m4 + re_validation2_m4 - 1.4*stdp_validation2_m4

gen hi_invlog_validation2_m4 = invlogit(hi_logodd_validation2_m4)
gen lo_invlog_validation2_m4 = invlogit(lo_logodd_validation2_m4)

brier mort720_postdc pr_validation2_m4

// standardized mortality ratio
total mort720_postdc if validation2==1 
	matrix list e(b)
	local validation2death=e(b)[1,1]
	di `validation2death'
total pr_validation2_m4 if validation2==1 
	matrix list e(b)
	local prvalidation2deathm4=e(b)[1,1]
	di `prvalidation2deathm4'	
total validation2 if validation2==1
	matrix list e(b)
	local validation2=e(b)[1,1]
	di `validation2'

di `validation2death'/`prvalidation2deathm4'	

// excess deaths 
di `validation2death' - `prvalidation2deathm4'	
di ((`validation2death' - `prvalidation2deathm4')/`validation2')*10000


preserve 
	keep if validation2==1 
	local ngroups 10
	xtile quantile_validation2_m4 = pr_validation2_m4, nq(`ngroups')
	bysort quantile_validation2_m4: sum pr_validation2_m4, de
	collapse (sum) mort720_postdc pr_validation2_m4 (count) numhosps=unique_hosp_count_id, by(quantile_validation2_m4)
	list numhosps mort720_postdc pr_validation2_m4, abb(20)
restore


*-----------------------------------------
* Comparing models 0 & 1 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort720_postdc xb_derivation_m0 xb_derivation_m1 	///
					if derivation==1

* Internal validation data set 
roccomp mort720_postdc xb_validation1_m0 xb_validation1_m1 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort720_postdc xb_validation2_m0 xb_validation2_m1 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 2 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort720_postdc xb_derivation_m1 xb_derivation_m2 	///
					if derivation==1

* Internal validation data set 
roccomp mort720_postdc xb_validation1_m1 xb_validation1_m2 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort720_postdc xb_validation2_m1 xb_validation2_m2 	///
					if validation2==1

*-----------------------------------------
* Comparing models 1 & 3 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort720_postdc xb_derivation_m1 xb_derivation_m3 	///
					if derivation==1

* Internal validation data set 
roccomp mort720_postdc xb_validation1_m1 xb_validation1_m3 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort720_postdc xb_validation2_m1 xb_validation2_m3 	///
					if validation2==1
					
*-----------------------------------------
* Comparing models 1 & 4 using roccomp
*-----------------------------------------

* Derivation data set 
roccomp mort720_postdc xb_derivation_m1 xb_derivation_m4 	///
					if derivation==1

* Internal validation data set 
roccomp mort720_postdc xb_validation1_m1 xb_validation1_m4 	///
					if validation1==1,						

* Validation 2021 data set 
roccomp mort720_postdc xb_validation2_m1 xb_validation2_m4 	///
					if validation2==1
	

log close