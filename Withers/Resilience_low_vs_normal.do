*****Symbols*****
* = analysis step
// = comment
**----------- = new section ------------------

* Importing dataset
import spss using "https://github.com/quemeb/USC_research/raw/main/Withers/EXPLORING%20MENTAL%20HEALTH%20AND%20RESILIENCY.sav"

**#------------------------------- DATA CLEANING ---------------------------------
* Dropping made up dicho variables
ds *_01
drop `r(varlist)'
drop if Ageyears < 18

* Making all variables lower case... 
rename _all, lower

* Labeling insurance
label define insuf 1 "None" 2 "Partial" 4 "Full Private" 10 "Full Public"
label values insu insuf 

// dropping more stuff
drop if race01 == 4
drop if religion01 == 3 | religion01 == 4 | religion01 == 5 | religion01 == 6


*------------------ Creating MENTALHEALTH score from scratch 
foreach var in a035 a036 a037 a038 a039 a040 a041 a042 a043 a044 a045 a046 a047 {
	replace `var'=0 if `var'==2
}
gen mentalhealth = (a035 + a036 + a037 + a038 + a039 + a040 + a041 + a042 + a043 + a044 + a045 + a046 + a047)
drop if missing(mentalhealth) //drop any missing outcomes 
codebook mentalhealth

* Creating binary Mental Health outcome 
gen mental_cat = 0
replace mental_cat = 1 if mentalhealth >= 1
* Labeling data
label variable mental_cat "Presense of MENTAL HEALTH Illness" 
label define mental_catf 0 "No Mental Illness" 1 "Presense of mental illness"
label values mental_cat mental_catf 
codebook mental_cat
drop a035 a036 a037 a038 a039 a040 a041 a042 a043 a044 a045 a046 a047 mental mentalhealth 

*------------------ back to data cleaning
drop heightincm weightinkg ver surveyversion n respondentid // variables provide no info
drop a187 a188 a189 a190 a191 a174 a175 a176 // investigators not interested in these vars
drop a180 // information already captured by a179
drop _v2 // regrouped 

*------------------ RESILIENCE  

gen res = (a156 + a157 + a158 + a159 + a160 + a161)/6
codebook res 
drop if missing(res)
*ttest res == resilience, unpaired // comparing to V's score
*hist res, frequency normal title("Distribution of Resilience Score")
// we will use our score because they are different 
drop resilience a156 a157 a158 a159 a160 a161 
* Generating Ordinal Resilience 
gen res_cat = 0
replace res_cat = 1 if res >= 3 
replace res_cat = 2 if res >= 4.31
*now we want to create labels
label variable res_cat "Resilience Categories" 
label define res_catf 0 "Low Resilience(1.00-2.99)" 1 "Normal Resilience(3.00-4.30)" 2 "High Resilience(4.31-5.00)" 
label values res_cat res_catf 
*tab res_cat r3, chi2 exact // comparing to V's - was different 
drop r3 res

*------------------ COVID_EFFECT 

gen cov_effect = (a162 + a163 + a164 + a165 + a166 + a167)
codebook cov_effect
*hist cov_effect, frequency normal title("Distribution of cov_effect")
*ttest cov_effect==covid_effect, unpaired // comparing to V's
drop covid_effect a162 a163 a164 a165 a166 a167

*------------------ COVID_FEAR 

gen cov_fear = (a168 + a169 + a170 + a171 + a172)
codebook cov_fear
//hist cov_fear, frequency normal title("distribution of cov_fear")
drop covid_fear a168 a169 a170 a171 a172

*------------------ COVID_PHYCH 

gen cov_psych = a173 + a177 + a178 + a179 + a181 + a182 + a183 + a184 + a185 + a186
codebook cov_psych
//hist cov_psych, frequency normal title("distribution of cov_psych")
drop covid_psych a173 a177 a178 a179 a181 a182 a183 a184 a185 a186

*------------------ WELLBEING

*generate wellbeing
gen cov_wellbeing = a211 + a212 + a213 + a214 + a215
codebook cov_wellbeing
codebook wellbeing
//hist cov_wellbeing, frequency normal title ("dist of cov_wellbeing")
//hist wellbeing, frequency normal title("Distribution of wellbeing")
signrank cov_wellbeing=wellbeing 
*make categories
gen well_cat = 0 if cov_wellbeing < 17.5
replace well_cat = 1 if cov_wellbeing >= 17.5 
*now we want to create labels
label variable well_cat "wellbeing Categories" 
label define well_catf 0 "good" 1 "poor" 
label values well_cat well_catf 
drop a211 a212 a213 a214 a215 wellbeing wellbeing01 cov_wellbeing


**# ------------------------------ NEW SUMMARY SCORES ----------------------
* -- Generating COVID_BURDEN 

*scaling variables
global burden "a192 a193 a194 a195 a196 a197 a198 a199 a200 a201 a202 a203 a204 a205 a206 a207 a208 a209 a210"
foreach var in $burden {
	replace `var' = 0 if `var' == 77
}
gen cov_burden = (a192 + a193 + a194 + a195 + a196 + a197 + a198 + a199 + a200 + a201 + a202 + a203 + a204 + a205 + a206 + a207 + a208 + a209 + a210)
codebook cov_burden
drop a192 a193 a194 a195 a196 a197 a198 a199 a200 a201 a202 a203 a204 a205 a206 a207 a208 a209 a210


* -- Generating COVID_CONTACT

*scaling variables 
foreach var in a002 a003 a005 {
    gen temp1 = (`var' == 1)
    gen temp2 = (`var' == 2)
    gen temp_missing = missing(`var')

    replace `var' = 1 if temp2
    replace `var' = 2 if temp1
    replace `var' = 0 if temp_missing

    drop temp1 temp2 temp_missing
}
gen cov_contact = a002 + a003 + a005
codebook cov_contact
drop a002 a003 a005


**# ------- Creating lists

vl set
vl list vlcontinuous
vl list vluncertain 
vl move (ageyears cov_effect cov_fear cov_psych cov_burden cov_contact) vlcontinuous
vl list vlcategorical
vl drop (mental_cat) // dropping since it is our outcome variable for macros
vl drop (res_cat)   // another outcome 

**** making it regular logistic 
drop if res_cat == 2
recode res_cat (1=0) (0=1)
label define res_catf1 0 "Normal Resilience" 1 "Low Resilience" 
label values res_cat res_catf1 

* ---------- LINEARITY CHECK

* height_cm - linear 
fp <height_cm>, scale center replace: logit res_cat <height_cm> 

* weight_kg - linear 
//lowess mental_cat weight_kg, logit
fp <weight_kg>, scale center replace: logit res_cat <weight_kg>

* ageyears - cosine 
//lowess res_cat ageyears, logit 
//fp <ageyears>, scale center replace: logit res_cat age18 <ageyears>
gen age_cosine = cos(ageyears)
//lowess res_cat age_cosine, logit 
//fp <age_cosine>, scale center replace: logit res_cat <age_cosine>
*drop ageyears 

* cov_effect - cosine  
//lowess res_cat cov_effect, logit 
//fp <cov_effect>, scale center replace: logit res_cat <cov_effect>
gen cov_effect_cos = cos(cov_effect)
//lowess res_cat cov_effect_cos, logit 
//fp <cov_effect_cos>, scale center replace: logit res_cat <cov_effect_cos>
drop cov_effect 

* cov_fear - linear
//lowess res_cat cov_fear, logit 
//fp <cov_fear>, scale center replace: logit res_cat <cov_fear>

* cov_psych - linear 
//lowess res_cat cov_psych, logit 
//fp <cov_psych>, scale center replace: logit res_cat <cov_psych>

* cov_burden - linear 
//lowess res_cat cov_burden, logit 
//fp <cov_burden>, scale center replace: logit res_cat <cov_burden>

* cov_contact - log 
//lowess res_cat cov_contact, logit 
//fp <cov_contact>, scale center replace: logit res_cat <cov_contact>
gen log_cov_contact = log(cov_contact)
//lowess res_cat log_cov_contact, logit 
//fp <log_cov_contact>, scale center replace: logit res_cat <log_cov_contact>
drop cov_contact 

vl rebuild 

**# Bivariant analysis for Mental Health 

* Categorical 
global cat_variables "$vlcategorical mental_cat"
global cat_prelim ""

foreach i in $cat_variables{
	// create a two-way table of the predictor and the response variable
	tab `i' res_cat, exact  
	// if the p-value is less than 0.25, add the variable to the list
	if r(p_exact) < 0.25 {
		global cat_prelim  "$cat_prelim `i' "
	}
}
// display the variables with p-value < 0.25
di "$cat_prelim"   // 

* Removed categoricals
global leftover_cat ""
foreach var in $cat_variables {
    local found 0
    foreach pre in $cat_prelim {
        if ("`var'" == "`pre'") {
            local found 1
        }
    }
    if (`found' == 0) {
        global leftover_cat "$leftover_cat `var'"
    }
}

// Display the unique variables
display "The removed categorical variables are: $leftover_cat"



* ---- Continouos ---
global cont_variables "$vlcontinuous age_cosine cov_effect_cos log_cov_contact"
global cont_prelim ""

foreach i in $cont_variables{
	// run the logistic regression
	logit res_cat `i', nolog 
	// get the p-value
    testparm `i'
	// if the p-value is less than 0.25, add the variable to the list
    if r(p) < 0.25 {
        // Note the addition of $ before prelim_vars
        global cont_prelim "$cont_prelim `i'"
    }
}
// display the variables with p-value < 0.25
di "$cont_prelim"

* Removed continuous 
global leftover_cont ""
foreach var in $cont_variables {
    local found 0
    foreach pre in $cont_prelim {
        if ("`var'" == "`pre'") {
            local found 1
        }
    }
    if (`found' == 0) {
        global leftover_cont "$leftover_cont `var'"
    }
}

// Display the unique variables
display "The removed categorical variables are: $leftover_cont"


* ---- Fixing categoricals 

global cat_prelim_expanded ""
foreach var in $cat_prelim {
	global cat_prelim_expanded "$cat_prelim_expanded (i.`var')"
}
di "$cat_prelim_expanded"

* prelim - stepwise models 
stepwise, pr(0.1): logit res_cat $cat_prelim_expanded $cont_prelim , nolog or

* checking for significance in the terms taken out earlier
global leftover_cat_expanded ""
foreach var in $leftover_cat {
	global leftover_cat_expanded "$leftover_cat_expanded (i.`var')"
}
di "$leftover_cat_expanded"

* ---------- Final stepwise 
sw, pr(0.1): logit res_cat cov_burden cov_psych log_cov_contact sex01 (i.residence01) (i.livingstatus01) $leftover_cat_expanded $leftover_cont, nolog or 
// 

* ---------- Preliminary model after stepwises 
logit res_cat log_cov_contact cov_burden cov_psych i._v3 ib2.sex01 i.residence01 ib3.livingstatus01, nolog or 


* --- combining categories 
test 1.livingstatus01 = 2.livingstatus01
gen living_combined = 0
replace living_combined = 1 if livingstatus01 == 1 | livingstatus01 == 2
label variable living_combined "Combined livingstatus" 
label define living_combinedf 0 "I live with other people than family" 1 "I live alone or with family"
label values living_combined living_combinedf 



* --- prelim final model

logit res_cat log_cov_contact cov_burden cov_psych _v3 sex01 i.residence01 i.living_combined, nolog or 

*# _______________________ INTERACTIONS ___________________________ 
// no interactions were found 
global potential_interactions "c.ageyear c.age_cosine race01"

*foreach i in $potential_interactions {
*	logit res_cat c.log_cov_contact##`i' c.cov_burden##`i' c.cov_psych##`i' c._v3##`i' c.sex01##`i' residence01##`i' living_combined##`i', nolog or 
*}

*# ______________ CONFOUNDERS ____________________

// Define a macro for potential confounders
global potential_confounders "age_cosine i.race01"
global confounders_cont ""

// Loop through each main independent variable and each potential confounder
foreach main in log_cov_contact cov_burden cov_psych _v3 sex01 living_combined {
    foreach conf in $potential_confounders {
        di "Testing for confounding effect of `conf' on `main'..."

        // Run the logistic regression model without the confounder
        logit res_cat `main'
        // Store the coefficient for the main variable
        scalar b1 = _b[`main']

        // Run the logistic regression model with the confounder
        logit res_cat `main' `conf'
        // Store the coefficient for the main variable
        scalar b2 = _b[`main']

        // Calculate the percentage change in the coefficient
       scalar perc_change = 100 * (b2 - b1) / b1
         di "Percent change in coefficient for `main' when including `conf': " float(perc_change) "%"
		if (abs(perc_change) > 10){
			global confounders_cont "$confounders_cont `conf'"
		}
    }
}
di "$confounders_cont"


foreach conf in $potential_confounders {
	logit res_cat i.residence01, nolog or
	scalar b1 = _b[2.residence01]
	scalar b2 = _b[3.residence01]
	
	logit res_cat i.residence01 `conf', nolog or
	scalar b11 = _b[2.residence01]
	scalar b22 = _b[3.residence01]


	scalar perc_change1 = 100 * (b1 - b11) / b1 
	scalar perc_change2 = 100 * (b2 - b22) / b2
	di "Percent change in coefficient for Rural: " float(perc_change1) "%"
	di "Percent change in coefficient for Suburban: " float(perc_change2) "%"
}
// counfounders = no confounding 


**# * ------------------ FINAL MENTAL MODEL --------------------
logit res_cat log_cov_contact cov_burden cov_psych ib2._v3 i.sex01 i.residence01 i.living_combined i.race01, nolog or


* ----------------- Model diagnostics ----------------------

estat gof, group(10) table 
// doesn't depart from goodness

predict predicted, p
predict delta_x2, dx2
predict delta_dev, ddeviance
predict delta_beta, dbeta

// change in Pearson GOF
scatter delta_x2 p

// change in Deviance GOF
scatter delta_dev p 

// change in cooks distance 
scatter delta_beta p 
twoway scatter delta_x2 p [aweight = delta_beta], msymbol(circle_hollow)

* --------------- DROPPING PROBLEMATIC OBSERVATIONS - SENSITIVITY ANALYSIS 
*drop if delta_x2 > 4 & delta_x2 != .
*logit res_cat log_cov_contact cov_burden cov_psych ib2._v3 i.sex01 i.residence01 i.living_combined i.race01, nolog or
*estat gof, group(10) table

*drop p delta_x2 delta_dev delta_beta 

*predict predicted, p
*predict delta_x2, dx2
*predict delta_dev, ddeviance
*predict delta_beta, dbeta

// change in Pearson GOF
*scatter delta_x2 p

// change in Deviance GOF
*scatter delta_dev p 

// change in cooks distance 
*scatter delta_beta p 
*twoway scatter delta_x2 p [aweight = delta_beta], msymbol(circle_hollow)

* ---------------- Predictions ---------------------------

logit res_cat log_cov_contact cov_burden cov_psych ib2._v3 i.sex01 i.residence01 i.living_combined i.race01, nolog or

* FINAL FINAL MODEL AFTER SENSITIVITY 
logit res_cat log_cov_contact cov_burden cov_psych i.sex01 i.residence01 i.race01, nolog or

estat classification
lroc 
lsens
predict p
cutpt res_cat p

logit res_cat log_cov_contact cov_burden cov_psych i.sex01 i.residence01 i.race01, nolog or
estat clas, cut(.26033308)


