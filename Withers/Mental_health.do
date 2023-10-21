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
drop if Ageyears < 18 | Ageyears > 75

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

egen res = rowmean(a156 a157 a158 a159 a160 a161)
codebook res 
drop if missing(res)
hist res, frequency normal title("Distribution of Resilience Score")
// we will use our score because they are different 
drop resilience a156 a157 a158 a159 a160 a161 
* Generating Ordinal Resilience 

gen res_cat = .
replace res_cat = 0 if res < 3.167 - .485713
replace res_cat = 1 if res >= 3.167 - .485713 & res != .
replace res_cat = 2 if res >= 3.167 + .485713 & res != .

*gen res_cat = 0
*replace res_cat = 1 if res >= 3 
*replace res_cat = 2 if res >= 4.31
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

gen cov_psych = (a173 + a177 + a178 + a179 + a181 + a182 + a183 + a184 + a185 + a186)
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


* ---------- LINEARITY CHECK

* height_cm - linear 
//lowess mental_cat height_cm, logit 
//fp <height_cm>, scale center replace: logit mental_cat <height_cm>

* weight_kg - linear 
//lowess mental_cat weight_kg, logit
//fp <weight_kg>, scale center replace: logit mental_cat <weight_kg>

* ageyears - linear 
//lowess mental_cat ageyears, logit 
//fp <ageyears>, scale center replace: logit mental_cat <ageyears>

* cov_effect - linear 
//lowess mental_cat cov_effect, logit 
//fp <cov_effect>, scale center replace: logit mental_cat <cov_effect>

* cov_fear - log 
//lowess mental_cat cov_fear, logit 
//fp <cov_fear>, scale center replace: logit mental_cat <cov_fear>
gen log_cov_fear = log(cov_fear)
//lowess mental_cat log_cov_fear, logit 
//fp <log_cov_fear>, scale center replace: logit mental_cat <log_cov_fear>
drop cov_fear 

* cov_psych - log
//lowess mental_cat cov_psych, logit 
//fp <cov_psych>, scale center replace: logit mental_cat <cov_psych>
gen log_cov_psych = log(cov_psych)
//lowess mental_cat log_cov_psych, logit 
//fp <log_cov_psych>, scale center replace: logit mental_cat <log_cov_psych>
drop cov_psych 

* cov_burden - linear 
//lowess mental_cat cov_burden, logit 
//fp <cov_burden>, scale center replace: logit mental_cat <cov_burden>

* cov_contact - linear 
//lowess mental_cat cov_contact, logit 
//fp <cov_contact>, scale center replace: logit mental_cat <cov_contact>


vl rebuild 

**# Bivariant analysis for Mental Health 

* Categorical 
global cat_variables "$vlcategorical res_cat"
global cat_prelim ""

foreach i in $cat_variables{
	// create a two-way table of the predictor and the response variable
	tab `i' mental_cat, exact  
	// if the p-value is less than 0.25, add the variable to the list
	if r(p_exact) < 0.25 {
		global cat_prelim  "$cat_prelim `i' "
	}
}
// display the variables with p-value < 0.25
di "$cat_prelim"   // dropped Resilience, Sex, Work-from-home 

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
global cont_variables "$vlcontinuous log_cov_fear log_cov_psych"
global cont_prelim ""

foreach i in $cont_variables{
	// run the logistic regression
	logit mental_cat `i', nolog 
	// get the p-value
    testparm `i'
	// if the p-value is less than 0.25, add the variable to the list
    if r(p) < 0.25 {
        // Note the addition of $ before prelim_vars
        global cont_prelim "$cont_prelim `i'"
    }
}
// display the variables with p-value < 0.25
di "$cont_prelim"	// dropped nothing

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


**# prelim - stepwise models 
stepwise, pr(0.1): logit mental_cat $cat_prelim_expanded $cont_prelim , nolog or

// need to explain the log-transformed meaning increase 
// ex. OR=2, an increase in independent variable by a factor of ~2.7 will double the odds of the mental_cat 

* checking for significance in the terms taken out earlier
global leftover_cat_expanded ""
foreach var in $leftover_cat {
	global leftover_cat_expanded "$leftover_cat_expanded (i.`var')"
}
di "$leftover_cat_expanded"

* ---------- Final stepwise 
sw, pr(0.1): logit mental_cat region01 (i.insu) (ib7.religion01) (ib3._v3) (i.well_cat) log_cov_psych cov_contact (ib3._v1) $leftover_cat_expanded $leftover_cont, nolog or 
// no lefovers came back significant 

* ---------- Preliminary model after stepwises 
logit mental_cat region01 i.insu i.religion01 _v3 i.well_cat log_cov_psych cov_contact ib3._v1, nolog or 

* --- combining categories 
test 4.insu = 10.insu
gen insu_full = 0
replace insu_full = 1 if insu == 2
replace insu_full = 2 if insu == 4 | insu == 10
label variable insu_full "Combined insurances" 
label define insu_fullf 0 "no insurance" 1 "Partial" 2 "Full Public/Private"
label values insu_full insu_fullf 


* --- prelim final model
logit mental_cat region01 i._v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01, nolog or 

*# _______________________ INTERACTIONS ___________________________ 
// no interactions were found 
*global potential_interactions "c.sex01 c.ageyears race01"

*foreach i in $potential_interactions {
*	logit mental_cat c.region01##`i' insu_full##`i' ib7.religion01##`i' c._v3##`i' c.well_cat##`i' c.log_cov_psych##`i' c.cov_contact##`i' ib3._v1##`i', nolog or 
*}

/// lrtest for sex01##religion01 and race01##_v1 
*logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 sex01, nolog or 
*estimates store A
*logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01##sex01, nolog or
*lrtest A 
// race01##_v1 
*logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 i.race01, nolog or 
*estimates store A
*logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1##race i.insu_full ib7.religion01, nolog or
*lrtest A


* --------- Model after interactions 
logit mental_cat region01 i.insu_full ib7.religion01 i._v3 well_cat log_cov_psych cov_contact ib3._v1, nolog or 

*# ______________ CONFOUNDERS ____________________

// Define a macro for potential confounders
global potential_confounders "ageyears sex01 i.race "
global confounders_cont ""

// Loop through each main independent variable and each potential confounder
*foreach main in region01 _v3 well_cat log_cov_psych cov_contact {
*    foreach conf in $potential_confounders {
*        di "Testing for confounding effect of `conf' on `main'..."

        // Run the logistic regression model without the confounder
*        logit mental_cat `main'
        // Store the coefficient for the main variable
*        scalar b1 = _b[`main']

        // Run the logistic regression model with the confounder
*        logit mental_cat `main' `conf'
        // Store the coefficient for the main variable
*        scalar b2 = _b[`main']

        // Calculate the percentage change in the coefficient
*       scalar perc_change = 100 * (b2 - b1) / b1
*         di "Percent change in coefficient for `main' when including `conf': " float(perc_change) "%"
*		if (abs(perc_change) > 10){
*			global confounders_cont "$confounders_cont `conf'"
*		}
*    }
*}
*di "$confounders_cont"
// confounders = race/_v3 (13%), race/well_cat (11.4%), race/log_cov_psych (9.8%), race/cov_contact (58%)


// Loop through each main independent variable and each potential confounder
	
*foreach conf in $potential_confounders {
*	logit mental_cat i.insu_full
*	scalar b1 = _b[1.insu_full]
*	scalar b2 = _b[2.insu_full]
	
*	logit mental_cat i.insu_full `conf'
*	scalar b11 = _b[1.insu_full]
*	scalar b22 = _b[2.insu_full]
	
*	scalar perc_change1 = 100 * (b1 - b11) / b1 
*	scalar perc_change2 = 100 * (b2 - b22) / b2
*	di "Percent change in coefficient for Partial: " float(perc_change1) "%"
*	di "Percent change in coefficient for Full: " float(perc_change2) "%"
*}
// counfounders = ageyears (17%, 18.5%), race (170%, 1.1%)

*foreach conf in $potential_confounders {
*	logit mental_cat ib7.religion01, nolog or
*	scalar b1 = _b[1.religion01]
*	scalar b2 = _b[2.religion01]
*	scalar b8 = _b[8.religion01]

*	logit mental_cat ib7.religion01 `conf', nolog or
*	scalar b11 = _b[1.religion01]
*	scalar b22 = _b[2.religion01]
*	scalar b88 = _b[8.religion01]

*	scalar perc_change1 = 100 * (b1 - b11) / b1 
*	scalar perc_change2 = 100 * (b2 - b22) / b2
*	scalar perc_change8 = 100 * (b8 - b88) / b8
*	di "Percent change in coefficient for Christian/Catholic: " float(perc_change1) "%"
*	di "Percent change in coefficient for Buddhist: " float(perc_change2) "%"
*	di "Percent change in coefficient for Other: " float(perc_change8) "%"
*}
// counfounders = ageyears (7%, 41%, 6%), race (113%, 141%, 62%)

*foreach conf in $potential_confounders {
*	logit mental_cat ib3._v1
*	scalar b1 = _b[1._v1]
*	scalar b2 = _b[2._v1]
*
*	logit mental_cat ib3._v1 `conf'
*	scalar b11 = _b[1._v1]
*	scalar b22 = _b[2._v1]

*	scalar perc_change1 = 100 * (b1 - b11) / b1 
*	scalar perc_change2 = 100 * (b2 - b22) / b2
*	di "Percent change in coefficient for Improved: " float(perc_change1) "%"
*	di "Percent change in coefficient for Worse: " float(perc_change2) "%"
*}
// counfounders = race (75%, 15%)


**# * ------------------ FINAL MENTAL MODEL --------------------
logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 i.race, nolog or 


* ----------------- Model diagnostics ----------------------

estat gof, group(10) table 
// doesn't depart from goodness

collin region01 _v3 well_cat log_cov_psych cov_contact
// no collinearity 

predict predicted, p
predict delta_x2, dx2
predict delta_dev, ddeviance
predict delta_beta, dbeta

// change in Pearson GOF
scatter delta_x2 p
list mental_cat region01 _v3 well_cat log_cov_psych cov_contact _v1 insu_full religion01 p delta_x2 if delta_x2 > 10 & delta_x2 != .

// change in Deviance GOF
scatter delta_dev p 
list mental_cat p delta_dev if delta_dev > 6 & delta_dev != .

// change in cooks distance 
scatter delta_beta p 
list mental_cat p delta_beta if delta_beta > 2 & delta_beta != .
twoway scatter delta_x2 p [aweight = delta_beta], msymbol(circle_hollow)

* --------------- DROPPING PROBLEMATIC OBSERVATIONS - SENSITIVITY ANALYSIS 
*drop if delta_x2 > 4 & delta_x2 != .
*logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 i.race, nolog or 
*estat gof, group(4) table

*drop p delta_x2 delta_dev delta_beta 

*predict predicted, p
*predict delta_x2, dx2
*predict delta_dev, ddeviance
*predict delta_beta, dbeta

// change in Pearson GOF
*scatter delta_x2 p
*list mental_cat region01 _v3 well_cat log_cov_psych cov_contact _v1 insu_full religion01 p delta_x2 if delta_x2 > 10 & delta_x2 != .

// change in Deviance GOF
*scatter delta_dev p 
*list mental_cat p delta_dev if delta_dev > 4 & delta_dev != .

// change in cooks distance 
*scatter delta_beta p 
*list mental_cat p delta_beta if delta_beta > 2 & delta_beta != .
*twoway scatter delta_x2 p [aweight = delta_beta], msymbol(circle_hollow)


// NOT RELEVANT SINCE WE DID A INFERENCE MODEL PREVIOUSLY, just for fun lol 
* ---------------- Predictions ---------------------------

logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 i.race, nolog or 

estat classification
lroc 
lsens
predict p
cutpt mental_cat p

logit mental_cat region01 _v3 well_cat log_cov_psych cov_contact ib3._v1 i.insu_full ib7.religion01 i.race, nolog or 
estat clas, cut(.13861314)



logit mental_cat region01 _v3 ib3._v1 well_cat cov_effect height_cm cov_contact sex01 log_cov_fear log_cov_psych ageyears cov_burden c.wfh c.religion01 

estat classification
lroc 
lsens
predict p
cutpt mental_cat p


**# --------------------------- RESILIENCE MODEL ------------------------------


*find cat and cont prelim

*find leftovers 

stepwise, pr(0.1): mlogit res_cat $cat_prelim_expanded $cont_prelim, nolog rrr













logit mental_cat height_cm, nolog or
estimates store A
logit mental_cat height_cm_1, nolog or
lrtest A

*drop old variables - innacurate ones
drop heightincm weightinkg


**#renaming variables...
*age
rename ageyears age
rename sex01 sex 
drop if age < 18 		//2 observations deleted
*sex
drop if sex == .		// 145 observations 
*country
drop r3

encode country1012, generate(country) // converting country1012
drop country1012




*********************** Mental Analysis *********************
*dropping missing mental
*142 observations dropped
**# Bookmark #1
drop if missing(mental) 

*recoding variables for analysis 
replace mental = 0 if mental == 2
label define mental_l 0 "No" 1 "Yes"
label values mental mental_l 


stepwise, pr(0.1) pe(0.09): logit mental $vlcontinuous
stepwise, pr(0.1) pe(0.09): logit mental $icategorical

stepwise, pr(0.05): logit mental $vlcontinuous



stepwise, pr(0.1) pe(0.09): logit mental $vlcontinuous i.$vlcategorical $vluncertain





*creating labels for clean height and weight
label variable height_cm "Height in cms"
label variable weight_kg "Weight in kgs"

*Work status collapsed
drop _v2
label variable Work_regp "Work or study status"

***** collinear data 



*create automatic lists
vl set

vl list vlcontinuous 
vl move (RespondentID N) vlother

vl list vlcategorical 
vl drop (mental R3)
vl substitute icategorical = i.vlcategorical

vl list vluncertain
vl move (Ageyears COVID_effect COVID_fear COVID_psych wellbeing) vlcontinuous
*vl rebuild 



************* Resilience Analysis ************** 

*drop missing Resilience
drop if missing(Resilience)

*generating resilience categorical
gen resilience_cat = 0
replace resilience_cat = 1 if Resilience >= 1.67
replace resilience_cat = 2 if Resilience >= 3.34

drop Resilience 




drop wellbeing01

*dropping missing missing variables 
drop if missing(Resilience)


*create automatic lists
vl set

vl list vlcategorical 
vl list vlcontinuous 
vl list vluncertain

*to move variables to continous
vl move (Ageyears COVID_effect COVID_fear COVID_psych wellbeing) vlcontinuous

*dropping mental out of our lists
vl drop (mental)

* need to drop this "perfect prediction factor"
vl drop (wellbeing01) 
vl drop (wellbeing)
vl rebuild 



lasso logit wellbeing01 Resilience $vlcategorical $vlcontinuous











******************* EXTRA STUFF *************

*splitting dataset 
set seed 123
splitsample, generate(sample) split(0.2 0.8)
label define svalues 2 "Training" 1 "Testing"
label values sample svalues 

**********************LASSO
lasso logit wellbeing01 Resilience $vlcategorical $vlcontinuous if sample == 2, rseed(123)
cvplot

********************** Elasticnet
elasticnet logit wellbeing01 Resilience $vlcategorical $vlcontinuous, rseed(123)


********************** Ridge Regression
elasticnet logit wellbeing01 Resilience $vlcategorical $vlcontinuous, rseed(123) alpha(0)

