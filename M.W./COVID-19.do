*****Symbols*****
* = analysis step
// = comment
**----------- = new section ------------------

* Importing dataset
import spss using "C:\Users\bryan\Desktop\EXPLORING MENTAL HEALTH AND RESILIENCY.sav"

**#------------------------------- DATA CLEANING ---------------------------------
* Dropping made up dicho variables
ds *_01
drop `r(varlist)'
drop if Ageyears < 18

* Making all variables lower case... 
rename _all, lower

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

*------------------ Generating RESILIENCE score from scratch 
gen res = (a156 + a157 + a158 + a159 + a160 + a161)/6
codebook res 
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
drop r3 

*------------------ COVID_EFFECT 

gen cov_effect = (a162 + a163 + a164 + a165 + a166 + a167)
codebook cov_effect
*hist cov_effect, frequency normal title("Distribution of cov_effect")
*ttest cov_effect==covid_effect, unpaired // comparing to V's
drop covid_effect a162 a163 a164 a165 a166 a167

*------------------ COVID_FEAR 

gen cov_fear = (a168 + a169 + a170 + a171 + a172)
codebook cov_fear
hist cov_fear, frequency normal title("distribution of cov_fear")
drop covid_fear a168 a169 a170 a171 a172

*------------------ COVID_PHYCH 

gen cov_psych = a173 + a177 + a178 + a179 + a181 + a182 + a183 + a184 + a185 + a186
codebook cov_psych
hist cov_psych, frequency normal title("distribution of cov_psych")
drop covid_psych a173 a177 a178 a179 a181 a182 a183 a184 a185 a186

*------------------ WELLBEING

*generate wellbeing
gen cov_wellbeing = a211 + a212 + a213 + a214 + a215
codebook cov_wellbeing
codebook wellbeing
hist cov_wellbeing, frequency normal title ("dist of cov_wellbeing")
hist wellbeing, frequency normal title("Distribution of wellbeing")
signrank cov_wellbeing=wellbeing 
*make categories
drop if missing(res) 
gen well_cat=0 if cov_wellbeing<17.5
replace well_cat=1 if cov_wellbeing>=17.5 
*now we want to create labels
label variable well_cat "wellbeing Categories" 
label define well_catf 0 "poor" 1 "good" 
label values well_cat well_catf 
drop a211 a212 a213 a214 a215 wellbeing wellbeing01



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
gen cov_contact = (a002 + a003 + a005)
codebook cov_contact
drop a002 a003 a005


*Generating categorical resilience
drop if missing(res)
gen res_cat = 0
replace res_cat = 1 if res >= 3
replace res_cat = 2 if res >= 4.31
*labeling categories
label variable res_cat "Resilience Categories"
label define res_catf 0 "Low Resilience (1.00-2.99)" 1 "Normal Resilience (3.00-4.30)" 2 "High Resilience (4.31-5.00)"
label values res_cat res_catf
codebook res_cat
codebook r3 // Victor's
tab res_cat r3, chi2 exact

*significant Fishers so drop r3
drop r3 




vl set


* ---------- LINEARITY CHECK

* height_cm - need to include 2 term polynomial 
lowess mental_cat height_cm, logit 
fp <height_cm>, scale center replace: logit mental_cat <height_cm>
gen log_height_cm = log(height_cm)
lowess mental_cat log_height_cm, logit 
fp <log_height_cm>, scale center replace: logit mental_cat <log_height_cm>

* weight_kg - log 
lowess mental_cat weight_kg, logit
gen log_weight_kg = log(weight_kg)
lowess mental_cat log_weight_kg, logit
fp <log_weight_kg>, scale center replace: logit mental_cat <log_weight_kg>

* 
fp <





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
