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

encode country1012, gen(country)

vl set
vl list vlcontinuous
vl list vluncertain 
vl move (ageyears cov_effect cov_fear cov_psych cov_burden cov_contact) vlcontinuous
vl move (country) vlcategorical 
vl list vlcategorical
vl drop (mental_cat) // dropping since it is our outcome variable for macros
vl drop (res_cat)   // another outcome 

global y "mental_cat"
global x "$vlcontinuous $vlcategorical res_cat"

* Split the dataset into training and testing
gen split = runiform() < 0.7

* Logistic regression
logit mental_cat $x if split
predict prob_logit if !split
gen pred_logit = (prob_logit > 0.5)

*knn 
knn train($x) trainoutcome(mental_cat) test($x) if !split, k(5) gen(pred_knn)
