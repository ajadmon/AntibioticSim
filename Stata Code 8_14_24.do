* Setting seed, specifying Stata version, and dropping programs in memory for reproducibility
	clear
	program drop _all
	set seed 1809917878
	version 18

* Importing final, deidentified cohort from Chanderraj, et al, JAMA IM, 2024
	import delimited "/Users/ajadmon/Dropbox (University of Michigan)/Shared with Me/IV_Cohort_Data/deindentified_iv_cohort.csv" , clear
	cd "/Users/ajadmon/Dropbox (University of Michigan)/Research/Bob and Rishi/Stata Code and Data"

* Encoding and labeling variables from imported dataset
	encode infection , generate (infection_encoded)
	encode sex , generate (sex_encoded)

	generate race_cat_encoded = .
	replace race_cat_encoded = 1 if race_cat == "African American"
	replace race_cat_encoded = 2 if race_cat == "Caucasian"
	replace race_cat_encoded = 3 if race_cat == "American Indian or Alaska Native"
	replace race_cat_encoded = 4 if race_cat == "Native Hawaiian and Other Pacific Islander" | race_cat == "Asian"
	replace race_cat_encoded = 5 if race_cat == "Other" | race_cat == "Other" | race_cat == "Patient Refused" | race_cat == "NA" | race_cat == "Unknown"

	generate race_cat_collapsed = .
	replace race_cat_collapsed = 1 if race_cat_encoded == 1
	replace race_cat_collapsed = 2 if race_cat_encoded == 2
	replace race_cat_collapsed = 3 if race_cat_encoded == 4
	replace race_cat_collapsed = 4 if race_cat_encoded == 3 | race_cat_encoded == 5
	label var race_cat_collapsed "Race Labels - Collapsed"
	label define race_cat_collapsed_lab 1 "Black" 2 "White" 3 "Asian or Pacific Islander" 4 "Other or Unknown"
	label values race_cat_collapsed race_cat_collapsed_lab

	generate zosyn_encoded = .
	replace zosyn_encoded = 0 if zosyn == "cefepime"
	replace zosyn_encoded = 1 if zosyn == "zosyn"

	destring charlson_score , force replace /*6 set as NA. Keeping them missing*/


	tabulate sex_encoded , generate(sex_encoded_)
	tabulate infection_encoded , generate(infection_encoded_)
	tabulate race_cat_collapsed , generate(race_cat_collapsed_)

* Fitting model for 90-day mortality. The covariate coefficients (for variables besides antibiotic choice, which we will "set") will be used in the simulation later.	
	logit ninety_day_outcome i.zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4


* Generating sampling weights to scale the population up to a national population. These assume that a) 1.7 million people are admitted nationally with sepsis (Rhee, et al, JAMA, 2017), b) that 50% of admitted patients are equally eligible for antianaerobic and anaerobe-preserving antibiotics (Novosad, et al, MMWR, 2017), and c) that patients in this study are representative of hospitalizes, septic patients nationally with regards to other factors (e.g., severity of illness, comorbidities, etc.). (1,700,000 * 0.5) / 7,569 = 112.3 
	generate sampleweights = 112.3
	generate identityweights = 1


* Generating alternative sample weights if anti-anaerobic antibiotic use is 20% higher nationally than in the study population, or if anti-anaerobic antibiotic use is 20% lower nationally than in the study population. This accounts for differences in antibiotic prescribing practice
	gen sampleweights_increased = sampleweights * 1.2 if zosyn_encoded == 1
	sum sampleweights
	scalar total_weights = r(sum)
	sum sampleweights_increased if zosyn_encoded == 1
	scalar total_weights_inc_zosyn = r(sum)
	scalar remain_weight_inc_zosyn = total_weights - total_weights_inc_zosyn
	summarize identityweights if zosyn_encoded == 0
	scalar untreated = r(sum)
	replace sampleweights_increased = remain_weight_inc_zosyn / untreated if zosyn_encoded == 0

	gen sampleweights_decreased = sampleweights * 0.8 if zosyn_encoded == 1
	sum sampleweights_decreased if zosyn_encoded == 1
	scalar total_weights_dec_zosyn = r(sum)
	scalar remain_weight_dec_zosyn = total_weights - total_weights_dec_zosyn
	replace sampleweights_decreased = remain_weight_dec_zosyn / untreated if zosyn_encoded == 0

*Saving datasets for simulations
	save data_for_sim.dta , replace
	keep ninety_day_outcome zosyn_encoded age sex_encoded* infection_encoded* sofa time_to_rx charlson_score race_cat_collapsed* sampleweights* identityweights
	save slimmed_data_for_sim.dta , replace


*Writing simulation program for baseline prescribing practice. We'll "set" the coeffecient for antianaerobic antibiotics across a range of plausible effect sizes, and use the other coefficients from the model above to fill in the rest of the model. We'll then use this model to estimate the number of deaths potentially averted if all patients who received anti-anaerobic antibiotics received anaerobe-preserving antibiotics instead under these assumptions. Model to be simmed: logit(ninety_day_outcome)=𝛽0+𝛽1(zosyn_encoded)+𝛽2(age)+𝛽3(sex_encoded_2)+𝛽4(infection_encoded_2)+𝛽5(infection_encoded_3)+𝛽6(infection_encoded_4)+𝛽7(infection_encoded_5)+𝛽8(sofa)+𝛽9(time_to_rx)+𝛽10(charlson_score)+ 𝛽11(race_cat_collapsed_2) + 𝛽12(race_cat_collapsed_3) + 𝛽13(race_cat_collapsed_4).

*Base model:
	program simstudy , rclass
		syntax,  [ oddsratio(real 1.14962) ] // inputs and default values

		quietly { 
			drop _all // getting rid of old data for each sim
			use slimmed_data_for_sim.dta , replace // extra variables eliminated
			bsample 7569
		svyset _n [weight=sampleweights]
		generate xb = ((-5.177891) + ln(`oddsratio')*zosyn_encoded + (.0207687*age) + (-.0477247*sex_encoded_2) + (.3912146*infection_encoded_2) + (.0609579*infection_encoded_3) + (.0761246*infection_encoded_4) + (.5709092*infection_encoded_5) + (.260431*sofa) + (.0207238*time_to_rx) + (.1301484*charlson_score) + (-.0096951*race_cat_collapsed_2) + (.2356021*race_cat_collapsed_3) + (.3448969*race_cat_collapsed_4))
		generate predicted_outcome = rlogistic(xb,1) > 0 // assigning outcomes based on covariates and the model above
		svy: logistic predicted_outcome i.zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4 // for predict later
		predict prediction_astreated // predictions with "actual" treatment status
		replace zosyn_encoded = 0 // predictions if untreated
		predict prediction_unexposed // predictions if noone is treated
		
		generate te_ifnotreated = prediction_astreated - prediction_unexposed // calculating ATT
		svy : total te_ifnotreated // this gives me the lives saved if all treated were untreated
	}
		return scalar excessdeaths = r(table)[1,1] /*returning the number of lives saved and storing it as excessdeaths*/
	end

*Model assuming 20% more use of anti-anaerobic antibiotics:
	program simstudy_increased , rclass
		syntax,  [ oddsratio(real 1.14962) ] // inputs and default values

		{
			drop _all // getting rid of old data for each sim
			use slimmed_data_for_sim.dta , replace // extra variables eliminated
			bsample 7569
		svyset _n [weight=sampleweights_increased]
		generate xb = ((-5.177891) + ln(`oddsratio')*zosyn_encoded + (.0207687*age) + (-.0477247*sex_encoded_2) + (.3912146*infection_encoded_2) + (.0609579*infection_encoded_3) + (.0761246*infection_encoded_4) + (.5709092*infection_encoded_5) + (.260431*sofa) + (.0207238*time_to_rx) + (.1301484*charlson_score) + (-.0096951*race_cat_collapsed_2) + (.2356021*race_cat_collapsed_3) + (.3448969*race_cat_collapsed_4))
		generate predicted_outcome = rlogistic(xb,1) > 0 // assigning outcomes based on covariates and the model above
		svy: logistic predicted_outcome i.zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4 // for predict later
		predict prediction_astreated // predictions with "actual" treatment status
		replace zosyn_encoded = 0 // predictions if untreated
		predict prediction_unexposed // predictions if noone is treated
		
		generate te_ifnotreated = prediction_astreated - prediction_unexposed // calculating ATT
		svy : total te_ifnotreated // this gives me the lives saved if all treated were untreated
	}
		return scalar excessdeaths = r(table)[1,1] /*returning the number of lives saved and storing it as excessdeaths*/
	end

*Model assuming 20% less use of anti-anerobic antibiotics:
	program simstudy_decreased , rclass
		syntax,  [ oddsratio(real 1.14962) ] // inputs and default values

		quietly {
			drop _all // getting rid of old data for each sim
			use slimmed_data_for_sim.dta , replace // extra variables eliminated
			bsample 7569
		svyset _n [weight=sampleweights_decreased]
		generate xb = ((-5.177891) + ln(`oddsratio')*zosyn_encoded + (.0207687*age) + (-.0477247*sex_encoded_2) + (.3912146*infection_encoded_2) + (.0609579*infection_encoded_3) + (.0761246*infection_encoded_4) + (.5709092*infection_encoded_5) + (.260431*sofa) + (.0207238*time_to_rx) + (.1301484*charlson_score) + (-.0096951*race_cat_collapsed_2) + (.2356021*race_cat_collapsed_3) + (.3448969*race_cat_collapsed_4))
		generate predicted_outcome = rlogistic(xb,1) > 0 // assigning outcomes based on covariates and the model above
		svy: logistic predicted_outcome i.zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4 // for predict later
		predict prediction_astreated // predictions with "actual" treatment status
		replace zosyn_encoded = 0 // predictions if untreated
		predict prediction_unexposed // predictions if noone is treated
		
		generate te_ifnotreated = prediction_astreated - prediction_unexposed // calculating ATT
		svy : total te_ifnotreated // this gives me the lives saved if all treated were untreated
	}
		return scalar excessdeaths = r(table)[1,1] /*returning the number of lives saved and storing it as excessdeaths*/
	end



/* Simulate across all of the settings, and collect results in a table. We'll "set" the effect of anti-anaerobic antibipotics on mortality to a range from OR 0.9755 to OR 1.3275 in equal steps (these correspond to the range of risk differences observed in published human studies, converted using the prevalence of mortality in this study population and the standard formula). We'll conduct 1000 repetitions under each setting and report means and bootstrapped confidence intervals*/
	collect clear
	collect create tables2 , replace
	quietly {
		*First, the base model with the same level of anti-anaerobic prescribing as in the study population*
	
		foreach oddsratiosetting of numlist 0.9755 1.015 1.054 1.093 1.132 1.171 1.210 1.249 1.289 1.3275 {
			simulate excessdeaths = r(excessdeaths) , reps(1000) : simstudy, oddsratio(`oddsratiosetting')
			bootstrap mean=r(mean) , reps(1000) : summarize excessdeaths
			generate mean = r(table)[1,1]
			generate ci_lb = r(table)[5,1]
			generate ci_ub = r(table)[6,1]
			noisily: display "Treatment effect: " `oddsratiosetting' " " r(table)[1,1] " (" r(table)[5,1] ", " r(table)[6,1] ")"
			collect get scenario = 1 oddsratiovalue = `oddsratiosetting' teffect=r(table)[1,1] ci_lb=r(table)[5,1] ci_ub=r(table)[6,1]
			}
	
	collect layout (cmdset) (result[scenario] result[teffect] result[ci_lb] result[ci_ub] ) (), name(tables2)

	*Second, assuming 20% greater anti-anaerobic prescribing as in the study population*
		foreach oddsratiosetting of numlist 0.9755 1.015 1.054 1.093 1.132 1.171 1.210 1.249 1.289 1.3275 {
			simulate excessdeaths = r(excessdeaths) , reps(1000) : simstudy_increased, oddsratio(`oddsratiosetting')
			bootstrap mean=r(mean) , reps(1000) : summarize excessdeaths
			generate mean = r(table)[1,1]
			generate ci_lb = r(table)[5,1]
			generate ci_ub = r(table)[6,1]
			noisily: display "Treatment effect: " `oddsratiosetting' " " r(table)[1,1] " (" r(table)[5,1] ", " r(table)[6,1] ")"
			collect get scenario = 2 oddsratiovalue = `oddsratiosetting' teffect=r(table)[1,1] ci_lb=r(table)[5,1] ci_ub=r(table)[6,1]
			}
	
	collect layout (cmdset) (result[scenario] result[teffect] result[ci_lb] result[ci_ub] ) (), name(tables2)

	*Second, assuming 20% less anti-anaerobic prescribing as in the study population*
		foreach oddsratiosetting of numlist 0.9755 1.015 1.054 1.093 1.132 1.171 1.210 1.249 1.289 1.3275 {
			simulate excessdeaths = r(excessdeaths) , reps(1000) : simstudy_decreased, oddsratio(`oddsratiosetting')
			bootstrap mean=r(mean) , reps(1000) : summarize excessdeaths
			generate mean = r(table)[1,1]
			generate ci_lb = r(table)[5,1]
			generate ci_ub = r(table)[6,1]
			noisily: display "Treatment effect: " `oddsratiosetting' " " r(table)[1,1] " (" r(table)[5,1] ", " r(table)[6,1] ")"
			collect get scenario = 3 oddsratiovalue = `oddsratiosetting' teffect=r(table)[1,1] ci_lb=r(table)[5,1] ci_ub=r(table)[6,1]
			}
			
	collect layout (cmdset) (result[scenario] result[oddsratiovalue] result[teffect] result[ci_lb] result[ci_ub] ) (), name(tables2)
	}

*We'll save the simulation results in a table for later viewing/reporting
	collect export simoutputdata.xlsx , name(tables2) replace
	import excel simoutputdata.xlsx, firstrow clear

*We'll next convert the odds ratios to risk differences for clarity of reporting
	generate ARI = -1*(0.205 - ((0.205 / (1 - 0.205)) * oddsratiovalue) / (1 + ((0.205 / (1 - 0.205)) * oddsratiovalue)))

*We'll generate the graph using the simulation results above and the absolute risk increases.
	graph twoway ///
	(line teffect ARI if scenario == 1, lcolor(navy) lpattern(solid)) ///
	(rarea ci_lb ci_ub ARI if scenario == 1, fcolor(navy%50) lcolor(navy%50)) ///
	(line teffect ARI if scenario == 2, lcolor(red) lpattern(solid)) ///
	(rarea ci_lb ci_ub ARI if scenario == 2, fcolor(red%50) lcolor(red%50)) ///
	(line teffect ARI if scenario == 3, lcolor(green) lpattern(solid)) ///
	(rarea ci_lb ci_ub ARI if scenario == 3, fcolor(green%50) lcolor(green%50)) ///
	, ytitle("Excess Deaths per Year") xtitle("Treatment Effect" "(Absolute Risk Difference, %)") ///
	xlabel(-0.005 "-0.5" 0.00 "0" 0.005 "0.5" 0.01 "1" 0.015 "1.5" 0.02 "2" 0.025 "2.5" 0.03 "3" 0.035 "3.5" 0.04 "4" 0.045 "4.5" 0.05 "5" , labsize(small)) ylabel(0(5000)40000, labsize(small)) ///
	legend(order(1 "Study population" "antibiotic practices" "" 3 "Anti-anaerobic antibiotics" "used 20% more often" 5 "Anti-anerobic antibiotics" "used 20% less often")) name("AttributableDeaths", replace) ///
	xline(0.03368, lpattern(dash) lcolor(black)) xline(0.039, lpattern(dash) lcolor(black)) xline(0.05, lpattern(dash) lcolor(black))

*We'll export the graph as a jpeg
graph export AttributableDeaths.jpg , as(jpg) replace

*Finally, we generate some summary data (e.g., frequency of antianaerobic antibiotics at baseline) and both unadjusted and adjusted treatment effect data for inclusion in the manuscript. We don't really use this data in the sim but report it for completeness.
	*Paper Text
	use slimmed_data_for_sim.dta , replace
	sum ninety_day_outcome
	tab zosyn_encoded
	teffects ipw (ninety_day_outcome) (zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4)

	teffects ipw (ninety_day_outcome) (zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4)

	logistic ninety_day_outcome i.zosyn_encoded
	logistic ninety_day_outcome i.zosyn_encoded c.age i.sex_encoded_2 i.infection_encoded_2 i.infection_encoded_3 i.infection_encoded_4 i.infection_encoded_5 c.sofa c.time_to_rx c.charlson_score i.race_cat_collapsed_2 i.race_cat_collapsed_3 i.race_cat_collapsed_4
	margins i.zosyn_encoded
		predict prediction_astreated // predictions with "actual" treatment status
		replace zosyn_encoded = 0 // predictions if untreated
		predict prediction_unexposed // predictions if noone is treated
		generate te_ifnotreated = prediction_astreated - prediction_unexposed // calculating ATT
		total te_ifnotreated
