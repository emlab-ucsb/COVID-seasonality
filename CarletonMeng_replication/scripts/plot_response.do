
//Description: This code makes figures from regressions run in country_regression.do

/********************************   HEADER ***********************************************/
	clear all
	set more off
	set trace off
	set tracedepth 1
	set matsize 11000 //10000
	set maxvar 32000
	set scheme s1color
	
	//set file paths
	global startDir "/Users/kylemeng/Dropbox/work/research/COVID/repo/replication"  
	global dataDir "$startDir/data"
	global tablesDir "$startDir/tables"
	global figuresDir "$startDir/figures"
	global tempDir "$startDir/temp"



/******************************** Load data ***********************************************/

//Load data & generate needed variables
	use $dataDir/country_date_COVID_temperature, replace
		
//Prepare variables
	gen D_cases_rate = D_cases/(pop2018/1e6) // cases per million
	egen cID=group(iso3)

	gen time1=date-21914
	gen time2=time1^2
	xtset cID date

//create lat bins 
	local lMin=-50
	local lInt=25
	local lMax=50
    egen latCat=cut(lat), at(`lMin'(`lInt')`lMax')
    replace latCat=`lMin' if lat<=`lMin'
    replace latCat=`lMax' if lat>=`lMax'
	
//create controls
	xi, pref(_t) noomit i.time1
	xi, pref(_ct) noomit i.cID*time1
	xi, pref(_lt) noomit i.latCat*time1
	xi, pref(_lt2) noomit i.latCat*time2


/******************************** Load estimates & predict ***********************************************/
	
* set up plotting data
*set parameters of the plot (min,max,"ommited temp")		
local min = -10
local max = 35
local obs = `max' - `min' + 1
local omit = 20
	
* predict 
drop if _n > 0
set obs `obs'
replace tavg_poly_1_pop = _n + `min' - 1

* loop over models you want to show results for	
foreach mod in "lin_dFE_noclim" "lin_dFE"  "quad_dFE" "lin_dFE_latbandtr2_noclim" "lin_dFE_latbandtr2" "quad_dFE_latbandtr2" "lin_dFE_5deglatbandtr2" { //"lin_dFE_conttr" "lin_dFE_latbandtr" "lin_dFE"  "quad_dFE"  "cubic_dFE" "quart_dFE" 
	
	estimates use "$tempDir/poisson_`mod'.ster"

	* nonlinear
	if regexm("`mod'","quad")==1 {
		local line = "_b[tavg_poly_1_pop] * (tavg_poly_1_pop-`omit') + _b[tavg_poly_2_pop]*(tavg_poly_1_pop^2 -`omit'^2)"
	} 
	else if regexm("`mod'","cubic")==1 {
		local line = "_b[tavg_poly_1_pop] * (tavg_poly_1_pop-`omit') + _b[tavg_poly_2_pop]*(tavg_poly_1_pop^2=`omit'^2) + _b[tavg_poly_3_pop]*(tavg_poly_1_pop^3 - `omit'^3)"
	} 
	else if regexm("`mod'","quart")==1 {
		local line = "_b[tavg_poly_1_pop] * (tavg_poly_1_pop-`omit') + _b[tavg_poly_2_pop]*(tavg_poly_1_pop^2=`omit'^2) + _b[tavg_poly_3_pop]*(tavg_poly_1_pop^3 - `omit'^3)+ _b[tavg_poly_4_pop]*(tavg_poly_1_pop^3 - `omit'^4)"
	} 
	
	* linear
	else {
	local line = "_b[tavg_poly_1_pop] * (tavg_poly_1_pop - `omit')"			
	}
	
	predictnl yhat_`mod' = `line', se(se_`mod') ci(lbci_`mod' ubci_`mod')
	}

/******************************** Plot ***********************************************/

	
* treating lat bands as main model, others plotted as robustness	
tw rarea ubci_lin_dFE_latbandtr2 lbci_lin_dFE_latbandtr2  tavg_poly_1_pop, fcolor(maroon%30) lcolor(white) || ///
	line yhat_lin_dFE_latbandtr2 tavg_poly_1_pop, color(maroon) lwidth(.8) || ///
	line yhat_lin_dFE_noclim tavg_poly_1_pop, color(blue*1.25) lpattern(longdash) || ///
	line yhat_lin_dFE tavg_poly_1_pop, color(ebblue*1.25) lpattern(dash_dot) || ///
	line yhat_quad_dFE tavg_poly_1_pop, color(purple*1.25) lpattern(shortdash_dot) || ///
	line yhat_lin_dFE_latbandtr2_noclim tavg_poly_1_pop, color(red*1.25) lpattern(longdash) || ///
	line yhat_quad_dFE_latbandtr2 tavg_poly_1_pop, color(pink*1.25) lpattern(dash_dot) || ///
	line yhat_lin_dFE_5deglatbandtr2 tavg_poly_1_pop, color(orange*1.25) lpattern(shortdash_dot) || ///
	, yline(0, lcolor(gs12) lpattern(shortdash)) xlabel(-10(5)35, nogrid) ///
	ylabel(,nogrid) ///
	ytitle("log(new COVID-19 cases per 1 million)") xtitle("Population weighted daily temperature (C)") ///
	legend(	label(1 "95% CI") ///
			label(2 "lat. trends") ///
			label(3 "No cntrls") ///
			label(4 "day FE") ///
			label(5 "poly2") ///
			label(5 "poly2") ///
			label(6 "lat. trends no cntrls") ///
			label(7 "lat. trends poly2") ///
			label(8 "5deg lat. trends") ///
			ring(0) position(2) col(2) bmargin(zero)) /// 
	scheme(plotplain)
	
 graph export "$figuresDir/temp_response_allspecs.pdf", replace	
 

 //Histogram 
	use $dataDir/country_date_COVID_temperature, replace
	
	bysort iso3: egen mean_tavg = mean(tavg_poly_1_pop)
	
	foreach cc in "USA" "CHN" "KOR" "ITA" "AUS" "IRN" "ZAF" {
		qui summ mean_tavg if iso3=="`cc'"
		loc m`cc' = `r(mean)'
		di "`cc' = " `m`cc''
		}
		
		// histogram to put under the plot later
		tw hist tavg_poly_1_pop, color(gs12) xline(`mUSA' `mKOR' `mCHN' `mITA' `mAUS', lcolor(emerald)) ///
		xlabel(-20(10)40  `mKOR' "KOR"  `mCHN' "CHN" `mITA' "ITA" `mAUS' "AUS") xtitle("Avg. temperature Jan 1, 2020 - March 15, 2020 (C)")
		graph export "$figuresDir/insample_tavg_2020.pdf", replace
	}
	
 }
