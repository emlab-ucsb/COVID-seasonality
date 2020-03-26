

//Description: This code processes ERA weather data
/********************************   HEADER ***********************************************/
	clear all
	set more off
	set trace on
	set tracedepth 1
	set matsize 11000 //10000
	set maxvar 32000

	//set file paths
	global startDir "/Users/kylemeng/Dropbox/work/research/COVID/repo/replication"  
	global dataDir "$startDir/data"
	global tablesDir "$startDir/tables"
	global figuresDir "$startDir/figures"
	global tempDir "$startDir/temp"

	
	//Parameters
	local b=-.1285677  
	local blb=-.2130995  
	local bub=-.0440359
	local censor=500	

/******************************** MAIN ***********************************************/


//Figure 1: country scatter pct change by July 2020

	use $dataDir/country_month_ERA_temperature, replace

	rename iso ISO3 
	merge ISO3 using $dataDir/country_db, unique sort
	tab _merge //all merged
	rename LAT lat
	rename ISO3 iso
	drop FIPS ISO2 UN NAME AREA POP2005 REGION SUBREGION LON
	drop _merge

	drop if name=="Antarctica"  //no people
	drop if name=="Mongolia" //really cold
	drop if name=="Kazakhstan" //cold 

	//Construct estimates
	gen DT=avg_m07T-avg_m03T
	gen effect=(exp(`b'*(DT))-1)*100
	gen effect_lb=(exp(`blb'*(DT))-1)*100
	gen effect_ub=(exp(`bub'*(DT))-1)*100

	//get means
	//North 
	sum effect if lat>=0, d
	loc xnorth_jul = `r(mean)'
	
	//South
	qui sum effect if lat<0
	loc xsouth_jul = `r(mean)'

	//top 20 countries
	gen top_country=name
	replace top_country="" if ///
	iso!=("CHN") & ///
	iso!=("ITA") & ///
	iso!=("IRN") & ///
	iso!=("KOR") & ///
	iso!=("ESP") & ///
	iso!=("FRA") & ///
	iso!=("DEU") & ///
	iso!=("USA") & ///
	iso!=("CHE") & ///	
	iso!=("NOR") & ///
	iso!=("SWE") & ///
	iso!=("NLD") & ///
	iso!=("DNK") & ///
	iso!=("GBR") & ///
	iso!=("JPN") & ///
	iso!=("BEL") & ///
	iso!=("AUT") & ///
	iso!=("QAT") & ///
	iso!=("SGP") & ///
	iso!=("AUS") 

	replace effect_lb=`censor' if effect_lb>`censor' 	//truncate bound

	//Figure 
	tw ///
	(rspike effect_lb effect_ub lat, horizontal lwidth(.2) lcolor(gs10)) ///
	(scatter lat effect if top_country=="", mcolor(black) msymbol(smcircle) msize(tiny) ///
	mlabel(top_country) mlabsize(vsmall) mlabp(3) mlabgap(45)) ///
	(scatter lat effect if top_country!="", mcolor(cranberry) msymbol(smcircle) msize(tiny)) ///
	, ///
	xline(0, lpattern(dash) lwidth(.2) lcolor(gs13)) ///
	xline(`xnorth_jul' `xsouth_jul', lpattern(solid) lwidth(.2) lcolor(eltblue)) ///
	ylabel(-60(20)80, nogrid) ytitle("Latitude") ///
	xlabel(-100(50)`censor', nogrid) 	///
	xtitle("Projected percent change in new 2019-nCoV cases from March 2020 to July 2020") ///
	legend(off) ///
	name(summer)	

	graph export $figuresDir/fig_projection_jul_by_country.pdf, as(pdf) replace


//Figure 1: country scatter pct change by January 2021

	use $dataDir/country_month_ERA_temperature, replace

	rename iso ISO3 
	merge ISO3 using $dataDir/country_db, unique sort
	tab _merge //all merged
	rename LAT lat
	rename ISO3 iso
	drop FIPS ISO2 UN NAME AREA POP2005 REGION SUBREGION LON
	drop if name=="Antarctica" // no people
	drop if name=="Mongolia" //really cold
	drop if name=="Kazakhstan" //cold 

//Construct estimates
	gen DT=avg_m01T-avg_m03T
	gen effect=(exp(`b'*(DT))-1)*100
	gen effect_lb=(exp(`blb'*(DT))-1)*100
	gen effect_ub=(exp(`bub'*(DT))-1)*100

	//get means
	//North 
	sum effect if lat>=0, d
	loc xnorth_jan = `r(mean)'
	
	//South
	qui sum effect if lat<0
	loc xsouth_jan = `r(mean)'

	//top 20 countries
	gen top_country=name
	replace top_country="" if ///
	iso!=("CHN") & ///
	iso!=("ITA") & ///
	iso!=("IRN") & ///
	iso!=("KOR") & ///
	iso!=("ESP") & ///
	iso!=("FRA") & ///
	iso!=("DEU") & ///
	iso!=("USA") & ///
	iso!=("CHE") & ///	
	iso!=("NOR") & ///
	iso!=("SWE") & ///
	iso!=("NLD") & ///
	iso!=("DNK") & ///
	iso!=("GBR") & ///
	iso!=("JPN") & ///
	iso!=("BEL") & ///
	iso!=("AUT") & ///
	iso!=("QAT") & ///
	iso!=("SGP") & ///
	iso!=("AUS") 

	replace effect_lb=`censor' if effect_lb>`censor' 	//truncate bound

//Figure 
	tw ///
	(rspike effect_ub effect_lb lat, horizontal lwidth(.2) lcolor(gs10)) ///
	(scatter lat effect if top_country=="", mcolor(black) msymbol(smcircle) msize(vsmall) ///
	mlabel(top_country) mlabsize(vsmall) mlabp(3) mlabgap(-35)) ///
	(scatter lat effect if top_country!="", mcolor(cranberry) msymbol(smcircle) msize(tiny)) ///
	, ///
	xline(0, lpattern(dash) lwidth(.2) lcolor(gs13)) ///
	xline(`xnorth_jan' `xsouth_jan', lpattern(solid) lwidth(.2) lcolor(eltblue)) ///
	ylabel(-60(20)80, nogrid) ytitle("Latitude") ///
	xlabel(-100(50)`censor', nogrid) ///
	xtitle("Projected percent change in new 2019-nCoV cases from March 2020 to January 2021") ///
	legend(off) ///
	name(winter)

	graph export $figuresDir/fig_projection_january_by_country.pdf, as(pdf) replace


//Figure 3: spaghetti country-level pct change plot

	use $dataDir/country_month_ERA_temperature, replace

	reshape long avg_m, i(iso) j(month, string)

	replace month=subinstr(month, "T", " ",.)
	destring month, force replace
	rename avg_m T

	egen cID=group(iso3)

	//April temperature
	bys cID: gen T3=T if month==3
	gsort cID -T3
	bys cID: replace T3=T3[1]
	
	//Get change in temperature
	sort cID month
	gen DT=0
	forvalues m=1(1)12 {
		bys cID: replace DT=T-T3 if m==`m'
	}

	//get estimate
	gen effect=(exp(`b'*(DT))-1)*100

	rename iso ISO3 
	merge ISO3 using $dataDir/country_db, uniqusing sort
	tab _merge //all merged
	rename LAT lat
	rename ISO3 iso
	drop FIPS ISO2 UN NAME AREA POP2005 REGION SUBREGION LON
	drop _merge

	//drop some countries because they're outliers
	drop if name=="Antarctica" // no people 
	drop if name=="Mongolia" //really cold
	drop if name=="Kazakhstan" //cold 
	
	//change xaxis
	gen new_month=month-2
	replace new_month=11 if new_month==-1 //january
	replace new_month=12 if new_month==0 //Feb
	drop month
	rename new_month month 

	sort cID month

	//get hemisphere mean
	preserve
	keep if lat>=0
	collapse (mean) effect, by(month)
	gen iso="northern"
	tempfile temp
	save "`temp'"
	restore

	merge iso month using "`temp'", uniqusing sort
	drop _merge

	preserve
	keep if lat<0
	collapse (mean) effect, by(month)
	gen iso="southern"
	tempfile temp
	save "`temp'"
	restore

	merge iso month using "`temp'", uniqusing sort
	drop _merge

//Figure 
	local width=.15
	local width_mean=.8
	local color_north "orange"
	local color_south "eltblue"
	local color_north_mean "red"
	local color_south_mean "ebblue"
	local trans=40
	local plotCommand ""

	//Northern hemisphere
	preserve
	keep if lat>=0
	levelsof cID, local(cList)
	local i=1
	foreach c in `cList' {
		local plotCommand "`plotCommand'(line effect month if cID==`c', lwidth(`width') lcolor(`color_north'%`trans') lpattern(solid))"
		local i=`i'+1
	}
	restore 
	//mean
	local plotCommand "`plotCommand'(line effect month if iso=="northern", lwidth(`width_mean') lcolor(`color_north_mean') lpattern(solid))"


	//Southern hemisphere
	preserve
	keep if lat<0
	levelsof cID, local(cList)
	local i=1
	foreach c in `cList' {
		local plotCommand "`plotCommand'(line effect month if cID==`c', lwidth(`width') lcolor(`color_south'%`trans') lpattern(solid))"
		local i=`i'+1
	}
	restore 
	//mean
	local plotCommand "`plotCommand'(line effect month if iso=="southern", lwidth(`width_mean') lcolor(`color_south_mean') lpattern(solid))"


	tw `plotCommand' ///
	, ///
	text(290 1.5 "2020") ///
	text(290 11.5 "2021") ///
	xlabel(1 "Mar" 2 "Apr" 3 "May" 4 "Jun" 5 "Jul" 6 "Aug" 7 "Sep" 8 "Oct" 9 "Nov" 10 "Dec" 11 "Jan" 12 "Feb", nogrid) ///
	ylabel(, nogrid) ///
	xtitle("Month") ///
	ytitle("Projected percent change in new 2019-nCoV cases") ///
	legend(off) 

	graph export $figuresDir/fig_projection_spaghetti_by_country.pdf, as(pdf) replace


