

//Description: This code creates main table 

/********************************   HEADER ***********************************************/
	clear all
	set more off
	set trace off
	set tracedepth 1
	set matsize 11000 //10000
	set maxvar 32000

	//set file paths
	global startDir "/Users/kylemeng/Dropbox/work/research/COVID/repo/replication"  
	global dataDir "$startDir/data"
	global tablesDir "$startDir/tables"
	global figuresDir "$startDir/figures"
	global tempDir "$startDir/temp"

	
/******************************** MAIN ***********************************************/

//Load data
	//insheet using $dataDir/country_date_COVID_temperature.csv, comma names clear
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

	label var tavg_poly_1_pop "Temperature"
	label var tavg_poly_2_pop "Temperature squared"
	label var prcp_poly_1_pop "Precipitation"
	label var humd_poly_1_pop "Specific humidity"

	drop if D_cases_rate<0 //only 3 observations, JHU team noted this was due to coding error

		
*********** Main table *********** 
	loc outfile = "$tablesDir/main.tex"
	
	*********** COLS 1:3 NO SPACExTIME CONTROLS ***********
	//2-way fixed effect: linear w/o climate controls
	xi: xtpoisson D_cases_rate tavg_poly_1_pop _ttime1*, fe robust
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	replace addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, No, Lat. bin trends, No, Lat. bin width, "-") ///
	drop(*ttime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_lin_dFE_noclim.ster", replace
	
	//2-way fixed effect: linear
	xi: xtpoisson D_cases_rate tavg_poly_1_pop prcp_poly_1_pop humd_poly_1_pop _ttime1*, fe robust
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, Yes, Lat. bin trends, No, Lat. bin width, "-") ///
	drop(*ttime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_lin_dFE.ster", replace
	
	//2-way fixed effect: quadratic 
	xi: xtpoisson D_cases_rate tavg_poly_1_pop tavg_poly_2_pop prcp_poly_1_pop humd_poly_1_pop _ttime1*, fe robust
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, Yes, Lat. bin trends, No, Lat. bin width, "-") ///
	drop(*ttime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_quad_dFE.ster", replace
	
	*********** COLS 4:6 LINEAR 5deg LAT BAND TRENDS ***********
	//2-way fixed effect plus latitude ban-level trends: linear w/o climate controls
	xi: xtpoisson D_cases_rate tavg_poly_1_pop  _ttime1* _ltlatXtime* _lt2latXtime*, fe robust  //_lt2latXtime*
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, No, Lat. bin trends, Quad, Lat. bin width, "`lInt'\$ ^\circ\$") ///
	drop(*ttime* *latXtime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_lin_dFE_latbandtr2_noclim.ster", replace
	
	//2-way fixed effect plus latitude ban-level quad trends: linear
	xi: xtpoisson D_cases_rate tavg_poly_1_pop  prcp_poly_1_pop humd_poly_1_pop _ttime1* _ltlatXtime* _lt2latXtime*, fe robust //_lt2latXtime*
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, Yes, Lat. bin trends, Quad, Lat. bin width, "`lInt'\$ ^\circ\$") ///
	drop(*ttime* *latXtime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_lin_dFE_latbandtr2.ster", replace


	//2-way fixed effect plus latitude ban-level quad trends: quadratic
	xi: xtpoisson D_cases_rate tavg_poly_1_pop tavg_poly_2_pop  prcp_poly_1_pop humd_poly_1_pop _ttime1* _ltlatXtime* _lt2latXtime*, fe robust
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, Yes, Lat. bin trends, Quad, Lat. bin width, "`lInt'\$ ^\circ\$") ///
	drop(*ttime* *latXtime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_quad_dFE_latbandtr2.ster", replace
	
	*********** COL 7 Narrower bandwidth ***********
	
	// need to remake lat bands
	cap drop latCat _ltlatXtime* _lt2latXtime*
	local lMin=-60
	local lInt=5
	local lMax=60
    egen latCat=cut(lat), at(`lMin'(`lInt')`lMax')
    replace latCat=`lMin' if lat<=`lMin'
    replace latCat=`lMax' if lat>=`lMax'
	xi, pref(_lt) noomit i.latCat*time1
	xi, pref(_lt2) noomit i.latCat*time1
		
	// climate controls, linear
	xi: xtpoisson D_cases_rate tavg_poly_1_pop prcp_poly_1_pop humd_poly_1_pop _ttime1* _ltlatXtime* _lt2latXtime* , fe robust // _lt2latXtime* 
	count if e(sample)
	local oNum: display %8.0f r(N)
	unique cID if e(sample)
	local cNum=r(sum)
	unique date if e(sample)
	local dNum=r(sum)
	outreg2 using "`outfile'", dec(3) label tex(frag) ///
	append addtext(Observations, `oNum', Countries, `cNum', Days, `dNum', ///
	Weather controls, Yes, Lat. bin trends, Quad, Lat. bin width, "`lInt'\$ ^\circ\$") ///
	drop(*ttime* *latXtime*) ctitle(" ") nonotes noobs noni
	estimates save "$tempDir/poisson_lin_dFE_5deglatbandtr2.ster", replace

	
	shell sed -i -e 's/VARIABLES \& /\&/' "`outfile'"
	rm "`outfile'-e"


