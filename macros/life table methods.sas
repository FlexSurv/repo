/*	
	life table methods (dev)
	(in development)

	5 Nov 2021

	11 Jan 2022
	added option for calculation of CIs that allows chosing either loghaz 
	or log transform (see CI parameter)

*/

*	rel_surv   	a macro to be used in the life table analysis of relative survival
				based on a much earlier version written by Dr. Paul Lambert, with 
				additions by Larry Elison (Statistics Canada)
				this version was developed with support from the Canadian Partnership
				Against Cancer (CPAC);
%macro rel_surv(
	infile = ,			 	/*	input dataset (may contain a ^where^ clause)			*/
	weight_lib = work,		/*	name of sas library containing weights					*/
	weight_var =,			/*	name of variable holding weights for standardisation	*/				
	standstrata =,			/*	name of variable to stratify for standardisation		*/
							/*	weight_var and stand_by are jointly sufficient for 		*/
							/*	and over-ride any other standardisation specification	*/
	age = age_yrs,			/*	age at diagnosis in completed years						*/
	sex = sex,				/*	sex variable (coding must agree with life table)		*/
	origin = ,				/*	date of diagnosis 										*/
	exit = 	,				/* 	date of death/censoring									*/
	entry = ,				/*	date of left truncation if required						*/
	yydx = ,				/*	year of diagnosis (required if survival duration 		*/
							/*	is supplied as death date variable.)					*/	
	scale = 365.24,			/*	to convert duration in days to duration in years		*/
	censor = censor(1),		/*	censor indicator variable, with value indicating 		*/
							/*	  censored observation eg. censor = c(1) in the case	*/
							/*	  where the variable c is 1 for a censored observation	*/
							/*	  can also be a comma-separated list of values, eg 0,4	*/
	stnd = ,				/*	weighting to be used for standardisation, if required 	*/
							/*	  use ICSS() or CDN()									*/
							/*	  ICSS takes values in (1,2,3,4,9)						*/
							/*	  CDN takes values in()									*/	
	intervals = ,			/*	life table intervals to be used in analysis				*/
							/*	monthly to 6 months, then every 3 months to 2 years		*/
							/*	these defaults are for the CPAC 2-yr analysis			*/
							/*	the specification intervals = failures is also allowed	*/
	popmort = stat.life_tab2014,
							/*	life table for your province							*/
	survprob = prob,		/*	variable in life table holding survival probability 	*/	
	ci = loghaz,			/*	loghaz for log cumulative excess hazard 				*/
							/*	that is, log(-log(R(t)))								*/
							/*	haz for calculation on the cumulative excess hazard 	*/
							/*	scale that is, log(R(t)	 								*/	
	mergeby = _age sex _year,
							/*	life table variables must be in correct order			*/ 
							/*	(age/sex/year) followed by any other stratifier,		*/
							/*	(eg, SES, province).  sql code will be generated for	*/
							/*	any other stratifiers 									*/
	maxage = 99,			/*	maximum age available in life tables					*/
							/*	subjects who age past this age will have their attained	*/
							/*	age fixed at this value									*/
	strat = ,				/*	stratification variable(s) 				 				*/								
	strat_fmt = ,			/*	formats for grouping stratification variables 			*/
							/*	  as a string of variable/format pairs					*/
	PatientID = ,			/*	to check for duplicated patients, and for debugging		*/
							/*	also needed for Pohar-Perme weights						*/
	crude = 0,				/*	set to 1 to see all survival estimates					*/
							/*	or to an expression that evaluates to 0 (false) or		*/
							/*	non-zero (true) that could be evaluated in a ^where^ clause	*/
	age_adj =,				/*	will be set to 0 if no age-standardisation.				*/
							/*	if specified in an age-standardisation run, will be		*/
							/*	retained unmodified, otherwise set to 1 (print all)		*/
	crude_estimates			/*	name of file holding crude estimates 					*/
			= grouped,		/*	default value is grouped								*/
	list = ederer,			/*	variables to output in the crude estimates step			*/
							/*	chose from 												*/
							/*	ederer  (Ederer 2, the default)							*/
							/*	pohar													*/
							/*	cuminc (cumulative probabilities of death)				*/
							/*	age-standardisation is only available for ederer and pohar	*/
							/*	stratifiers, interval label and numbers of cases		*/
							/*	are always listed.  Custom list is also a possibility	*/
	std_estimates 			/*	name of file holding age-standardised estimates			*/
			= std_final);	/*	default value is std_final								*/

*	-----------------------------------------;
*	begin error checking of passed parameters;
*	-----------------------------------------;

*	required paramaters:
	case dataset must exist and contain variables passed in 
		patientID, age, sex, diagnosis date (origin), death date (exit) , censoring
		entry (if present).  

	population life table must exist and contain variablees sex, _age, _year, prob
	intervals parameter may not be null and must start at 0 (because of Pohar)
	
	standardisation may be based on one of the following mutually-exclusive methods)
		a) specified variable and standardisation stratifier (using weight_var and standstrata)
		b) requesting ICSS or Canadian weights (stnd parameter)
		method a) takes precedence over method b);
		

*	failure on these tests causes immediate end of error processing;	
	%if "&infile." eq  "" %then %do;
		%err_mess( Infile parameter is required ); 
		%goto fin;
	%end;
	%let filename = %scan(&infile.,1,'(');
	%if not %sysfunc(exist(&filename.)) %then %do;
		%err_mess( case file &filename. not found ); 
		%goto fin;
	%end;
	%if "&popmort." eq "" %then %do;
		%err_mess( population life table parameter is required); 
		%goto fin;
	%end;
	%if not %sysfunc(exist(&popmort.)) %then %do;
		%err_mess( life table file &popmort. not found  ); 
		%goto fin;
	%end;

*	must be at least 3 variables in the mergeby parameter;
	%let wc = %sysfunc(countw(&mergeby.,' '));
	%if &wc. lt 3 %then %do;
		%err_mess( at least 3 variables must be specified (age/sex/year) ); 
		%goto fin;
	%end;
	
*	either one or none of method a) and b) above may be specified;
	%if "&weight_var." ^= "" and "&stnd." ^= "" %then %do;
		%err_mess(supply only one of weight_var or stnd);
		%goto fin;
	%end;
	%if &weight_var ^= and &standstrata = %then %do;
		%err_mess(standstrata parameter is required for use of weight_var parameter);
		%goto fin;
	%end;

	%if &yydx. ^= and &exit. = %then %do;
		%err_mess(yydx must be accompanied by the exit paramter);
		%goto fin;
	%end;

	%let ci = %lowcase(&ci.);
	%if %sysfunc(findw(loghaz haz, &ci.)) eq 0 %then %do;
		%err_mess(CI parameter accepts only loghaz or haz. You specified &ci.);
		%goto fin;
	%end;
	
	%if &exit. ^= and &yydx. = and &origin. = %then %do;
		%err_mess(if specifying dates, both origin and exit dates are required (do not specify yydx parameter));
		%goto fin;
	%end;

	%let err = 0;
*	all other potential errors are processed independently;

*	check for required parameters;
	%if &sex. eq  %then %err_mess(sex parameter is required ); 
	%if &age. eq  %then %err_mess(age parameter is required); 
	%if &exit. eq  %then %err_mess(exit parameter is required);
	%if "&censor." eq "" %then %err_mess(censoring parameter is required ); 
	%let cens_var = %scan(&censor.,1,'(');
	%let cens_val = %scan(&censor.,2,'()');
	%if &cens_val. eq  %then %err_mess(censoring value is required); 
	%if &intervals. eq  %then %err_mess(intervals parameter is required); 
	%if &intervals. ne  and %lowcase(&intervals.) ne failures %then %do;
		%if %scan(&intervals., 1) ne 0 %then %err_mess(intervals must start at 0); 
	%end;
	%if "&PatientID." eq "" %then %err_mess(patient ID parameter is required); 

*	check life table file for required variables;
	%let dsid=%sysfunc(open(&popmort.));
	%do c = 1 %to &wc.;
		%if not %sysfunc(varnum(&dsid., %scan(&mergeby.,&c.))) 
			%then %err_mess(%str(variable %scan(&mergeby.,&c.) is not in &popmort.  )); 			
		%let merge&c. = %scan(&mergeby.,&c.);
	%end;
	%if not %sysfunc(varnum(&dsid., &survprob.)) %then %err_mess(variable &survprob. is not in &popmort. ); 
	%let rc = %sysfunc(close(&dsid.));
	
* 	Check for empty case file;
	data _null_;
  		if eof then do;
  		 	call symput('err',1);
  		 	put "ERROR: *** Error in paramater specification ***";
     		put "ERROR: specified data file is empty";
 		end;
 		stop;
  		set &infile. end=eof;
	run;
	%if &err. %then %goto fin;
	
*	check case file for specified variables;
	%let dsid=%sysfunc(open(&filename.));
	%if not %sysfunc(varnum(&dsid.,&sex.)) 
		%then %err_mess(variable &sex. specified by sex parameter is not in &filename.); 
 	%if "&yydx." ne "" %then %do; 
		%if not %sysfunc(varnum(&dsid.,&yydx. ))  
 			%then %err_mess(variable &yydx. specified by yydx parameter is not in &filename.);  
	%end;
	
	%if "&PatientID." ne "" %then %do;
		%if not %sysfunc(varnum(&dsid.,&PatientID.)) 			 
 			%then %err_mess(variable name &PatientID. specified by PatientID parameter is not in &filename.); 
	%end;

 	%if "&entry." ne "" %then %do;
		%if not %sysfunc(varnum(&dsid.,&entry.)) 			 
 			%then %err_mess(variable name &entry. specified by entry parameter is not in &filename.); 
	%end; 
 
 	%if "&origin." ne "" %then %do;
		%if not %sysfunc(varnum(&dsid.,&origin.)) 			 
 			%then %err_mess(variable name &origin. specified by origin parameter is not in &filename.); 
	%end; 
 
 	%if "&weight_var." ne "" %then %do;
		%if not %sysfunc(varnum(&dsid.,&weight_var.)) 			 
 			%then %err_mess(variable name &weight_var. specified by weight_var parameter is not in &filename.); 
	%end; 

 	%if "&standstrata." ne "" %then %do;
		%if not %sysfunc(varnum(&dsid.,&standstrata.)) 			 
 			%then %err_mess(variable name &standstrata. specified by standstrata parameter is not in &filename.); 
	%end;

	%if not %sysfunc(varnum(&dsid.,&age.)) 
		%then %err_mess(variable &age. specified by age parameter is not in &filename.); 
	%if not %sysfunc(varnum(&dsid.,&exit.)) 
		%then %err_mess(variable &exit. is not in &filename.); 
	%if not %sysfunc(varnum(&dsid.,&cens_var.)) 
		%then %err_mess(variable name &cens_var. specified by censoring parameter is not in &filename. ); 
	%if &wc. gt 3 %then %do c = 4 %to &wc.;
		%if not %sysfunc(varnum(&dsid., %scan(&mergeby.,&c.))) 			
			%then %err_mess(variable name &merge&c. is not in &filename.); 
	%end;
	
	%let rc = %sysfunc(close(&dsid.));
	
*	standardisation choices (Canadian or international weights);
*	only certain values for standard group are acceptable, depending on the method chosen;
*	also, check that age is not an explicit stratifier;
*	check for weight files at this point;
	%if "&weight_var" = "" %then %do;
		%if "&stnd." = "" %then %let age_adj  = 0;
		%if "&stnd." ne "" %then %do;
			%if "&age_adj" = "" %then %let age_adj = 1; 
			%let stand = %upcase(%scan(&stnd.,1,'('));
			%if %sysfunc(findw(ICSS CDN,&stand.)) = 0 
				%then %err_mess(you must specify either ICSS or CDN standardisation);
			%else %do;
				%let stand_val = %scan(&stnd.,2,'()');
					%if &stand. = ICSS and not %sysfunc(exist(&weight_lib..ICSS)) %then %do;
						%err_mess(ICSS weights file not found ); 
						%goto fin;
					%end;

				%if &stand. = ICSS and %sysfunc(findw(1 2 3 4,&stand_val.)) = 0 
					%then %err_mess(Value must be either 1, 2, 3 or 4 (prostate only)  );
				%else %if &stand. = CDN %then %do; 			/*	build list of allowed CCS codes	*/
					%if not %sysfunc(exist(&weight_lib..CDN)) %then %do;
						%err_mess( Canadian weights file not found ); 
						%goto fin;
					%end;
				proc sql noprint;
					select distinct ccr_fine format 4. into :class
					separated by ' '
					from &weight_lib..CDN;
				quit;
				%if %sysfunc(findw(&class.,&stand_val.)) = 0
					%then %err_mess(Value must be in the list of allowed CCS_fine codes for CDN standardisation) ;
			%end;
		%end;
		
		%if "&strat." ne "" and %sysfunc(findw(&strat., &age.)) ne 0 
			%then %err_mess(Age cannot be an explicit stratifier when age standardising);
	%end;
	%else %if "&weight_var." ^= "" and "&age_adj." = "" %then %let age_adj = 1; 
		
	%end;
	
*	check for incompatable combinations of parameters;
	%let list = %lowcase(&list.);
	%if %sysfunc(findw(ederer pohar,&list.)) = 0 and  "&stnd." ^= "" 
		%then %err_mess(standardisation only available for Ederer and Pohar options);

*	if any errors were detected, then skip all processing and go to end of macro;
	%if &err. %then %goto fin;
*	end of error checking;

*	build age and age format macro strings. Age must be kept separately from
	any other stratification variables, if this is age standardisation run;

*	note specification of age recode, based on standard and value selected
	for most sites, a standard age grouping is used
	the ICSS and CDN weights use a different grouping for prostate
	the CDN weights also define more age groups for breast (ICSS does not);
	%let stand_fmt = ;
	%let age_recode = ;
	%let stand_strat = ;
	%if "&stnd." ne "" %then %do;
		%let group_count = 5;
		%let stand_strat = &age.;
		%let age_recode = &age. = input(put(&age.,z4.), agegrp_o.);
		%let stand_fmt =  &age. agegrpn.; 

		%if &stand. eq ICSS %then %do;
			%if &stand_val = 4 %then %let age_recode = &age. = input(put(&age.,z4.), agegrp_prostate.);
			%let weight_file = icss;
			%let ref_var = standard;
		%end;

		%else %if &stand. eq CDN %then %do;
			%if &stand_val = 700 %then %do;
				%let age_recode = &age. = input(put(&age.,z4.), agegrp_o.);
/*				%let group_count = 6;*/
			%end;
			%else %if &stand_val = 901 %then %let age_recode = &age. = input(put(&age.,z4.), agegrp_prostate.);
			%let weight_file = CDN;
			%let ref_var = ccr_fine;
		%end;
	%end;
	%else %if "&standstrata" ^= "" %then %let stand_strat = &standstrata.;
	
*	variables for analytic file;
	%let vars = d &cens_var. &sex. &age. &stand_strat. &origin. &exit. &strat. &weight_var. &standstrata.;

*	definition of origin, entry and exit macro strings
	also join for merge to life table estimates;
	
*	using dates?;
	%if &yydx. ^= %then %do;
		%let dates = 0;
		%let date_code =;
		%let orig = zero;
		%let origin = zero;
		%let lt_join = a.&yydx.;
		%let vars = &vars. &yydx. zero;
	%end;	
	%else %do;
		%let dates = 1;
		%let orig = &origin.;
		%let date_code = %str(_dxdate = &origin.; format _dxdate date9.;);
		%let vars = &vars _dxdate;
		%let lt_join = year(a._dxdate);
	%end; 
	
	%let ext = &exit.;
	%if "&entry." eq "" %then %let ent = &origin.;
	%else %do;
		%let ent = &entry.;
		%let vars = &vars &ent.;
	%end;

*	if extra mergeby variables in life table, build code for proc sql, add variables to
	list to keep on case file;
	%let sql_merge = ;
	%let sql_vars = ;
	%if &wc. gt 3 %then %do c = 4 %to &wc.;
		%let sql_vars = &sql_vars. %scan(&mergeby.,&c.);
		%let sql_merge = &sql_merge. and a.%scan(&mergeby.,&c.) = b.%scan(&mergeby.,&c.);
	%end;
	%let vars = &vars. &sql_vars.;

*	now build format string used in summation of life table quantities
	and printing results;
	%let format = ;							
	%if "&strat_fmt." ^= "" or "&stnd." ^= "" %then  %let format = format &strat_fmt. &stand_fmt. ;

*	add patient ID to kept variables.  Needed for Pohar-Perme weights;
	%let vars =  &vars. &PatientID.;

*	choice of variables to list in crude estimates;
	%let list_var = interval &strat. &stand_strat. ;

*	define default list elements for EDE POH CUMINC;
	%let scale00 = ;
	%if &list. = ederer %then %do;
		%let list_var = &list_var. l d w cp lo_cp hi_cp cp_star cr se_cr lo_cr hi_cr se_flag;
		%let as_cr = cr;
		%let as_cp = cp;
		%let as_cp_se = se_cp;
		%let as_cr_se = se_cr;
	%end;

	%else %if &list. = pohar %then %do;
		%let list_var = &list_var. l d w cp_p cp_star_w cr_p se_ppe lo_cr_p hi_cr_p se_flag_ppe;
		%let as_cr = cr_p;
		%let as_cp = cp_p;
		%let as_cp_se = se_cp;		*	bit of a crock, this.  should have something specific to PPE;
									*	but we dont actually report confidence limits on age-standardised OBSERVED;
		%let as_cr_se = se_ppe;
	%end;

	%else %if &list. = cuminc %then %do;
		%let list_var = &list_var. l d w cp cr cgc lo_cgc hi_cgc cgo lo_cgo hi_cgo;
	%end;

	%else %let list_var = &list_var. &list.;

*	we need an indicator for the first interval record for a stratum.
	if there is no stratification, then this is just the first record.
	if there is stratification, then we need to select the variable that 
	immediately preceeds the interval number, that is, the second-last variable in the following list;
*	need a variation for the same sort of test in the age standardisation routine;
	
	%let int1 = _n_ = 1;		*	for the processing of crude estimates;
	%let int2 = _n_ = 1;		*	for the processing of standardised estimates;
	%if "&stand_strat." ne "" or "&strat." ne "" %then %let int1 = first.%scan(&strat. &stand_strat. fu, -2) ;
	%if "&strat." ne "" %then  %let int2 = first.%scan(&strat. fu, -2) ;
	

*	order by string for age adjusted report;
	%if "&strat." eq "" %then %let ord_by = fu;
	%else %let ord_by = %sysfunc(translate(&strat.,","," ")) , fu;			*	proc sql ^order by^ statement requires commas;

*	---------------------------------------------------------------;
*	end of processing of parameters and values passed in macro call;
*	---------------------------------------------------------------;

*	now we can get down to business.  build standard case file for analysis;
data _cases_;
	set &infile. ; 
/*	set &infile. end = last; */

*	standard death indicator;
	d = &cens_var. not in( &cens_val.);			*	event indicator (1 = event happened);

*	may need a zero for duration-based analysis;
	zero = .;
	if not &dates. then zero = 0;	
	
*	otherwise, we are in a date-based analysis will need a copy of diagnosis date;
	&date_code.;


*	variables needed for analysis;
	keep &vars. ; 

run;

*	check for duplicated subjects;
%let dups = 0;
proc sql noprint;
		select count(*) as count into : dups from _cases_ 
		group by &patientID. 
		having count >1;
quit;
%if &dups. gt 0 %then %do;
	%err_mess(some subjects appear more than once in input file ); 
	%goto fin;
%end;

*	Split the data to obtain one observation for each life table interval for each individual.
	The scale must be transformed to years;
	%lexis (
		closed = left,
		data = _cases_ , 
		out = _cases_,
		breaks = &intervals.,      
		origin = &orig.,
		entry = &ent.,
		exit = &ext.,
		fail = d,
		scale = &scale.,
		right = right,
		risk = y,				/*	person-years survived in interval	*/
		lrisk = ln_y,
		lint = length,			/*	length of interval					*/
		cint = w,
		os_right = os_right,
		os_left = os_left,
		nint = fu);

*	append survival probability from life tables, based on attained age 
	and attained period for each interval.  further merge by other
	life table variables if requested;
proc sql noprint;
	create table _cases_new as select 
	a.*,
	b.&survprob.
	from _cases_ a left join &popmort. b
		on min(&maxage.,int(a.&age. + a.left)) = b.&merge1.
		and a.&sex = b.&merge2.
		and int(&lt_join. + a.left) = b.&merge3.
		&sql_merge.

	where fu > 0 
/* 		and _st = 1 */

	order by &PatientID., fu ;
	
quit;

*	compute weighted estimates for PPE ;
proc sort data = _cases_new (keep = &PatientID. fu &survprob. length _st)
	nodupkeys out = _Case_weights;
	by &PatientID. fu ;
run;

data _Case_weights;
	set _Case_weights end = last;
	by &PatientID. fu;

	if _n_ = 1 then do;
		n_fu = 1;						*	need the max number of intervals for a later step;
		n_missing = 0;
	end;	

	retain n_fu  WeighEnd n_missing;

	n_fu = max(n_fu, fu);

	if first.&PatientID. then WeighEnd = 1;
	
	if &survprob. ^=. then do;
		p_star 		= &survprob.**length;
		WeighPohar 	= 1/(WeighEnd*sqrt(p_star));	*	use survival prob of mid-point of interval for current interval;
		WeighEnd 	= WeighEnd*p_star;				*	cumulative prob. of survival at end of interval;
	end;
	else if _st = 1 then n_missing =  n_missing+1;

	keep &PatientID. fu WeighPohar n_missing _st;
	if last then do;
		call symputx('nmiss',n_missing);	*	cases did not match to popmort;
		call symputx('Max_fu', n_fu);		*	save max number of intervals in local macro string;	
	end;
run;

*	check for unmatched cases;
	%if &nmiss. gt 0 %then %do;
		%err_mess(&nmiss. cases did not match to population life table &popmort.  ); 
		%goto fin;
	%end;

*	compute life table components at the individual level
	p-y prior to a period window are excluded at this point;
data _cases_;
	merge _cases_new (in=a where = (_st and Y > 1e-10))
	  _Case_weights;
	by &PatientID. fu ;

	if a;
	
	interval = put(left,4.2) || ' - ' || left(put(right,4.2));

	if &survprob. ^= . then do;
  		p_star = &survprob.**length;				* Adjust for interval lengths other than 1 year;
  		d_star = -log(&survprob.)*y;				* expected number of deaths;
	end;

*	if this is a standardisation run based on one of the internal standards, 
	then re-code age at this point;
	%if "&stnd."  ne "" %then &age_recode.;;

*	for PPE, y is adjusted to 1/2 potential interval length if death or censoring in interval;
*	if late entry, then use enter and exit times to get the observable interval length.  
	this will only happen for the first valid interval for a patient ;
	y_p = length;
	if  &ent. ^= os_left then y_p = (os_right-&ent.)/&scale.;
	y_p = ifn(d or w, 0.5*y_p, y_p);  

*	weighted versions of deaths, person-time at risk, contribution to expected hazard of death;
	d_weigh 	= d*WeighPohar;
	d_s_weigh 	= -1*log(&survprob.)*y_p*WeighPohar;
	y_weigh 	= y_p*WeighPohar;
	d2_w		= d*(WeighPohar**2);

	label 
    	d_star  	= 'Expected number of deaths'
    	d       	= 'observed deaths'
    	w       	= 'censored'
    	y       	= 'Person-time (years)*at risk'
    	length  	= 'Interval length*(potential not actual)'
    	ln_y   	 	= 'ln(person-time at risk)'
    	p_star  	= 'Expected survival*probability'
    	interval	= 'Life table*interval'
    	fu      	= 'Follow-up*interval'
		d_weigh 	= 'weighted number*of deaths'
		y_p			= 'Person-time (years)*at risk during*the interval (PPE)'
		d_s_weigh 	= 'weighted expected*deaths'
		y_weigh		= 'weighted number*at risk';
run;

*	Collapse the data to produce the life table estimates;
*	stratification by age is managed by specific macro strings for the case of standardisation
	otherwise, the macro strings are null; 
proc summary data = _cases_ nway;
  	var d w p_star y d_star ;
  	id interval length;
  	class &strat. &stand_strat. &weight_var. fu left right;       
  	output out = &crude_estimates. 
		(drop = _type_ rename = (_freq_ = l)) 
    	sum(d w y d_star d_weigh d_s_weigh  y_weigh d2_w) 	
							= d w y d_star d_weigh d_s_weigh  y_weigh d2_w
    	mean(p_star) 		= p_star;
  		&format.;
run;

*	lifetable quantities;

%let RS_mess =;
data &crude_estimates.; 
  	set &crude_estimates. (where = (fu ^=.));
	by &strat. &stand_strat. fu ;

*	arrays for covariance estimates needed for variance of cummulative crude prob death
	one for each interval. all need to be retained until a new stratum is reached;
	array Pi_p (*) p1-p&Max_fu.;
	array se_t (*) se_t1-se_t&Max_fu.;
	array ex   (*) es_t1-es_t&Max_fu.;
	array px   (*) p_t1-p_t&Max_fu.;
	array nx   (*) n_t1-n_t&Max_fu.;

	length se_flag se_flag_ppe $20 ;

  	retain cp cp_star cr 1 cgc cgo cr_p cp_p cp_star_w v_l_w se_temp
		p1-p&Max_fu. se_t1-se_t&Max_fu. es_t1-es_t&Max_fu. 
		p_t1-p_t&Max_fu. n_t1-n_t&Max_fu. c_v_gxc c_v_gxo fu1;

*	start of stratum (there will not always be an interval 1) - initialise everything in sight;
  	if &int1. then do;
		fu1			= 0;		*	interval order ;
    	cp 			= 1;		*	Cumulative observed survival;
    	cp_star 	= 1;		*	Cumulative expected survival;
    	cr 			= 1;		*	Cumulative relative survival;
    	se_temp 	= 0;
		cgc 		= 0;		*	cumulative crude prob death due to cancer;
		cgo 		= 0;		*	cumulative crude prob death due to other causes;
		cp_p 		= 1;		*	cumulative (weighted) observed survival;
		cp_star_w 	= 1;		*	cumulative (weighted) expected survival;
		cr_p 		= 1;		*	cumulative PPE;
		v_l_w 		= 0;		*	cumulative variance component for PPE variance;
		c_v_gxc		= 0;		*	variance of crude prob death due to cancer;
		c_v_gxo		= 0;		*	variance of crude prob death due to other;
		do _i_ = 1 to &Max_fu.;
			Pi_p(_i_) = 0;
			se_t(_i_) = 0;
			ex(_i_)	  = 0;
			px(_i_)   = 0;
			nx(_i_)   = 0;
		end;
  	end;

	p 			= exp(-(d/y)*length);	*	hazard definition of survival probability;
	var_lambda	= length**2*d/y**2;		*	variance of hazard;
  	l_prime 	= l - w/2;				*	Effective number at risk;
	ns 			= l_prime - d;			*	Number surviving the interval;

*	other interval-specific calculations;
*	test for conditional (interval-specific) RS that is greater than 1
	set a flag, so that a message can be written later;
  	r 				= p/p_star;						*	Interval-specific relative survival;
  	cp 				= cp*p;							*	Cumulative observed survival;
  	cp_star 		= cp_star*p_star;				*	Cumulative expected survival;
 	if p > p_star then call symput("rs_mess","y");

	cr 				= cp/cp_star;					*	Cumulative relative survival;
 	ln_y 			= log(y);
  	d_star_group 	= l_prime*(1 - p_star);
  	diff 			= d-d_star;
  	excess 			= diff/y;
  	
 /*	not certain that this makes any sense.  Excess can be negative, but poisson			*/
 /*	limits are bounded by zero. 														*/
 /*		caveat emptor  16 jul 2017														*/
	if excess > 0 then do;
		%cipoiss(0.95, diff, ll,lu);
		excess_uci 		=  lu/y;
		excess_lci 		=  ll/y;
	end;
	
	se_temp 		= se_temp + var_lambda;			*	Component of the variance of the cumulative survival;
	se_p			= p*sqrt(var_lambda);
  	if r ^= 0 then  
		se_r 		= se_p/p_star;	
  	se_cp 			= cp*sqrt(se_temp);
  	se_cr 			= se_cp/cp_star;


*	6 Dec 2017:
	Cumulative probabilities of death.  (CIF as in Cronin and Feuer 2002);
	fu1 = fu1+1;

*	observed and expected survival probabilities for this interval;
	px(fu1)	= p;
	ex(fu1)	= p_star;		

*	variance-equivalence estimate of effective number at risk (Enzo Coviello, 23 Dec, 2014);
	nx(fu1) 	= d/2 + sqrt((d/2)**2 + (y/length)**2);

*	temporary calculations for interval-specific crude prob death (due to cancer/other);
	c_cgc 		= se_temp - var_lambda + ((p/(p_star-p))**2)*var_lambda;
	c_cgo 		= se_temp - var_lambda + ((p/(p_star+p))**2)*var_lambda;

*	calculations for crude probabilities of death;
	gxc			= cp*(1/p)*(1-p/p_star)*(1-.5*(1-p_star));	*	gxc in Cronin & Feuer;
	gxo			= cp*(1/p)*(1-p_star)*(1-.5*(1-p/p_star));	*	gxo in Cronin & Feuer;
	Pi_p(fu1) 	= cp;
	se_t(fu1)	= se_temp;

*	variance calculations for crude probabilities of death;
	se_gxc		= gxc*sqrt(c_cgc);			*	se of interval crude prob of death due to cancer;
	se_gxo		= gxo*sqrt(c_cgo);			*	se of interval crude prob of death due to other causes;

*	from p 1733 of Cronin & Feuer;
	if fu1 = 1 then do;
		c_v_gxc	= se_gxc**2;				*	variance of first interval;
		c_v_gxo	= se_gxo**2;
		var_gxc = c_v_gxc;
		var_gxo = c_v_gxo;
	end;
	else do;
		c_v_gxc	= c_v_gxc + se_gxc**2;		*	cumulative sums of variances;
		c_v_gxo	= c_v_gxo + se_gxo**2;
		cov_c = 0;
		cov_o = 0;
		do _k_ = 2 to fu1;
			do _j_ = 2 to _k_; 
				if _j_ = 2 then do;
					PP_1 = 1;
					ne = 0;
				end;
				else do;
					PP_1 = Pi_p(_j_-2);
					ne = se_t(_j_-2);
				end;
				cov_c = cov_c + (
					Pi_p(_k_-1)*PP_1
					*(1-.5*(1-ex(_j_-1)))*(1-.5*(1-ex(_k_)))
					*(1-px(_k_)/ex(_k_))*(1-px(_j_-1)/ex(_j_-1))
					*(ne-(1-px(_j_-1))/((ex(_j_-1)-px(_j_-1))*nx(_j_-1)))); 

				cov_o = cov_o + (
					Pi_p(_k_-1)*PP_1
					*(1-ex(_j_-1))*(1-ex(_k_))
					*(1-.5*(1-px(_k_)/ex(_k_)))*(1-.5*(1-px(_j_-1)/ex(_j_-1)))
					*(ne+(1-px(_j_-1))/((ex(_j_-1)+px(_j_-1))*nx(_j_-1)))); 
			end;
		end;

		var_gxc = c_v_gxc + 2*cov_c;		*	variances + 2*covariances;
		var_gxo = c_v_gxo + 2*cov_o;
	end;				 

*	confidence limits on interval-specific crude probabilities;
	lo_gxc		= gxc-1.96*se_gxc;
	hi_gxc		= gxc+1.96*se_gxc;
	lo_gxo		= gxo-1.96*se_gxo;
	hi_gxo		= gxo+1.96*se_gxo;

*	cumulative crude probabilities and confidence limits;
	cgc 		= cgc + gxc;				*	Gxc in Cronin and Feuer;
	cgo 		= cgo + gxo;				*	Gxo in Cronin and Feuer;
		
*	If either cgc or cgo fall outside of the natural domain (0,1), then the 
	confidence limits are undefined;
	if 1 > cgc > 0 and var_gxc > 0 then do;
		lo_cgc		= cgc**exp(-1.96*sqrt(var_gxc)/(cgc*log(cgc)));
		hi_cgc		= cgc**exp(1.96*sqrt(var_gxc)/(cgc*log(cgc)));
	end;
	if 1 > cgo > 0 and var_gxo > 0 then do;
		lo_cgo		= cgo**exp(-1.96*sqrt(var_gxo)/(cgo*log(cgo)));
		hi_cgo		= cgo**exp(1.96*sqrt(var_gxo)/(cgo*log(cgo)));
	end;
*	end of CIF calculations;

*	PPE calculations  (PPE weighted y can be 0 in situations of late entry)   19 June 2020;
	if y_weigh > 0 then do;
		os_w 			= exp(-length*(d_weigh/y_weigh));			*	interval observed survival (weighted);
		es_w 			= exp(-length*(d_s_weigh/y_weigh));			*	interval expected survival (weighted);
		ns_w 			= exp(-length*(d_weigh-d_s_weigh)/y_weigh);	*	interval net survival;
		cp_p 			= cp_p*os_w;								*	cumulative observed survival (weighted);
		cp_star_w	 	= cp_star_w*es_w;
		cr_p 			= cr_p*ns_w;
		v_os_w 			= (os_w**2)*d2_w/(y_weigh**2);		*	variance of interval observed survival (weighted);
		v_ns_w			= (ns_w**2)*d2_w/(y_weigh**2);		*	variance of interval PPE;
		v_l_w			= v_l_w + (length**2)*d2_w/(y_weigh**2);	*	variance of cumulative relative surival (PPE);	
		se_ppe			= cr_p*sqrt(v_l_w);					*	se of cumulative PPE relative survival;

*	calculations for interval net survival CI;
		se_pw 		= os_w*length*sqrt(d_weigh)/y_weigh;
		if "&ci." = "loghaz" then do;
			if 0 < os_w and os_w ^= 1 then do;
				se_lh_pw	= -length*sqrt(d_weigh)/(y_weigh*log(os_w));
				lo_lh_pw	= log(-log(os_w))+1.96*se_lh_pw;
				hi_lh_pw	= log(-log(os_w))-1.96*se_lh_pw;
				lo_pw		= exp(-exp(lo_lh_pw));
				hi_pw		= exp(-exp(hi_lh_pw));
				lo_rw		= lo_pw/es_w;
				hi_rw		= hi_pw/es_w;
			end;
			else do;
				lo_rw		= .;
				hi_rw		= .;
			end;
		end;
		if "&ci." = "haz" then do;
			lo_pw		= os_w*exp(-se_pw*1.96/os_w);
			hi_pw		= os_w*exp(se_pw*1.96/os_w);
			lo_rw		= lo_pw/es_w;
			hi_rw		= hi_pw/es_w;
		end;
		
	end;

*	Pohar-Perme cumulative net survival;
	if "&ci." = "loghaz" then do;			*	loghaz transform;
*	if cr_p >= 1 then confidence limits are undefined;
*	if cr_p is very small, or very close to 1, then the computation below involves division by zero;
*	use the ci = log specification, or request Ederer II estimates of relative survival;
		if 0 < cr_p < 1  then do;
			if abs( log(-log(cr_p)) - se_ppe*1.96/(cr_p*log(cr_p)) ) < 250 then do;
				lo_cr_p 		=  exp(-exp(log(-log(cr_p)) - se_ppe*1.96/(cr_p*log(cr_p)))); 
				hi_cr_p 		=  exp(-exp(log(-log(cr_p)) + se_ppe*1.96/(cr_p*log(cr_p))));   
				se_lh_cpp 		=  sqrt( se_ppe**2/(cr_p*log(cr_p))**2 );
			end;
		end;
		else do;	*	confidence limits undefined;
			lo_cr_p = .;
			hi_cr_p = .;
		end;
	end;
	else if "&ci." = "haz" then do;		*	log transform;
		lo_cr_p	= cr_p*exp(-se_ppe*1.96/cr_p);
		hi_cr_p	= cr_p*exp(se_ppe*1.96/cr_p);
	end;

*	interval Ederer II estimates;
	if "&ci." = "loghaz" then do;			*	loghaz transform;
*	Ederer II estimates;	
	  	if p > 0 and  p ^= 1 then do;			*	First for the interval specific estimates;
	    	se_lh_p = sqrt(se_p**2/(p*log(p))**2);			*	SE on the log-hazard scale using Taylor series approximation;

*	Confidence limits on the log-hazard scale;
	    	lo_lh_p = log(-log(p)) + 1.96*se_lh_p;
	    	hi_lh_p = log(-log(p)) - 1.96*se_lh_p;

*	Confidence limits on the survival scale (observed survival);
	    	lo_p = exp(-exp(lo_lh_p));
	    	hi_p = exp(-exp(hi_lh_p));

*	Confidence limits for the corresponding relative survival rate;
	    	lo_r = lo_p/p_star;
	    	hi_r = hi_p/p_star;
		end;
	end;

	if "&ci." = "haz" then do;			*	log transform;
	    	lo_p = p*exp(-se_p*1.96/p);
	    	hi_p = p*exp(se_p*1.96/p);

	    	lo_r = lo_p/p_star;
	    	hi_r = hi_p/p_star;
	end;

*	the cumulative estimates;
	if "&ci." = "loghaz" then do;
	    if  1.0e-150 < cp < 1 then do;
	    	se_lh_cp = sqrt( se_cp**2/(cp*log(cp))**2);        * SE on the log-hazard scale using Taylor series approximation;

*	Confidence limits on the log-hazard scale;
	    	lo_lh_cp = log(-log(cp)) + 1.96*se_lh_cp;
	    	hi_lh_cp = log(-log(cp)) - 1.96*se_lh_cp;

*	Confidence limits on the survival scale (observed survival);
	    	lo_cp = exp(-exp(lo_lh_cp));
	    	hi_cp = exp(-exp(hi_lh_cp));

*	Confidence limits for the corresponding relative survival rate;
	    	lo_cr = lo_cp/cp_star;
	    	hi_cr = hi_cp/cp_star;
	    end;
		else do;
	        	lo_cp = .;
	    		hi_cp = .;
		    	lo_cr = .;
	    		hi_cr = .;
		end;
	end;

*	using log transform;
	if "&ci." = "haz" then do;
	    lo_cp = cp*exp(-se_cp*1.96/cp);
	   	hi_cp = cp*exp(se_cp*1.96/cp);

	   	lo_cr = lo_cp/cp_star;
	   	hi_cr = hi_cp/cp_star;
	end;
*	end of interval EII estimates;

*	evaluation of variability of survival estimate(s);
		se_flag = 'SE unavailable';
		if se_cr ^= . then do;
	 		if se_cr > 0.1 then se_flag = 'Supress';
	   		else if se_cr > 0.05 then se_flag = 'Potentially unstable';
	        else  se_flag = 'No comment';
		end;
		
		se_flag_ppe = 'SE unavailable';
		if se_ppe ^= . then do;
	 		if se_ppe > 0.1 then se_flag_ppe = 'Supress';
	   		else if se_ppe > 0.05 then se_flag_ppe = 'Potentially unstable';
	        else  se_flag_ppe = 'No comment';
		end;
		
		

 *	Formats and labels;
    format lo_p hi_p lo_r hi_r lo_cp hi_cp lo_cr hi_cr lo_cr_p hi_cr_p lo_rw hi_rw 
		cgc cgo lo_cgc hi_cgc lo_cgo hi_cgo 8.5;
    label 
		lo_p 		= 'Lower 95%*CI P'
		hi_p 		= 'Upper 95%*CI P'
		lo_r 		= 'Lower 95%*CI R'
		hi_r 		= 'Upper 95%*CI R'
		lo_cp 		= 'Lower 95%*CI CP'
		hi_cp 		= 'Upper 95%*CI CP'
		lo_cr 		= 'Lower 95%*CI CR'
		hi_cr 		= 'Upper 95%*CI CR'
		interval		= 'Interval'
		fu			= 'Interval'
		l			= 'Alive at start'
		l_prime		= 'Effective*number at risk'
		ns			= 'Number surviving*the interval'
		d			= 'Deaths'
		w			= 'Withdrawals'
		p			= 'Interval-specific*observed survival'
		cp			= 'Cumulative*observed survival'
		r			= 'Interval-specific*relative survival'
		cr			= 'Cumulative*relative survival'
		p_star		= 'Interval-specific*expected survival'
		cp_star		= 'Cumulative*expected survival'
	    ln_y         = 'ln(person-time)'
	    d_star_group = 'Expected deaths (approximate)'
	    excess       = 'Empirical excess hazard'
		excess_uci 	 = 'Upper 95% CI*excess hazard'
		excess_lci 	 = 'Lower 95% CI*excess hazard'
		se_p         = 'Standard error of P'
	    se_r         = 'Standard error of R'
	    se_cp        = 'Standard error of CP'
	    se_cr        = 'Standard error of CR'
	    se_flag		 = 'Suggested*Interpretation'
	    se_flag_ppe	 = 'Suggested*Interpretation'
		cgc			 = 'Cumulative incidence*death due*to cancer'
		cgo			 = 'Cumulative incidence*death due*to other causes'
		lo_cgc		 = 'Lower 95%*CI cgc'
		hi_cgc		 = 'Upper 95%*CI cgc'
		lo_cgo		 = 'Lower 95%*CI cgo'
		hi_cgo		 = 'Upper 95%*CI cgo'
		ns_w		 = 'Interval-specific*net survival (PPE)'
		lo_rw		 = 'Lower 95% CI*interval net survival'
		hi_rw		 = 'Upper 95% CI*interval net survival'
		cp_p 		 = 'cumulative (weighted)*observed survival'
		cp_star_w 	 = 'cumulative (weighted)*expected survival'
		cr_p 		 = 'cumulative net*survival (PPE)'
		lo_cr_p 	 = 'Lower 95%*CI PPE'
		hi_cr_p 	 = 'Upper 95%*CI PPE';

  	drop se_temp p1-p&Max_fu. se_t1-se_t&Max_fu. es_t1-es_t&Max_fu. 
		p_t1-p_t&Max_fu. n_t1-n_t&Max_fu. _i_ _k_ _j_
		v_l_w  var_lambda  cov_c cov_o pp_1 ne os_w es_w 
		v_ns_w ;
run;

*	note if interval-specific estimates are greater than 1;
%if "&RS_mess." ne "" %then %do;
	%put Conditional relative/net Survival is greater than 1 for some intervals;
	%if "&ci." eq "loghaz" %then 
		%put correspoding confidence intervals are undefined and have been set to missing;
%end;


*	Print crude rates estimates (if requested);
proc print data = &crude_estimates. noobs label split = '*';
	where &crude.;
  	var &list_var.;
  	&format.;
run;
 
 
*	if no age standardisation is being requested, then we are finished;
%if "&stnd." ne "" or "&weight_var." ne "" %then %do;

*	age standardisation;

*	apply weights to grouped data if not already available; 
%if "&stnd." ^= "" %then %do;
	proc sort data = &weight_lib..&weight_file. ( 
			rename = (agegroup = &age.) 
			where = (&ref_var. = &stand_val.))
		out = weights (keep= weight &age.)
		nodupkeys;
		by &age.;
	run;

	proc sort data = &crude_estimates.;
		by &age.;
	run:
%end;

*	either append weights (implicit merge on age) or
	accept weights already assigned
	also need a count of stratification strata;
%if "&stnd." ^= "" %then %do;
	%let w_merge = 
	%str(merge &crude_estimates. (in = m1) weights; by &age.; if m1;); 
	%end;
%else %do;
	%let w_merge = 
	%str(set &crude_estimates. (rename = (&weight_var. = weight)));
	proc sql noprint;
		select count(*) into :group_count 
		from (select distinct &standstrata. from &crude_estimates.);
	quit;
%end;
	
data cr_weight;
   	&w_merge;

   	as_cr = weight*&as_cr.;   
	as_cp = weight*&as_cp.;   
   	as_se_cp_sqr = (weight*&as_cp_se.)**2;		*	weighted variance (observed);
   	as_se_cr_sqr = (weight*&as_cr_se.)**2;		*	weighted variance (relative);
run;


proc sort data = cr_weight;
	by &strat. &stand_strat. fu ;
run; 

*	if there are no patients alive at the end of an interval, it can happen in two ways:
	The remaining patients all die in the interval.  In this case the survival proportion 
		is zero.  Such an interval can still contribute to the age-standardisation computations
		condition to test for is l = d
	The patients entering the innerval are all censored in the interval.  In this case, the
		survival proportion is 1, and is undefined for subsequent intervals.

	in the first instance, set the survival estimate to zero for subsequent non-valid intervals,
	and copy the variance estimate from the previous interval;

data cr_weight0;
	set cr_weight;
	by &strat. &stand_strat. fu ;

	retain  l_as_se_cr_sqr l_as_se_cp_sqr  prev_missing;
	
	l_as_se_cp_sqr 	= as_se_cp_sqr;
	l_as_se_cr_sqr 	= as_se_cr_sqr;

	
	if first.&stand_strat. then prev_missing = 0;
	
	if l = w then do;
		output;
		prev_missing = 1;
	end;
	else if not prev_missing then output;
	
	if l=d and last.&stand_strat. and not prev_missing then do;
		do i = fu+1 to &Max_fu.;
			fu = i;
			as_cr = 0;
			as_cp = 0;
			as_se_cp_sqr 	= l_as_se_cp_sqr;	*	use variance estimate from previous interval;
			as_se_cr_sqr 	= l_as_se_cr_sqr;	*	use variance estimate from previous interval;
			output;
		end;
	end;

	drop  i l_as_se_cp_sqr l_as_se_cr_sqr;
run;

proc sort data = cr_weight0;
	by &strat. fu &stand_strat. ;
run; 


*	Test to see if there are RS estimates for all age groups, in each combination of
	stratifier, and follow-up interval;
proc means data = cr_weight0 noprint nway;
	class &strat. fu  ;
	output out = how_many
		(drop = _freq_ _type_)
		n(fu) = strat_count;
run;

data cr_weight_1;
	merge cr_weight0 ( in = m1 ) 
	how_many;
	by &strat. fu;

	if m1;
	if &int2. then missing_flag = 0;
	
	retain missing_flag;

*	If at least one age-specific value is missing, then set all to missing for this interval
	and all subsequent intervals;
	if strat_count ^= &group_count. or missing_flag = 1 then do;                
		as_cr = .;          
		as_cp = .;
		as_se_cp_sqr = .;
		as_se_cr_sqr = .;
		missing_flag = 1;
	end;

	drop strat_count missing_flag;
run;

proc summary data = cr_weight_1 nway;
	class &strat. fu ; 
	var as_cr as_cp as_se_cr_sqr;
	output out = rs_ageadj  
		sum( as_cr ) = cr_ageadj 
		sum( as_cp ) = cp_ageadj 
		sum(as_se_cp_sqr) = var_p_ageadj
		sum(as_se_cr_sqr) = var_r_ageadj;
run; 

*	look up interval labels;
proc sql;
	create table rs_ageadj2 as
		select a.*, b.interval
		from rs_ageadj a left join (
			select distinct fu, interval from cr_weight) b
			on a.fu = b.fu
	order by &ord_by.;  
quit;


  *****************************************************************************
  **  Combine age-specific and standardized data. Calculate the age-specific **
  **  components of the variance estimate for age-standardized survival      **
  **  estimates. If the standard  error of an age-specific observed survival **
  **  is undefined (i.e. observed survival is 0 or 1) then set the           ** 
  **  variance estimate component to zero. If the age-standardized observed  **
  **  survival estimate itself is undefined (see directly above - missflag=1)**
  ** then all variance components will be set to zero.                       **
  *****************************************************************************;

  *****************************************************************************
  ** Determine age standardized CIs (undefined if the age standardized point **
  ** estimate is 0, 1, or undefined itself).                       	         **
  *****************************************************************************;
data &std_estimates.;
	length se_flag $20. ;
	set rs_ageadj2 (in = a rename = 
		(cp_ageadj = as_obs
		cr_ageadj = as_rel));

	if var_p_ageadj ^= . then se_cp_ageadj = sqrt(var_p_ageadj);
	if var_r_ageadj ^= . then se_cr_ageadj = sqrt(var_r_ageadj);
	
*	standardised observed survival is missing if >= 1;
	if as_obs < 1 and se_cr_ageadj ^= . then do;
		lo_obs = as_obs**exp( -1.96*abs(se_cp_ageadj/(as_obs*log(as_obs))));   
		hi_obs = as_obs**exp( 1.96*abs(se_cp_ageadj/(as_obs*log(as_obs))));   
	end;

*	using CI settings as of 26 Jan 2021;
	if "&ci." = "loghaz" then do;
		if 0< as_rel and as_rel ^= 1 then do;
			hi_rel = as_rel**exp( -1.96*abs(se_cr_ageadj/(as_rel*log(as_rel))));   
			lo_rel = as_rel**exp( 1.96*abs(se_cr_ageadj/(as_rel*log(as_rel)))); 
		end;
		else do;
			hi_rel = .;
			lo_rel = .;
		end; 
	end;
	if "&ci." = "haz" then do;
		if 0< as_rel and as_rel ^= 1 then do;
			hi_rel = as_rel*exp(1.96*se_cr_ageadj/as_rel) ; 
			lo_rel = as_rel*exp(-1.96*se_cr_ageadj/as_rel) ; 
		end;
		else do;
			hi_rel = .;
			lo_rel = .;
		end; 
	end;

	ci_rel = put( lo_rel, 4.2 ) || ' - ' || put( hi_rel, 4.2 );

	if se_cr_ageadj ^= . then do;
		if se_cr_ageadj > 0.1 then se_flag = 'Supress';
		else if se_cr_ageadj > 0.05 then se_flag = 'Potentially unstable';
		else  se_flag = 'No comment';
	end;
	else do ;
		se_flag = ' ';
	end;

	right = input(scan(interval,-1,' '),6.2);
	left = input(scan(interval, 1,' '),6.2);
	label
	right			= 'right end of interval'
	left			= 'left end of interval'
	as_obs			= 'Age Standardised*Observed Survival'
	as_rel 			= 'Age Standardised*Rel Survival'
	ci_rel 			= '95% CI*Rel Survival'
	se_cp_ageadj	= 'SE Obs Survival'
	se_cr_ageadj	= 'SE Rel Survival'
	se_flag 		= 'suggested*interpretation'
	lo_rel 			= 'lower 95% CI*Rel Surv'
	hi_rel			= 'upper 95% CI*Rel Surv';
		
	keep &strat. interval left right as_obs as_rel ci_rel se_cp_ageadj se_cr_ageadj se_flag lo_rel hi_rel;
run;

proc print data = &std_estimates. noobs label split='*';
	where &age_adj.;
	var &strat. interval as_obs as_rel se_cr_ageadj lo_rel hi_rel se_flag ;
    format as_obs lo_rel hi_rel as_rel 6.4;
run;

%end;

*	clean up temporary files;
proc datasets library = work nolist force nowarn;
	delete _case_weights _cases_new discrd rs_ageadj2 cr_weight0 cr_weight1 
		how_many weights cr_weight;
quit;
run;

*	end of processing (also error return);
%fin:
%mend rel_surv;		*	end of relative survival macro;


/**************************************************************************
Author: Bendix Carstensen, 1999-2002
Update: Paul Dickman, BxC, November 2003
Bug-fix: BxC, December 2007:
         If the origin= argument had missing values erroneous output would
         be generated (too much risk time). Now remedied so that
         observations with missing values of origin are excluded.
This macro is in: http://www.biostat.ku.dk/~bxc/Lexis/Lexis.sas
Example program:  http://www.biostat.ku.dk/~bxc/Lexis/xLexis.sas 
***************************************************************************/

/*	added closed = right (default) to retain original functionality that
	lines up with Stata stsplit
	closed = left will add a small amount of time to any event that 
	falls on a break value.  This will have the effect of putting that
	observation into the correct interval, as if the intervals were
	closed left and open right (even though they are not)
	corresponds to the correction made in Stata strs to handle the
	same situation
	
	reduced output so redundant intervals are not output
	30 June 2017  Ron Dewar
*/


%macro Lexis ( data = ,       /* Data set with original data,             */
                              /*    defaults to _last_                    */
                out = ,       /* Where to put the result,                 */
                              /*    defaults to &data.                    */
              entry = entry,  /* Variable holding the entry date          */
               exit = exit,   /* Variable holding the exit date           */
               fail = fail,   /* Variable holding the exit status         */
                              /* If any of the entry, exit or fail        */
                              /*    variables are missing the person is   */
                              /*    discarded from the computations.      */
             breaks = ,       /* Specification of the cutpoints on        */
                              /*    the transformed scale.                */
                              /*    Syntax as for a do statement.         */
                              /*    The ONLY Mandatory argument.          */
 							  /*  the specification ^failures^ generates  */
 							  /*  breaks at each distinct failure time	  */	
               cens = 0,      /* Code for censoring (may be a variable)   */
              scale = 1,      /* Factor to transform from the scale       */
                              /*    of entry and exit to the scale        */
                              /*    where breaks and risk are given       */
             origin = 0,      /* Origin of the transformed scale          */
               risk = risk,   /* Variable recieving the risk time         */
              lrisk = lrisk,  /* Variable recieving the log(risk time)    */
               left = left,   /* Variable recieving left  endpoint of int */
              other = ,       /* Other dataset statements to be used such */
                              /*     as: %str( format var ddmmyy10. ; )   */
                              /*     or: %str( label risk ="P-years" ; )  */
               disc = discrd, /* Dataset holding discarded observations   */
           /*-------------------------------------------------------------*/
           /* Variables for making life-tables and other housekeeping:    */
           /* These will only appear in the output dataset if given here  */
           /* The existence of these arguments are tested in the macro so */
           /* they cannot have names that are also logical operators such */
           /* as: or, and, eq, ne, le, lt, gt.                            */
           /*-------------------------------------------------------------*/
              right = ,       /* Variable recieving right endpoint of int */
               lint = ,       /* Variable recieving interval length       */
            os_left = ,       /* Variable recieving left  endpoint of int */
           os_right = ,       /* Variable recieving right endpoint of int */
            os_lint = ,       /* Variable recieving interval length       */
                              /*    - the latter three on original scale  */
               cint = ,       /* Variable recieving censoring indicator   */
                              /*    for the current input record          */
               nint = ,       /* Variable recieving index of follow-up    */
                              /*       interval                           */
                 st = _st,    /* flag to indicate if subject is at risk   */
                			  /* in this interval						  */
             closed = right   /* intervals are closed on right by default */
            				  /* life table intervals used by %rel_surv   */
            				  /* are open on right (closed left)          */
              );

%if "&breaks." = "" %then %put ERROR: breaks MUST be specified. ; 
%if %upcase("&breaks.") = "FAILURES" %then %do;
%put Generating breaks for each distinct failure time;
proc sql noprint;
	select distinct ( &exit.  -  &origin. ) / &scale. 
		format 9.2 into : breaks
		separated by ', '
		from &data.
		where &fail.
		order by ( &exit.  -  &origin. ) / &scale.;
quit;
%end;

*	use right closed intervals (default) or left closed (life tables);
%if %substr(%lowcase(&closed.),1,1) = r %then %let r_close = 1;
%else %let r_close = 0;


%if &data.  = %then %let data = &syslast. ;
%if &out. 	= %then %do ;
                    %let out=&data. ;
                    %put
NOTE: Output dataset not specified, input dataset %upcase(&data.) will be overwritten. ;
                  %end ;

data &disc. &out. ;
  set &data. ;
  if ( nmiss ( &entry., &exit., &fail., &origin. ) gt 0 ) then output &disc. ;
  else do;
/*                goto next ; */
/*           end ; */
  * Labelling of variables ;
  label &entry.  = 'Entry into interval' ;
  label &exit.   = 'Exit from interval' ;
  label &fail.   = 'Failure indicator for interval' ;
  label &risk.   = 'Risktime in interval' ;
  label &lrisk.  = 'Natural log of risktime in interval' ;
  label &left.   = 'Left endpoint of interval (transformed scale)' ;
%if    &right.^= %then  label &right. = 'Right endpoint of interval (transformed scale)' ; ;
%if     &lint.^= %then  label &lint. = 'Interval width (transformed scale)' ; ;
%if  &os_left.^= %then  label &os_left. = 'Left endpoint of interval (original scale)' ; ;
%if &os_right.^= %then  label &os_right. = 'Right endpoint of interval (original scale)' ; ; 
%if  &os_lint.^= %then  label &os_lint. = 'Interval width (original scale)' ; ;
%if     &cint.^= %then  label &cint. = 'Indicator for censoring during the interval' ; ;
%if     &nint.^= %then  label &nint. = 'Sequential index for follow-up interval' ; ;
  &other. ;
  drop _entry_ _exit_ _fail_
       _origin_ _break_
       _cur_r _cur_l _int_r _int_l
       _first_ _cint_ _nint_;

/*
Temporary variables in this macro:

  _entry_  holds entry date on the transformed timescale
  _exit_   holds exit  date on the transformed timescale
  _fail_   holds exit  status
  _break_  current cut-point
  _origin_ origin of the time scale
  _cur_l   left  endpoint of current risk interval
  _cur_r   right endpoint of current risk interval
  _int_l   left  endpoint of current break interval
  _int_r   right endpoint of current break interval
  _first_  indicator for processing of the first break interval
  _cint_   indicator for censoring during the interval
  _nint_   sequential index of interval
   
If a variable with any of these names appear in the input dataset it will
not be present in the output dataset.
*/

  _origin_ = &origin. ;
  _entry_  = ( &entry. - _origin_ ) / &scale. ;
  _exit_   = ( &exit.  - _origin_ ) / &scale. ;
  _fail_   = &fail. ;
  _cur_l   = _entry_ ;
  _first_  = 1 ;

  do _break_ = &breaks. ;
     if _first_ then do ;
        _nint_=-1;
        _cur_l = max ( _break_, _entry_ ) ;
        _int_l = _break_ ;
     end ;
     _nint_ + 1;
     _first_ = 0 ;
     _int_r = _break_ ;
     _cur_r = min ( _exit_, _break_ ) ;
/*     if _cur_r gt _cur_l then do ;*/
	 &st. = 0;
     if _cur_r gt _cur_l then &st. = 1 ;
/*
Endpoints of risk interval are put back on original scale.
If any of left or right are specified the corresponding endpoint
of the break-interval are output.
*/
        &entry.  = _cur_l * &scale. + _origin_ ;
        &exit.   = _cur_r * &scale. + _origin_ ;
        &risk.   = _cur_r - _cur_l ;
		if &risk. > 0 then &lrisk.  = log ( &risk. ) ;

/*	if not &r_close. and _exit_ = _break_ then _exit_ = _exit_ + 0.000001;*/
        &fail.   = _fail_ * ( _exit_ eq _cur_r ) +
                   &cens. * ( _exit_ gt _cur_r ) ;
        _cint_   = not( _fail_ ) * ( _exit_ eq _cur_r ) ;            

        %if     &left.^= %then &left.     = _int_l ; ;
        %if    &right.^= %then &right.    = _int_r ; ;
        %if     &lint.^= %then &lint.     = _int_r - _int_l ; ;
        %if  &os_left.^= %then &os_left.  = _int_l * &scale. + _origin_ ; ;
        %if &os_right.^= %then &os_right. = _int_r * &scale. + _origin_ ; ;
        %if  &os_lint.^= %then &os_lint.  = ( _int_r - _int_l ) * &scale. ;;
        %if     &cint.^= %then &cint.     = _cint_ ; ;
        %if     &nint.^= %then &nint.     = _nint_ ; ;

        if not( &st. = 0 and _int_l > _exit_ ) then output &out. ;
/*     end ;*/
     _cur_l = max ( _entry_, _break_ ) ;
     _int_l = _break_ ;
  end ;
/*   next: ; */
end;
run ;

%mend Lexis ;


%macro err_mess(message);
	%put ;
	%put ERROR: *** Error in parameter specification ***;
	%put ERROR: *** &message.;
	%put ERROR: ***;
	%let err = 1;
%mend err_mess;

/*
************************************************************************
*  CIPOISS: A SAS macro for exact Poisson confidence limits            *
*                                                                      *
*  This SAS macro calculates the exact confidence limits for           *
*  a Poisson count.                                                    *
*                                                                      *
*  Algorithm:        Daly, Leslie,  "Simple SAS macros for the         *
*                    calculation of exact binomial and Poisson         *
*                    Confidence limits"                                *
*                    Comput Biol Med, 22(5):351-361,  1992.            *
*                                                                      *
*     Leslie Daly,
*     Lecturer in Medical Statistics,
*     Department of Public Health Medicine and Epidemioilogy,
*     University College Dublin,
*     Earlsfort Terrace,  Dublin 2, IRELAND
*
*     Email:  LDALY@IVEAGH.UCD.IE
*
*  Functions called: LOG                                               *
*                    GAMINV                                            *
*                                                                      *
*  Usage:            %CIPOISS(cl,x,ll,lu)                              *
*                                                                      *
*  Input parameters: cl - confidence level as a proportion             *
*                         (e.g.  0.95)                                 *
*                    x  - observed number of events                    *
*                                                                      *
*  Output parameters: ll - Lower confidence limit                      *
*                     lu - Upper confidence limit                      *
*                                                                      *
*  Usage notes:      Missing values for the output parameters are      *
*                    generated if (1) X is less than zero, or if       *
*                    (2) CL is not between 0.0 and 1.0.                *
* 																	   *
* 	note: small values of x <= .006 cause numerical problems		   *
* 	R Dewar  18 july 2017											   *
*                                                                      *
************************************************************************
*/
%macro cipoiss(cl,x,ll,lu);
if (&x) lt 0 or not(0 lt (&cl) lt 1) then
   do;
      &ll = .;
      &lu = .;
   end;
else
   do;
      if (&x) gt 0.006 then
         do;
            &ll = gaminv((1 - (&cl))/2, (&x));
            &lu = gaminv((1 + (&cl))/2, (&x) + 1);
         end;
      else if 0 <= (&x) <= 0.006 then
         do;
            &lu = -log(1 - (&cl));
            &ll = 0;
         end;
   end;
%mend;
