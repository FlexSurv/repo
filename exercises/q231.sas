/*
==================
 EXERCISE 231
 
 Modelling excess mortality using flexible parametric models I
 
 sql sgplot sgpanel loess tabulate
 
 %stset %rcsgen %stpm2 %predict
 
reviewed:  5 May 2020
=====================

*/
/* Start timer */
%let _timer_start = %sysfunc(datetime());

options fmtsearch = (data.colon_formats);

*	(a) Load Data, merge expected mortality and fit inital model ;
data colon;
	set data.colon;	
	id = _n_;
	
	if surv_mm > 60.5 then do;
		surv_mm = 60.5;
		status = 0;
	end;
	
	female = sex = 2;

	surv = surv_mm/12;
	agegrp2 = agegrp = 1;
	agegrp3 = agegrp = 2;
	agegrp4 = agegrp = 3;
run;

*	append mortality rates and sort by patient ID;
proc sql;
	create table colon_mm as select a.*,
		b.rate
		from colon a left join data.popmort b
		on floor(a.age + a.surv) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv) = b._year
	order by id;
quit;

*	use all causes of death as an event.  Other status codes are censored events;
%stset(colon_mm, status(1 2), surv, id);

*	fit base model with no covariates;
%stpm2( , df=5, scale=hazard, bhazard=rate, options = eform );


*	(b) evaluate residuals wrt age;
%predict(mg1, martingale);
*	merge to _events_ to pick up age at diag;


title 'martingale fit plot against age at diagnosis';
proc loess data = _events_ plots(only maxpoints = none)  = fit ;
	model mg1=age ;
run;

*	(c) Add linear effect of age to model;
%stpm2( age, df=5, scale=hazard, bhazard=rate, options = eform  );

*	(d) Excess mortality rate ratio as a function of age;
%predict(hr, hrnum = age:. , hrdenom = age:50, options = ci);

proc sort data = _events_ nodupkeys out = sorted; by age;run;
title 'Hazard ratio by age (reference age is 50)';
proc sgplot data = sorted;
	band x = age upper=hr_uci lower = hr_lci;
	series x = age y = hr;
	xaxis label='age at diagnosis';
	yaxis type=log label='excess hazard ratio';
	
run;
	
*	(e) Calculate martingale residuals ;	
%predict(mg2, martingale);

title 'martingale plot against age at diagnosis';
proc loess data = _events_ plots(only maxpoints = none)  = fit ;
	model mg2=age;
run;

*	(f) Generate splines for age and fit model; 		
%rcsgen (age, gen=rcsage, df=4, orthog=1);
%let age_knots = &save_knots.;
*	save orthog matrix and knots dataset;
proc datasets lib = work nolist force;
	delete agemat;
	change Tmat = agemat;
quit;
run;

%stpm2( rcsage1 rcsage2 rcsage3 rcsage4, scale=hazard, df=5, bhazard=rate);

%predict(mg3, martingale);

title 'martigale plot against age at diagnosis';
proc loess data = _events_ plots(only maxpoints = none)  = fit ;
	model mg3=age;
run;

*	g) Predicting hazard and survival functions;
%range(temptime,0, 5, 200);

*	set up a temporary dataset with a single row;
data single; age = 1; run;

*	macro to compute spline values for a specified age and store them in a series of macro strings;
%macro sp_set (ag, _v, _deg);
*	ag 		= scalar age to generate spline values for
	_v		= variable stub to hold spline values generated
	_deg	= number of variables to generate;
	
	%global &_v.1 &_v.2 &_v.3 &_v.4 &_v.5;

	%rcsgen(age, gen = &_v. ,  knots = &age_knots., tmatrix = agemat, set = single, scalar = &ag.);
	proc sql noprint;
		%do _ind = 1 %to &_deg;
			select &_v.&_ind. into :&_v.&_ind. from single ;
		%end;
	quit;
%mend;

*	macro to predict survival and hazard functions at specified ages;
%macro haz_surv;
	%do iage = 40 %to 80 %by 20;
		%sp_set(&iage., a, 4);
		%predict(haz&iage., hazard, 
			at= rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., timevar = temptime, per=1000);
		%predict(surv&iage., survival, 
			at= rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., timevar = temptime);
	%end;
%mend;

%haz_surv;

title 'Excess hazard by age';
proc sgplot data = _events_;	
	where age ^=.;
	series x = temptime y = haz40/legendlabel = '40';
	series x = temptime y = haz60/legendlabel = '60';
	series x = temptime y = haz80/legendlabel = '80';

xaxis label='years since diagnosis';
	yaxis type=log label='excess mortality rate (per 1,000 py)';
run;

title 'Relative survival by age at diagnosis';
proc sgplot data = _events_;	
	where age ^=.;
	series x = temptime y = surv40/legendlabel = '40';
	series x = temptime y = surv60/legendlabel = '60';
	series x = temptime y = surv80/legendlabel = '80';

	xaxis label='years since diagnosis';
	yaxis label='relative survival';
run;


*	(h) One year relative survival as a function of age;
data _events_;
	set _events_;
 	temptime=1;
run;

%predict (s1, survival, timevar = temptime , options=ci);

*	sort for plotting;
proc sort data = _events_ out = sorted;by age;run;

title 'One year relative survival as a function of age';
proc sgplot data = sorted;
	band x = age lower=s1_lci upper=s1_uci;	
	series x = age y = s1;
	xaxis label='age atdiagnosis';
	yaxis values=(0 to 1 by .2) label='1 year relative survival';
run;


*	(i) Five year relative survival as a function of age;
data _events_; set _events_;temptime = 5;run;

%predict (s5, survival, timevar = temptime, options=ci );

proc sort data = _events_ out = sorted;by age;run;

title 'Five year relative survival as a function of age';
proc sgplot data = sorted;
	band x = age lower=s5_lci upper=s5_uci;	
	series x = age y = s5;
	xaxis label='age atdiagnosis';
	yaxis values=(0 to 1 by .2) label='5 year relative survival';
run;

*	(j) Conditional relative survival;
*	survival at 5 years, among those who have survived 1 year;
data _events_;
	set _events_;
	conditional = s5/s1;
run;
	
proc sort data = _events_ out = sorted; by age; run;
title '5 year relative survival by age, for patients alive at 1 year';
proc sgplot data = sorted;
	series x = age y = conditional;
	xaxis label='age atdiagnosis';
	yaxis values=(0 to 1 by .2) label='5 year conditional relative survival ';
run;

*	(k) Obtain the hazard ratio as a funtion of age with age 50 as the reference;
%macro partpred2;
	%sp_set(50, a, 4);
	data hr_by_age;run;
	%do iage = 15 %to 90;
		%sp_set(&iage., s, 4);
		%predict(hr, hrnum = rcsage1:&s1. rcsage2:&s2. rcsage3:&s3. rcsage4:&s4. , 
			hrdenom = rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., ifp = _study_id_ = 1, options=ci);
		data hr_by_age; 
			set hr_by_age 
				_events_ (in=a keep = _study_id_ hr hr_uci hr_lci where = (_study_id_ = 1)); 
			if a then do;
				age = &iage.;
			end;
		run;
	%end;
%mend;


%partpred2;

title 'excess hazard ratio by age at diagnosis';
proc sgplot data = hr_by_age;
	band x = age upper=hr_uci lower = hr_lci/legendlabel = '95% CI';
	series x = age y = hr/ legendlabel = 'Excess Hazard ratio';
	xaxis  label='age at diagnosis';
	yaxis values = (.5 1 2 4 8) type=log  logbase= 2  label='excess hazard ratio (log scale)';
	
run;

*	(l) Sensitivity to the number of knots;
*	macro to evaluate the relationship between age and HR, evaluated with three different models
*	that differ only in the number of splines used for description of age
*	also retain the fit statistics to be able to report them together;
%macro knot_sens;
	data hr_by_age_k;run;	
	data fits; run;
	%do kns = 3 %to 5;
		%rcsgen (age, gen=rcsage, df=&kns., orthog=1);
		%let age_knots = &save_knots.;
		proc datasets lib = work nolist force;
			delete  agemat;
			change Tmat = agemat;
			quit;
		run;
		%let covstr = ;
		%do iat = 1 %to &kns.;
			%let covstr = &covstr. rcsage&iat.;
		%end;
		%stpm2( &covstr., scale=hazard, df=5, bhazard=rate, options = eform );

		data fits; set fits _fit_ (in=a); if a then df = &kns.;run;
		%sp_set(50, a, &kns.);
		%let refstr = ;
		%do iat = 1 %to &kns.;
			%let refstr = &refstr. rcsage&iat.:&&a&iat.;
		%end;

		%do iage = 15 %to 90;
			%let atstr = ;
			%sp_set(&iage, s, &kns.);
			%do iat = 1 %to &kns.;
				%let atstr = &atstr. rcsage&iat.:&&s&iat. ;
			%end;
			%predict(hr, hrnum = &atstr. , hrdenom = &refstr., ifp = _study_id_ = 1, options=ci);
			data hr_by_age_k; set hr_by_age_k 
				_events_ (in=a keep = _study_id_ hr hr_uci hr_lci where = (_study_id_ = 1)); 
				if a then do;
					age = &iage.;
					kns = &kns.;
				end;
			run;
		%end;
	%end;
%mend;

%knot_sens;

proc sort data = hr_by_age_k; by kns age;run;

title 'Excess hazard ratios by age at diagnosis';
title2 'for 3 different numbers of age splines';
proc sgpanel data = hr_by_age_k ;
	panelby kns /novarname  columns = 3;
	band x = age upper=hr_uci lower = hr_lci/ legendlabel = '95% CI';
	series x = age y = hr/legendlabel = 'Excess hazard ratio';
	colaxis label='age at diagnosis';
	rowaxis type=log logbase= 2 label='excess hazard ratio';
run;

*	compare fit statistics;;
*	ignore the keyword 'sum' in table;
proc tabulate data=fits format =8.1;
   class descr df;
   var value;

    table descr = '',df='degrees of freedom'*(value = '')
    	/ box = 'Fit Statistic';

   title 'comparison of fit statistics';
run;

/* Stop timer */
data _null_;
  dur = datetime() - &_timer_start;
  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';
run;
