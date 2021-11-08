/*
==================
 EXERCISE 232
 
 Modelling excess mortality using flexible parametric models II
 
 sgplot sgpanel sql
 %stset %stpm2 %rcsgen %predict

reviewed: 6 May 2020
====================


*/

/* Start timer */
/*%let _timer_start = %sysfunc(datetime());*/

options fmtsearch = (data.colon_formats);

*	(a) Load Data, merge expected mortality;
data colon;
	set data.colon;
	where age <= 90;
	
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

proc sql;
	create table colon_mm as select a.*,
		b.rate
		from colon a left join data.popmort b
		on floor(a.age + a.surv) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv) = b._year
		order by id;
quit;

%stset(colon_mm, status(1 2), surv, id);

*	(b) Fit flexible parametric model using splines for age;		
%rcsgen (age, gen=rcsage, df=4, orthog=1);
%let age_knots = &save_knots.;
proc datasets lib = work nolist force;
	delete agemat;
	change Tmat = agemat;
quit;
run;

%stpm2( rcsage1 rcsage2 rcsage3 rcsage4, scale=hazard, df=5, bhazard=rate);
	
*	(c) Time-dependent effect of age 2df for each age spline, total of 8 new parameterss;
%stpm2( rcsage1 rcsage2 rcsage3 rcsage4, scale=hazard, df=5, bhazard=rate,
	tvc = rcsage1 rcsage2 rcsage3 rcsage4, dftvc = 2);
	
*	LR test of nested models;
*	enter -2LogLikelihood value ;
data lrt_pval;
	lrt = abs(36107-35805);
	df = 8;							*	in this case;
	p_value = 1 - probchi(LRT,df);
run;

title 'likelihood ratio test';
proc print data=lrt_pval;
	title1 "LR test statistic and p-value";
run;


title ;
* 	(d) Predicing hazard and survival functions;
%range(temptime,0, 5, 201);

data single; age = 1; run;
%macro sp_set (ag, _v, _deg);
*	ag 		= scalar age to generate spline values for
	_v		= variable stub to hold spline values generated
	_deg	= number of variables to generate;
	
	%global &_v.1 &_v.2 &_v.3 &_v.4 &_v.5;

	%rcsgen(age, gen = &_v. ,  knots = &age_knots., tmatrix = agemat, 
		set = single, scalar = &ag., orthog=1);
	proc sql noprint;
		%do _ind = 1 %to &_deg;
			select &_v.&_ind. into :&_v.&_ind. from single ;
		%end;
	quit;
%mend;

%macro haz_surv;
	%do iage = 40 %to 80 %by 20;
		%sp_set(&iage, a, 4);

		%predict(haz&iage., hazard, 
			at= rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., timevar = temptime, per=1000);
		%predict(surv&iage., survival, at= rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., timevar = temptime);
		
	%end;
%mend;

%haz_surv;


title 'Excess mortality for selected ages at diagnosis';
proc sgplot data = _events_;	
	where age ^=.;
	series x = temptime y = haz40/legendlabel = '40';
	series x = temptime y = haz60/legendlabel = '60';
	series x = temptime y = haz80/legendlabel = '80';

	xaxis label='years since diagnosis';
	yaxis type=log label='excess mortality rate (per 1,000 py)';
run;

title 'Relative survival for selected ages at diagnosis';
proc sgplot data = _events_;	
	where age ^=.;
	series x = temptime y = surv40/legendlabel = '40';
	series x = temptime y = surv60/legendlabel = '60';
	series x = temptime y = surv80/legendlabel = '80';

	xaxis label='years since diagnosis';
	yaxis label='relative survival';
run;


*	(e) One year relative survival as a function of age;
*	this forces avery one of the 15,000+ records to have their 1-year survival estimated
	if there are 100 patients with identical covariate patters, then you are estimating
	the same value 99 times unnecessarily;

data _events_;
	set _events_;
 	temptime=1;
run;

%predict (s1, survival, timevar = temptime, options = ci );

proc sort data = _events_ out = sorted;by age;run;

title 'Relative survival at 1 year by age at diagnosis';
proc sgplot data = sorted;
	band x = age lower=s1_lci upper=s1_uci;	
	series x = age y = s1;
	xaxis label='age atdiagnosis';
	yaxis values=(0 to 1 by .2) label='1 year relative survival';
run;


*	(f) Five year relative survival as a function of age;
data _events_; set _events_;temptime = 5;run;

%predict (s5, survival, timevar = temptime, options = ci );


proc sort data = _events_ out = sorted;by age;run;

title 'Relative survival at 5 years by age at diagnosis';
proc sgplot data = sorted;
	band x = age lower=s5_lci upper=s5_uci;	
	series x = age y = s5;
	xaxis label='age atdiagnosis';
	yaxis values=(0 to 1 by .2) label='5 year relative survival';
run;

*	(g) Conditional relative survival
	parameter 'tcond' gives the conditioning time (in years);

%predict (s15, survival ,  timevar = temptime, tcond = 1, options = ci );

proc sort data =_events_ out = sorted nodupkey; by age;run;

title 'Relative survival at 5 years, conditional on 1 year survival';
title2 'by age, obtained by estimation';
proc sgplot data = sorted;

	band x = age upper = s15_uci lower = s15_lci/legendlabel = '95% CI';	
	series x = age y = s15/legendlabel = '5 year survival (conditional 1 year)';

	xaxis label='age at diagnosis';
	yaxis values=(0 to 1 by .2) label='5 year conditional relative survival ';
run;

*	(h) Obtain the hazard ratio as a funtion of age with age 50 as the reference;
*	note the new temporary time variable;

%macro partpred2;
%range(tempt, 0, 5, 101);
	%sp_set(50, a, 4);

	%let age_list = 40 60 70 80;
	
	%local _ind next_name;
	%do _i= 1 %to %sysfunc(countw(&age_list.));
   		%let iage = %scan(&age_list., &_i.);
		%sp_set(&iage, s, 4);
		%predict(hr&iage., hrnum = rcsage1:&s1. rcsage2:&s2. rcsage3:&s3. rcsage4:&s4. , 
		hrdenom = rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., 
		timevar = tempt, options = ci);
	%end;
%mend;

%partpred2;

title 'excess hazard ratio by age, relative to patients aged 50';
*	macro to plot the 4 graphs;
%macro plot_four;
	%let age_list = 40 60 70 80;
	
	%do _i= 1 %to %sysfunc(countw(&age_list.));
   		%let iage = %scan(&age_list., &_i.);
		proc sgplot data = _events_ ;
			band x = tempt upper=hr&iage._uci lower = hr&iage._lci;
			series x = tempt y = hr&iage.;
			xaxis label='time since diagnosis';
			yaxis label='excess mortality rate ratio'
				type=log;
		run;
	%end;
%mend;

%plot_four;

*	(i) survival differences and hazard differences;
%macro partpred3;
	%sp_set(50, a, 4);
	data hd_by_age;run;
	data sd_by_age;run;
	%let age_list = 40 60 70 80;
	
	%local _ind next_name;
	%do _i= 1 %to %sysfunc(countw(&age_list.));
   		%let iage = %scan(&age_list., &_i.);
		%sp_set(&iage, s, 4);
		%predict(hd, hdiff1 = rcsage1:&s1. rcsage2:&s2. rcsage3:&s3. rcsage4:&s4. , 
		hdiff2 = rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., per=1000,
		timevar = tempt, options = ci);
		%predict(sd, sdiff1 = rcsage1:&s1. rcsage2:&s2. rcsage3:&s3. rcsage4:&s4. , 
		sdiff2 = rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., 
		timevar = tempt, options = ci);
		data hd_by_age; set hd_by_age _events_ (in=a); 
			if a then do;
				age = &iage.;
			end;
		run;
		data sd_by_age; set sd_by_age _events_ (in=a); 
			if a then do;
				age = &iage.;
			end;
		run;
	%end;
%mend;


%partpred3;

title 'excess mortality differences for seleted ages compared to reference (50 years)';
proc sgpanel data = hd_by_age;
	panelby age;
	band x = tempt upper=hd_uci lower = hd_lci;
	series x = tempt y = hd;
	colaxis label='time since diagnosis';
	rowaxis label='difference in excess mortality rate (per 1,000 p-yr)'
		values =(-100 0 100 200 400 600 800);
run;

title 'relative survival differences for selected ages, compared to reference (50 years)';
proc sgpanel data = sd_by_age;
	panelby age;
	band x = tempt upper=sd_uci lower = sd_lci;
	series x = tempt y = sd;
	colaxis label='time since diagnosis';
	rowaxis label='difference in relative survival'
		values =(-0.25 to 0.05 by .05);
run;

*	(j) Sensitivity to the number of knots for time-dep effects;
*	be sure to use the 'tempt' time variable, as 'temptime' is non-missing for all subjects;
%macro knot_sens;
	data hr_by_k;run;	
	data fits; run;
	%sp_set(70, s, 4);
	%sp_set(50, a, 4);
	%do kns = 1 %to 3;
		%stpm2(rcsage1 rcsage2 rcsage3 rcsage4, scale=hazard, df=5, bhazard=rate,
			tvc = rcsage1 rcsage2 rcsage3 rcsage4, dftvc=&kns.);
		data fits; set fits _fit_ (in=a); if a then df = &kns.;run;

		%do _i= 1 %to %sysfunc(countw(&age_list.));
   			%let iage = %scan(&age_list., &_i.);
			%predict(hr, hrnum = rcsage1:&s1. rcsage2:&s2. rcsage3:&s3. rcsage4:&s4.  , 
				hrdenom = rcsage1:&a1. rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., 
				timevar = tempt, options = ci);
			data hr_by_k; set hr_by_k _events_ (in=a); 
				if a then do;
					kns = &kns.;
				end;
			run;
		%end;
	%end;
%mend;

%knot_sens;

title 'excess hazard ratio age 70 vs age 50, from different versions of TVC models';
proc sgpanel data = hr_by_k;
	panelby kns /columns = 3 novarname;
	band x = tempt upper=hr_uci lower = hr_lci;
	series x = tempt y = hr;
	colaxis label='time since diagnosis';
	rowaxis type=log label='excess hazard ratio';
run;

title 'excess hazard ratio age 70 vs age 50, from different versions of TVC models';
proc sgplot data = hr_by_k;
	series x = tempt y = hr/group = kns;
	series x = tempt y = hr_lci/group = kns lineattrs=(pattern=dash);	
	series x = tempt y = hr_uci/group = kns lineattrs=(pattern=dot);	
	xaxis label='time since diagnosis';
	yaxis type=log label='excess hazard ratio';
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

/*	stop timer and report duration	*/
/*data _null_;*/
/*  dur = datetime() - &_timer_start;*/
/*  put 30*'-' / ' TOTAL DURATION:' dur time13.2 / 30*'-';*/
/*run;*/
