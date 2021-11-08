/*
==================
 EXERCISE 121
 
 cox regression, compare with poisson regression; interactions

 lifetest; %smooth(); phreg; %lexis() ; genmod

reviewed:  4 May 2020
=====================
*/


options fmtsearch = (data.melanoma_formats);

*	truncate follow-up time to 10 years;
data melanoma10;
	set data.melanoma (where = (stage=1));
	
	entry = 0;			*	need an entry variable for the split using %lexis later;
	
	ca = status = 1;	*	deceased (status = 1)  censored (all other values);
	
	if surv_mm >120 then do;
		surv_mm = 120;
		if ca = 1 then ca = 0;
	end;
	
	surv_mm = surv_mm/12;	*	stata example has a scale parameter of 12;
run;

*	hazard plot by year8594;
ods graphics on;
proc lifetest notable data = melanoma10 
	outsurv=surv_by_year
	plots= (hazard loglogs);
	time surv_mm*ca(0);
	strata year8594;
run;
ods graphics off;

*	smoothed hazard plot on the log scale;
%smooth (data=surv_by_year, strat = year8594, time=surv_mm, survival=survival, yscale=log);


*	cumulative hazard plot (from stratified cox regression);
proc phreg data = melanoma10 noprint;
	model surv_mm*ca(0) = ;
	strata year8594;
	baseline out= ph_out cumhaz = ch;
run;

title 'cumulative hazard plot from Cox regression';
proc sgplot data = ph_out;
	series x = surv_mm y = ch /group = year8594;
	yaxis label='cumulative hazard';
	xaxis label='time since diagnosis';
run;
title ;

*	Cox regression;
proc phreg data = melanoma10;
	model surv_mm*ca(0) = year8594;
	hazardratio year8594;

run;


*	adjust for sex, calendar period and age;
proc phreg data = melanoma10;
	class agegrp year8594 sex (ref=last)/ref=first ;
	model surv_mm*ca(0) = year8594 sex agegrp;
	hazardratio year8594;
	output out = ph_test ressch = year8594_0 sex;
run;

proc sgplot data = ph_test;
	where ca=1;
	loess x = surv_mm y = year8594_0;
	refline 0 ;
	yaxis label='Schoenfeld residuals';
	xaxis label='time';
run;


*	similar for age group;

*	hazard plot ;
ods graphics on;
proc lifetest notable data = melanoma10 
	outsurv=surv_by_age
	plots = (hazard, loglogs);
	time surv_mm*ca(0);
	strata agegrp;
run;
ods graphics off;

*	smoothed hazard plot on the log scale;
%smooth (data=surv_by_age, strat= agegrp, time=surv_mm, survival=survival, yscale=log);


*	Cox regression;
proc phreg data = melanoma10;
	class agegrp (ref = first);
	model surv_mm*ca(0) = agegrp;
	hazardratio agegrp;
run;

*	adjust for sex, calendar period and age;
ods graphics off;
proc phreg data = melanoma10;
	class  year8594 sex agegrp /ref=first;
	model surv_mm*ca(0) = agegrp year8594 sex ;
	
*	test of PH assumption;
	assess  ph/resample ;
	output out = ph_test2 ressch=agegrp_1 agegrp_2 agegrp_3;
run;
ods graphics on;

proc sgplot data = ph_test2;
	where ca=1;
	loess x = surv_mm y = agegrp_1;
	refline 0 ;
run;

proc sgplot data = ph_test2;
	where ca=1;
	loess x = surv_mm y = agegrp_2;
	refline 0 ;
run;

proc sgplot data = ph_test2;
	where ca=1;
	loess x = surv_mm y = agegrp_3;
	refline 0 ;
run;

* 	Alternative 1: Cox regression using programming statements option;
proc phreg data = melanoma10;
	class  year8594 sex agegrp /ref=first;
	model surv_mm*ca(0) = agegrp agegrp*tvc year8594 sex/type3 ;
	tvc = (surv_mm>=2);
run;

*	alternative 2: split the data at 2 years;

%lexis(data = melanoma10, out = melanoma10_split,
		origin = 0,
		exit = surv_mm,
		breaks = %str(0, 2, 10),
		scale = 1,
		fail = ca,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,
		lrisk = ln_y,
		nint = fu);
		
proc print data = melanoma10_split noobs;
	where id <= 10 and _st = 1;
	var id _t0 surv_mm fu;
run;


proc phreg data = melanoma10_split;
	where _st = 1;
	class  year8594 sex agegrp (ref=first) fu;
	model surv_mm*ca(0) = agegrp  fu|agegrp year8594 sex/type3 ;
run;

* 	Effect of exposure for each level of the modifier;

proc phreg data = melanoma10_split;
	class  year8594 sex agegrp (ref=first) fu;
	model surv_mm*ca(0) =  agegrp fu*agegrp year8594 sex/type3 ;

run;


* Fit an analogous Poisson regression model;

*	truncate follow-up time to 10 years;
%lexis(data = melanoma10, out = melanoma10_finer,
		origin = 0,
		entry = entry,
		exit = surv_mm,
		breaks = 0 to 10 by 1,
		scale = 1,
		fail = ca,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,
		lrisk = ln_y,
		nint = fu);
		
data melanoma10_finer;
	set melanoma10_finer;
	where _st = 1;
	
	fuband = 0;
	if fu>2 then fuband = 1;
run;

proc genmod data = melanoma10_finer;
	class fu  agegrp sex year8594/ref=first;
	model ca = fu fuband agegrp sex year8594 fuband*agegrp/dist = p offset = ln_y;
	estimate "fuband0age45-59"  agegrp 1 0 0 -1/exp e;
	estimate "fuband0age60-74"  agegrp 0 1 0 -1/exp;
	estimate "fuband0age75+"    agegrp 0 0 1 -1/exp;
	
	estimate "fuband2age45-59" agegrp 1 0 0 -1 fuband*agegrp 1 0 0 -1 /exp e;
	estimate "fuband2age60-74" agegrp 0 1 0 -1 fuband*agegrp 0 1 0 -1/exp ;
	estimate "fuband2age75+"   agegrp 0 0 1 -1 fuband*agegrp 0 0 1 -1/exp;	
run;

