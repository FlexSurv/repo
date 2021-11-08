/*
==================
 EXERCISE 120
 
 cox regression, compare with poisson regression

	genmod
 	lifetest
 	phreg
 	sgplot

 	%lexis
 	%rcsgen

reviewed: 4 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

/* Data set used */

*	truncate follow-up time to 10 years;
data melanoma10;
	set data.melanoma (where = (stage=1));
		
	ca = status = 1;
	entry = 0;
	
	if surv_mm >120 then do;
		surv_mm = 120;
		if ca = 1 then ca = 0;
	end;
run;

*	Cox regression;
proc phreg data = melanoma10;
	model surv_mm*ca(0) = year8594;
	hazardratio year8594;
run;

*	log-rank test;
proc lifetest data = melanoma10 notable plots = none;
	time surv_mm*ca(0);
	strata year8594;
run;

*	covariate pattern for output of baseline hazard function; 
data cov;
	year8594 = 0;
	sex = 1;
	agegrp = 1;
	format agegrp agegrpa.;
run; 
* 	Cox regression including sex and age
	type3 statistics request Wald test;
proc phreg data = melanoma10;
	class agegrp (ref=first);
	model surv_mm*ca(0) = year8594 sex agegrp/ risklimits=wald type3;
	hazardratio year8594;
	hazardratio agegrp/diff = ref;
	hazardratio sex;

	baseline out = cox_out covariates = cov cumhaz = cumhaz;

	ods output parameterestimates = cox;
run;

*	log of the hazard function, estimated from differences in cumulative hazard;
*	save the mid-point of the interval as the survival time;
data cox_haz;
	set cox_out;
	if _n_ = 1 then last_cum = 0;
	retain last_cum;

	haz = cumhaz-last_cum;
	if haz > 0 then ln_cox = log(haz);
	last_cum = cumhaz;

	surv_mm = max(0, surv_mm-.5);

	keep surv_mm haz ln_cox;
run;

*	LR test of nested models;
proc phreg data = melanoma10;
	class agegrp (ref=first);
	model surv_mm*ca(0) = year8594  sex/risklimits=wald type3;
run;

*	enter -2LogLikelihood value from 'with covariates' column in results;
data lrt_pval;
	lrt = abs(15731.812-15588.962);
	df = 3;							*	in this case;
	p_value = 1 - probchi(LRT,df);
run;

proc print data=lrt_pval;
	title1 "LR test statistic and p-value";
run;


*	Comparison of Poisson and Cox regression models;
*split follow-up time to one-year intervals;
%lexis(data = melanoma10, out = melanoma10_split1,
		breaks = 0 to 120 by 12,      
		origin = 0,
		exit = surv_mm,
		scale = 1,
		fail = ca,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,
		lrisk = ln_y,
		nint = fu,
		st = _st);
	
		
proc genmod data = melanoma10_split;
	where _st = 1;
	class fu agegrp id sex year8594/ref=last;
	model ca = fu year8594 sex agegrp /dist=poisson offset=ln_y ;
	repeated subject = id;
	estimate "sex" sex 1  -1/exp;
	estimate "year8594" year8594 -1  1/exp;
	estimate "age 45-59" agegrp -1 1 0 0/exp;
	estimate "age 60-74" agegrp -1 0 1 0/exp;
	estimate "age 75 +" agegrp -1 0 0 1/exp;
	ods output estimates = poisson;
run;

data comp ;
	set poisson (keep = label lbetaestimate stderr 
		where = (substr(label,1,3) = 'Exp'))
		cox(in=a keep = parameter classval0 hazardratio stderr);
		
*	phreg does not transform the stderr for the hazard ratio,
	while stcox does.  Use the delta method to get a stderr estimate
	for the exponentiated HR;
	if a then stderr = stderr*hazardratio;

run;
proc print data = comp;

run;

*	split data into intervals defined by distinct failure times;
%lexis(data = melanoma10, out = melanoma10_fine,
		origin = 0,
		exit = surv_mm,
		breaks = failures,
		scale = 1,
		fail = ca,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,
		lrisk = ln_y,
		nint = fu,
		st = _st);

proc genmod data = melanoma10_fine;
	where _st = 1;
	class fu agegrp id sex year8594/ref=last;
	model ca = fu year8594 sex agegrp /dist=poisson offset=ln_y ;
	repeated subject = id;
	estimate "sex" sex 1  -1/exp;
	estimate "year8594" year8594 -1  1/exp;
	estimate "age 45-59" agegrp -1 1 0 0/exp;
	estimate "age 60-74" agegrp -1 0 1 0/exp;
	estimate "age 75 +" agegrp -1 0 0 1/exp;
	ods output estimates = poisson_fine;
run;

data comp2 ;
	set poisson (keep = label lbetaestimate stderr 
		where = (substr(label,1,3) = 'Exp'))
		poisson_fine (keep = label lbetaestimate stderr 
		where = (substr(label,1,3) = 'Exp'))
		cox(in=a keep = parameter classval0 hazardratio stderr);
		
*	phreg does not transform the stderr for the hazard ratio,
	while stcox does.  Use the delta method to get a stderr estimate
	for the exponentiated HR;
	if a then stderr = stderr*hazardratio;

run;
title "Hazard ratios and standard errors for various models";
proc print data = comp;
run;
	
title ;
*	split data at one-month intervals;
%lexis(data = melanoma10, out = melanoma_monthly,
		breaks = 0 to 120 by 1,      
		origin = 0,
		exit = surv_mm,
		scale = 1,
		fail = ca,
		right = _t,
		cens = 0,
		risk = y,
		lrisk = ln_y,
		nint = fu,
		st = _st);

data melanoma_monthly;
	set melanoma_monthly (where = (_st = 1));
	fu = fu-1;
run;
		
*	generate spine variables for time, rather than estimating so many paramters;
%rcsgen(fu, df = 4, set = melanoma_monthly, orthog=1);

proc genmod data = melanoma_monthly;
	class  agegrp  sex year8594 id/ref=last;
	model ca = rcs1-rcs4 year8594  sex agegrp /dist=poisson offset=ln_y ;
	repeated subject=id;
	output out = melanoma_pred xbeta = xb ;
run;

*	do we really need to plot each point 100 times over? No;
proc sort data = melanoma_pred 
	(where =(year8594 = 0 and sex = 1 and agegrp = 1))
	out = for_plot (keep = xb ln_y fu) nodupkeys;
	by fu;
run;

*	append baseline plot from same PHREG model for the covariate
	pattern chosen;
data for_plot;
	set for_plot (in=a) cox_haz;
	if a then ln_haz = xb-ln_y;		*	offset is included in the linear predictor from genmod;
run;

title 'spline fit of log hazard function vs smoothed Cox baseline log hazard';
proc sgplot data = for_plot;
	series y = ln_haz x = fu/ legendlabel = 'Spline fit';
	loess y = ln_cox x = surv_mm/nomarkers lineattrs = (pattern = dot)
		legendlabel = 'Smoothed Cox';
	yaxis label='Log Hazard';
	xaxis values=(0 to 150 by 50) label='months since diagnosis';
run;

		
