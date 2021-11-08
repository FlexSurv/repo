/*
==================
 EXERCISE 201
 
 Life table estimates of relative survival  

reviewed:  4 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;
data melanoma;set data.melanoma (where = (stage=1)); run;

* a) compute life table survival estimates with a variety of interval requests;
title 'life table estimates, annual intervals';
*	Lifetable estimates, annual intervals;
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat = year8594, crude=1);
	
title 'life table estimates, 6 month intervals';
*	(b) Lifetable estimates, 6 month intervals;
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), intervals = 0 to 10 by .5,
	popmort = data.popmort, mergeby = _age sex _year,  
	strat = year8594, crude=1);
	
title 'life table estimates, 3 month intervals for first year, followed by annual';
*	(c) Lifetable estimates, 3 month intervals for first year then annual;
%rel_surv(infile = melanoma, patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), 
	intervals = %str(0, 0.25, 0.5, 0.75, 1 to 10 by 1),
	popmort = data.popmort, mergeby = _age sex _year, 
	strat = year8594, crude=1);
	
title 'life table estimates, annual intervals to 20 years';
*	(d) Lifetable estimates, annual intervals up to 20 years;
%rel_surv(infile = melanoma, patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), intervals = 0 to 20 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat = year8594, crude=1);
	
*	 (e) Plot estimates of cumulative relative survival;
title 'cumulative relative survival (E II)';
proc sgplot data = grouped;
	series x = right y = cr/group=year8594 markers;
	xaxis label = 'End of interval';
	yaxis label = 'Cuulative Reiative Survival (Ederer II)';
run;

*	Alternative approach for producing a graph of estimates;
proc sgpanel data = grouped;
	panelby year8594 /novarname ;
	scatter x = right y = cr/
		yerrorlower = lo_cr yerrorupper = hi_cr;
		colaxis label = 'Years from diagnosis';
		rowaxis label = 'Relative survival'
			values=(.4 to 1.0 by .2);
run;

*	Another alternative approach for producing a graph of estimates;
proc sgplot data = grouped;

	scatter x = right y = cr/group=year8594
		yerrorlower = lo_cr yerrorupper = hi_cr;
		xaxis label = 'Years from diagnosis';
		yaxis label = 'Relative survival'
			values=(.4 to 1.0 by .2);
run;

*	(f) Plot estimates of interval-specific relative survival;
title 'interval-specific relative survival (Ederer 2)';
proc sgplot data = grouped;
	series x = right y = r/group = year8594 markers;
	yaxis label='Interval-specific relative survival';
	xaxis label='End of interval';
run;

title 'interval-specific net survival (Pohar-Perme)';
proc sgplot data = grouped;
	series x = right y = ns_w/group = year8594 markers;
	yaxis label='Interval-specific relative survival';
	xaxis label='End of interval';
run;

*	(g) Comparing 2 approaches to estimating expected (and relative) survival;
title "using the 'list' option to compare E2 with Pohar (net) survival"; 
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, scale = 365.24,
	censor = status(0,4), intervals = 0 to 20 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	list = l d w cr cr_p,
	strat = year8594, crude=1);

*	(h) Estimating rates (supress printed output);
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, crude = 0);
	
data rates;
	set grouped;
		
*	first use an approximation to person-time at risk;
	obs_rate1 = d / (l-d/2-w/2);
	
*	now using exact person-time;
	obs_rate2 = d / y;
	
*	now estimate the rate by transforming the estimated survival probability;
	obs_rate3 = -log(p);
	
*	now the probability of death (which is not strictly a rate);
	q=1-p;
	
	format obs_rate1 obs_rate2 obs_rate3 q 6.4;
	
*	Estimate the expected mortality rate;
	exp_rate = d_star/y;
	
*	Estimate the excess mortality rate;
	excess_rate = (d-d_star)/y;
run;

title 'estimated interval rates using different methods';
proc print data = rates noobs split = '*';
	var right d w y obs_rate1 obs_rate2 obs_rate3 q;
	label right = 'end-point*of interval'
		d = 'deaths'
		w = 'censored'
		y = 'person-time'
		obs_rate1 = 'actuarial*(approximate*person-time)'
		obs_rate2 = 'exact*person-time'
		obs_rate3 = 'hazard*transformation'
		q = 'death*probability';
run;

title 'excess mortality (Ederer II approach)';
proc print data = rates noobs split = '*';
	var right d d_star y obs_rate2 exp_rate excess_rate;
	label right = 'end-point*of interval'
		d = 'deaths'
		d_star = 'expected*deaths'
		y = 'person-time'
		obs_rate2 = 'exact*person-time'
		exp_rate = 'expected*mortality*rate'
		excess_rate = 'excess*mortality*rate';

run;

*	Plot excess mortality rate as a function of time since diagnosis;
title 'excess mortality';
proc sgplot data = rates;
	series x = right y = excess_rate;
	xaxis label="Time since diagnosis"
		values =(0 to 10 by 2);
	yaxis label="Excess rate"
		values=(0 to .05 by .01);
run;

*	compute attained age (CPAC macro does not retain this value);
data _cases_;
	set _cases_;
	
*	attained age at death or censoring;
	_age= floor(age + input(scan(interval,1,' '),5.2));

	label _age = 'age at death/censoring';
	
run;

title 'attained age (at death or censoring) compared to age at diagnosis';
proc means data = _cases_ mean maxdec = 1;
	class right;
	var age _age;
run;
	
	
