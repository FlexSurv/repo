/*
==================
 EXERCISE 240
 
 Age standardization life table methods I
 
 tabulate print plot
 
 %rel_surv

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

* 	load melanoma data and assign weights;
data melanoma;
	set data.melanoma (where = (stage=1));
	select (agegrp);
		when (0) iwei = 0.2751034;
		when (1) iwei = 0.296164;
		when (2) iwei = 0.2888304;
		when (3) iwei = 0.1399022;
		otherwise;
	end;
run;

*	(a) Crude estimates duration in months;
%rel_surv(infile = melanoma, age = age, sex = sex, patientid=id,
	yydx = yydx, exit = surv_mm, scale = 12, censor = status(0 4),
	intervals = 0 to 15 by 1,
	popmort = data.popmort, 
	crude = 1);
	
*	crude estimates using dates of diagnosis and death/censoring;
%rel_surv(infile = melanoma, age = age, sex = sex, patientid=id,
	origin = dx,  exit = exit, censor = status(0 4),
	intervals = 0 to 15 by 1,
	popmort = data.popmort, 
	crude = 1);

	
*	(b) Age-specific estimates (using duration), stratified by age group;
%rel_surv(infile = melanoma, age = age, sex = sex, patientid=id,
	yydx = yydx, exit = surv_mm, scale = 12, censor = status(0 4),
	intervals = 0 to 15 by 1,
	strat = agegrp,
	popmort = data.popmort, 
	crude = 1);

*	Age-standardised 10-year RSR 'by hand';
proc freq data = melanoma noprint;
	table agegrp/out  = weights;
run;

data rsr10;
	merge grouped weights;
	by agegrp;
	
	if right = 10;
	n0 = count;
	weight = percent/100;
	x = cr*weight;
run;
	
proc tabulate data = rsr10;
	class agegrp;
	var x n0 weight cr;
	table agegrp = 'Age group' all, 
	(n0 cr weight x)* sum = ''  cr*mean  = 'Age standardised' ;
	format n0 5.0 cr weight x 5.3;
run;	

*	(c) Age-standardised with user-supplied weights;
%rel_surv(infile = melanoma, age = age, sex = sex, patientid=id,
	yydx = yydx, exit = surv_mm, scale = 12, censor = status(0 4),
	intervals = 0 to 10 by 1/12, 
	weight_var = iwei, standstrata = agegrp, 
	age_adj = mod(right,1)=0,
	popmort = data.popmort);
	
*	(d) Pohar Perme estimate;
*	truncate at 10 years;
data melanoma10;
	set melanoma;

	if (exit-dx)/365.5 > 10 then do;
		status = 0;
		exit = dx+ 10*365.5;
	end;
run;

*	lifetable estimates of Ederer 2 and Pohar
	restriction on printing selects intervals with
	the integer endpoints (1, 2, ... 10)
	use of the 'list' parameter to select specific output fields;
%rel_surv(infile = melanoma, age = age, sex = sex,
	origin = dx,  exit = exit, censor = status(0 4),
	intervals = 0 to 10 by 1/12, 
	list = right cr cr_p,
	popmort = data.popmort, patientid = id,
	crude =  mod(right,1)=0);
	
*	corresponding plot of all interval estimates;
title 'life table estimates of relative survival';
proc sgplot data = grouped;
	series x = right y = cr/legendlabel='Ederer II';
	series x = right y = cr_p/legendlabel='Pohar Perme';
	yaxis values = (.7 to 1.0 by .05) label = 'cumulative relative survival' ;
	xaxis label = 'Time since diagnosis (years)';
run;
	
	
