/*
==================
 EXERCISE 241
 
 Age standardization life table methods II

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

*	(a) estimates by time period;
%rel_surv(infile = data.melanoma, age = age, sex = sex, 
	exit = surv_mm, censor = status(0 4),
	yydx = yydx, scale = 12,
	intervals = 0 to 10 by 1,
	strat = year8594,
	popmort = data.popmort, patientid = id,
	crude =  1);
	

*	(b) crude and age standardised estimates by time period and age;
%rel_surv(infile = data.melanoma, age = age, sex = sex,
	exit = surv_mm, censor = status(0 4),
	yydx = yydx, scale = 12,
	intervals = 0 to 10 by 1,
	strat = agegrp year8594,
	popmort = data.popmort, patientid = id, crude = 1);
	
	
*	generate age group weights from distribution of cases;
proc freq data = data.melanoma noprint;
	table agegrp/out  = weights;
run;

*	merge weights on to interval estimates;
data rsr10;
	merge grouped weights;
	by agegrp;
	
	if right = 10;
	n0 = count;
	weight = percent/100;
	x = cr*weight;
run;

*	sort and print;
proc sort data = rsr10; by year8594;run;
	
proc tabulate data = rsr10;
	by year8594;
	class agegrp year8594;
	var x n0 weight cr;
	table agegrp = 'Age group' all, 
	(n0 cr weight x)* sum = ''  cr*mean  = 'Age standardised' ;
	format n0 5.0 cr weight x 5.3;
run;	
	

*	(c)  standardised rsr by user-supplied weights;
data melanoma;
	set data.melanoma;
	select (agegrp);
		when (0) iwei = 0.3039627;
		when (1) iwei = 0.2955711;
		when (2) iwei = 0.2927739;
		when (3) iwei = 0.1076923;
		otherwise;
	end;
run;

*	monthly intervals to 10 years followup;
*	save standardised estimates for later comparison, supress printed output;
%rel_surv(infile = melanoma, age = age, sex = sex, 
	exit = surv_mm, censor = status(0 4),
	yydx = yydx, scale = 12,
	intervals = 0 to 10 by 1/12,
	strat = year8594,
	popmort = data.popmort, patientid = id,
	standstrata = agegrp, weight_var = iwei,
	std_estimates = std_e2, 
	age_adj = 0);
	
*	Pohar-Perme Net survival estimates;
%rel_surv(infile = melanoma, age = age, sex = sex, 
	exit = surv_mm, censor = status(0 4),
	yydx = yydx, scale = 12,
	intervals =  0 to 10 by 1/12,
	strat = year8594,
	list = pohar,					/*	requests that Pohar estimate be standardised (default is E2)	*/
	popmort = data.popmort, patientid = id,
	standstrata = agegrp, weight_var = iwei,
	age_adj=0,
	std_estimates = std_pohar);

*	combine estimates from E2 and PP run for comparison;
data std_comp;
	merge std_e2 (rename = (as_rel = cr_e2 lo_rel = cr_e2_lo hi_rel = cr_e2_hi))
		std_pohar (rename = (as_rel = cns_pp lo_rel = cns_pp_lo hi_rel = cns_pp_hi));
	by year8594 interval;
		
	if first.year8594 then fu = 0;
	fu = fu+1;
	retain fu;
run;

*	print report
	labels saved from the standardisation are not specific to method of estimation;
proc print data = std_comp noobs split = '*' label;
	where mod(right,1) = 0;							/*	select only interger interval endpoints	*/
	by year8594;
	id year8594;
	var  right cns_pp cns_pp_lo cns_pp_hi cr_e2 cr_e2_lo cr_e2_hi;
	label right 		= 'End of*Interval'
		cns_pp 		= 'Pohar Net*survival'
		cns_pp_lo 	= 'Lower 95%*CI'
		cns_pp_hi 	= 'Upper 95%*CI'
		cr_e2 		= 'Ederer 2*Survival'
		cr_e2_lo 	= 'Lower 95%*CI'
		cr_e2_hi	= 'Upper 95%*CI';
run;

*	plot of estimates;
title 'comparison of Pohar and Ederer II standardised estimates';
proc sgpanel data = std_comp;
	panelby year8594/novarname;

	series x = right y = cns_pp/ legendlabel = 'Pohar-Perme';
	series x = right y = cr_e2/ legendlabel = 'Ederer II';

	rowaxis values = (.6 to 1 by .2) label = 'Age standardised survival proportion';
	colaxis label = 'Time since diagnosis (years)';
run;

*	(d)	plot SE estimates for comparison;
*	offset time point to be able to see the SE estimates separatel;
*	select annual intervals for displaly;
data se_plot;
	set std_comp;
	 where mod(right,1)=0;

	right_adj = right + .3;
run;

*	plot of se estimates;
title 'comparison of Pohar and Ederer II SE estimates at annual intervals';
title2 'PPE has been offset for visibility';
proc sgpanel data = se_plot;

	panelby year8594/novarname;

	scatter x = right_adj y = cns_pp/ 
		yerrorupper = cns_pp_hi yerrorlower = cns_pp_lo legendlabel = 'Pohar-Perme';
	scatter x = right y = cr_e2/ 
		yerrorupper = cr_e2_hi yerrorlower = cr_e2_lo legendlabel = 'Ederer II';

	rowaxis values = (.6 to 1 by .2) label = 'Age standardised survival proportion';
	colaxis label = 'Time since diagnosis (years)';
run;

