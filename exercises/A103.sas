/*
	A103 conditional net survival

	cohort version first, then the (slightly more complicated) period approach
	this uses the Pohar estimates, but would be the steps for the Ederer 2 estimates

	June 2020

	exercise from Paul Dickmans's online tutorial on conditional survival

	see  https://www.pauldickman.com/software/stata/  for full list of tutorials


*/


options fmtsearch = (data.colon_formats);

*	in this example dataset, we need the following variables:
	age			age in years at diagnosis
	sex			sex (coded as in the popmort file)
	dx			date of diagnosis
	yydx		year of diagnosis
	exit		date of death/censoring
	status		status at death/censoring;

data colon;
	set data.colon;
	id = _n_;				*	need an ID variable;
	surv = exit - dx;		*	survival time in days;
run;


/****************************************************************/
/*																*/
/*		cohort approach to analysis								*/
/*																*/
/****************************************************************/


*	net survival using life table methods
	monthly intervals to one year, anually afterwards;
title 'Net survival, Colon cancer, cohort approach';
title2 '(using dates)';
%rel_surv(infile = colon, 
	age 		= age, 
	origin 		= dx,
	exit 		= exit, 
	censor 		= status(0,4), 		/*	censored events	??????	*/
	intervals 	= %str(0 to 1 by 1/12,  2 to 10 by 1),
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right in(1 5)); * does this mean, return survival estimate at 1 and 5 years *******************?;

*	as above, using duration;
title 'Net survival, Colon cancer, cohort approach';
title2 '(using duration)';
%rel_surv(infile = colon, 
	age 		= age, 
	yydx 		= yydx,
	exit 		= surv, 
	censor 		= status(0,4), 		
	intervals 	= %str(0 to 1 by 1/12,  2 to 10 by 1),
	list		= pohar,
	popmort 	= data.popmort,
	patientid 	= id,
	crude		= right in(1 5));

*	ad hoc computation of conditional survival at 5 years, given survival to 1 year;
data cond_5_1;
	set grouped; * where was 'grouped' made *******************?;
	where right in (1 5);

	if right = 1 then cond = cr_p;
	retain cond;

	if right = 5 then do;
		cr_5_1 = cr_p/cond;
		output;
	end;
run;

title ' Net survival at 5 years, cohort approach, conditional on 1 year survival';
title2 'ad hoc method';

proc print data = cond_5_1 noobs label split = '*';
	var  cond cr_p cr_5_1;
	format cond cr_p cr_5_1 6.4;
	label cond 	= 'survival at*1 year'
		cr_p 	= 'survival at*5 years'
		cr_5_1	= 'survival at*5 years*given 1 year*survival';
run;
	
*	data setup for conditional life table analysis;
data colon_1;
	set colon;
	date_enter 		= dx+365.24;		*	1 year after diagnosis;
	time_enter		= 365.24;			*	ditto;
run;

title ' Net survival at 5 years, cohort approach, conditional on 1 year survival';
title2 '(using dates)';
%rel_surv(infile = colon_1  , 
	age 		= age, 
	origin 		= dx,
	exit 		= exit, 
	entry 		= date_enter,
	censor 		= status(0,4), 	
	intervals 	= 0 to 10 by 1,
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right = 5);

title ' Net survival at 5 years, cohort approach, conditional on 1 year survival';
title2 '(using duration)';
%rel_surv(infile = colon_1  , 
	age 		= age, 
	yydx		= yydx,
	exit 		= surv, 
	entry		= time_enter,
	censor 		= status(0,4), 	
	intervals 	= 0 to 10 by 1,
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right = 5);



/****************************************************************/
/*																*/
/*		period (1993 -1995 approach to analysis					*/
/*																*/
/****************************************************************/

data colon;
	set data.colon;
	id = _n_;				*	need an ID variable;
	surv = exit - dx;		*	survival time in days;

	period_start_date	= '01-jan-1993'd;				*	start of period window;
	period_start_time	= max(0, '01-jan-1993'd-dx);	*	patient enters risk at this time;
run;

*	net survival using life table methods
	monthly intervals to one year, anually afterwards;
title 'Net survival, Colon cancer, period approach (1993 - 1995)';
title2 '(using dates)';
%rel_surv(infile = colon, 
	age 		= age, 
	origin 		= dx,
	exit 		= exit, 
	entry		= period_start_date,
	censor 		= status(0,4), 		/*	censored events		*/
	intervals 	= %str(0 to 1 by 1/12,  2 to 10 by 1),
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right in(1 5));

*	as above, using duration;
title 'Net survival, Colon cancer, period approach (1993 - 1995)';
title2 '(using duration)';
%rel_surv(infile = colon, 
	age 		= age, 
	yydx 		= yydx,
	exit 		= surv, 
	entry		= period_start_time,
	censor 		= status(0,4), 		
	intervals 	= %str(0 to 1 by 1/12,  2 to 10 by 1),
	list		= pohar,
	popmort 	= data.popmort,
	patientid 	= id,
	crude		= right in(1 5));

*	ad hoc computation of conditional survival at 5 years, given survival to 1 year;
data cond_5_1;
	set grouped;
	where right in (1 5);

	if right = 1 then cond = cr_p;
	retain cond;

	if right = 5 then do;
		cr_5_1 = cr_p/cond;
		output;
	end;
run;

title ' Net survival at 5 years, period approach, conditional on 1 year survival';
title2 'ad hoc method';

proc print data = cond_5_1 noobs label split = '*';
	var  cond cr_p cr_5_1;
	format cond cr_p cr_5_1 6.4;
	label cond 	= 'survival at*1 year'
		cr_p 	= 'survival at*5 years'
		cr_5_1	= 'survival at*5 years*given 1 year*survival';
run;
	

*	data setup for conditional life table analysis (adjust entry into risk period);
data colon_1;
	set colon;
	adj_date_enter 		= max(period_start_date, dx+365.24);
	adj_time_enter		= max(period_start_time, 365.24);
run;

title ' Net survival at 5 years, period approach, conditional on 1 year survival';
title2 '(using dates)';
%rel_surv(infile = colon_1  , 
	age 		= age, 
	origin 		= dx,
	exit 		= exit, 
	entry 		= adj_date_enter,
	censor 		= status(0,4), 	
	intervals 	= 0 to 10 by 1,
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right = 5);

title ' Net survival at 5 years, period approach, conditional on 1 year survival';
title2 '(using duration)';
%rel_surv(infile = colon_1  , 
	age 		= age, 
	yydx		= yydx,
	exit 		= surv, 
	entry		= adj_time_enter,
	censor 		= status(0,4), 	
	intervals 	= 0 to 10 by 1,
	popmort 	= data.popmort,
	list		= pohar,
	patientid 	= id,
	crude		= right = 5);

