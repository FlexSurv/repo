/*
==================
 EXERCISE 204
 
 Period and cohort analysis II
 
%rel_surv 

reviewed:  30 Sept 2019
==================

*/

options fmtsearch = (data.melanoma_formats);


*	Load the Melanoma data, keep those with localized stage;

/*	make sure to only count person-time before 1984		*/
data m0;
	set data.melanoma (where = (stage=1 and yydx <= 1983));

	if exit > '31-dec-1983'd then do;
		status = 0;
		exit = '31-dec-1983'd;
	end;
run;


*	(a) : traditional cohort estimates;
title 'traditional cohort estimates, cases diagnosed prior to 1984';

%rel_surv(infile = m0,  patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);


*	(b) : traditional cohort estimates;
title 'traditional cohort analysis cases diagnosed 1977 - 1983';
%rel_surv(infile = m0 (where = (1977 <= yydx <=1983)), patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);

/* time since diagnosis as the timescale */ 
/* restrict person-time at risk to that within the period window (01jan1994-31dec1995) */
data m1;set data.melanoma (where = (stage=1 and yydx <= 1983));
	
	p_entry = '01-jan-1983'd;	*	start of period window (as a date);
	p_entry2 = '01-jan-1982'd;	*	start of period window (as a date);
	p_exit = '31-dec-1983'd;	*	end of period window (as a date);
		
*	censor observations at end of period window
	rel_surv manages the window start (late entry) internally;
	if exit > p_exit then do;
		exit = p_exit;
		status = 0;
	end;
run;

*	(c) : period estimates (period window 1 jan 1983 - 31 Dec 1983);
title 'period estimates (window is 1983 calendar year)';
%rel_surv(infile = m1, patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, entry = p_entry, scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);
	
*	(d) : period estimates (period window 1 jan 1982 - 31 Dec 1983);
title 'period estimates (window is 1982 - 1983 calendar years)';
%rel_surv(infile = m1, patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, entry = p_entry2,
	scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);

*	get melanoma data without year of diagnosis restriction;

data melanoma;
	set data.melanoma (where = (stage=1));
run;


*	(e) : Actual relative survival estimates (patients diagnosed in 1983);
title 'patients diagnosed in 1983';
%rel_surv(infile = melanoma (where = (yydx=1983)),
	age = age, sex=sex, origin=dx, exit = exit, patientid=id,
	scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);

*	(f) : Actual relative survival estimates (patients diagnosed in 1984);
title 'patients diagnosed in 1984';
%rel_surv(infile = melanoma (where = (yydx=1984)),
	age = age, sex=sex, origin=dx, exit = exit, patientid=id,
	scale = 365.24,
	censor = status(0,4), intervals = 0 to 15 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude=1);
	
*	(g) : Calculate all measures of relative survival and plot them on the same graph;
%macro surv_method;

	%do year = 1981 %to 1990;

*	cohort;
		data mel;set melanoma;
			exit_date = mdy(12,31,&year.);
			if exit > exit_date then do;
				exit = exit_date;
				status = 0;
			end;
		run;
		%rel_surv(infile = mel (where = (yydx<=&year.)), 
			age = age, sex=sex, origin=dx, exit = exit, patientid=id,
			scale = 365.24, censor = status(0,4), intervals = 0 to 15 by 1,
			popmort = data.popmort, mergeby = _age sex _year);

		%if &year. eq 1981 %then %let setfile = set grouped (in=a where = (fu= 5));
		%else %let setfile = set for_plot grouped (in=a where = (fu= 5));
	
		data for_plot;
			&setfile.;
			if a then do;
				method = 1;
				year = &year.;
			end;
			keep year method cr lo_cr hi_cr;
		run;
		
*	period;
		data mel;set melanoma;
			exit_date = mdy(12,31,&year.);
			entry_date = mdy(1,1,&year.);
			if exit > exit_date then do;
				exit = exit_date;
				status = 0;
			end;
		run;
		%rel_surv(infile = mel (where = (yydx<=&year.)), patientid=id,
			age = age, sex=sex, origin=dx, exit = exit, entry = entry_date,
			scale = 365.24, censor = status(0,4), intervals = 0 to 15 by 1,
			popmort = data.popmort, mergeby = _age sex _year);
	
		data for_plot;set for_plot grouped (in=a where = (fu= 5));
			if a then do;
				method = 2;
				year = &year.;
			end;
			keep year method cr lo_cr hi_cr;
		run;
		
*	Actual relative survival;
		%rel_surv(infile = melanoma (where = (yydx=&year.)),
			age = age, sex=sex, origin=dx, exit = exit, patientid=id,
			scale = 365.24, censor = status(0,4), intervals = 0 to 15 by 1,
			popmort = data.popmort, mergeby = _age sex _year);
	
		data for_plot;set for_plot grouped (in=a where = (fu= 5));
			if a then do;
				method = 3;
				year = &year.;
			end;
			keep year method cr lo_cr hi_cr;
			format method meth.;
		run;

	%end;
%mend;

proc format;
	value meth
	1 = 'Cohort'
	2 = 'Period'
	3 = 'Actual';
run;
	
%surv_method;

*	(h)	plot results;
proc sort data = for_plot (where =(method^=.));
	by method;
	run;
	
title '5-year survival';
proc sgplot data = for_plot ;
	by method ;
	scatter x = year y = cr/ 
		yerrorlower = lo_cr yerrorupper = hi_cr
		legendlabel='95% confidence interval';
	series x = year y = cr/legendlabel='Relative Survival';
	
		yaxis values=(.5 to 1 by .1) label='5-year RS';
		xaxis values = (1980 to 1990 by 2) label= 'Year of Diagnosis';

	format method meth.;

run;
	
proc sgpanel data = for_plot;
	panelby method /novarname columns = 3;
	scatter x = year y = cr/ 
		yerrorlower = lo_cr yerrorupper = hi_cr
		legendlabel='95% confidence interval';
	series x = year y = cr/legendlabel='Relative Survival';
	
	rowaxis values=(.5 to 1 by .1) label='5-year RS';
	colaxis values = (1980 to 1990 by 2) label= 'Year of Diagnosis';

run;	
	
