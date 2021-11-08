/*
==================
 EXERCISE 203
 
 Period and cohort analysis  I


reviewed:  4 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);


*	Load the Melanoma data, keep those with localized stage;

*	PERIOD ANALYSIS;
/* time since diagnosis as the timescale */ 
/* restrict person-time at risk to that within the period window (01jan1994-31dec1995) */
data period;
	set data.melanoma (where = (stage=1)); 
	
	p_entry = '01-jan-1994'd;	*	start of period window (as a date);
	p_exit = '31-dec-1995'd;	*	end of period window (as a date);
		
*	censor observations at end of period window
	rel_surv manages the window start (late entry) internally;
	if exit > p_exit then do;
		exit = p_exit;
		status = 0;
	end;
	
run;

*	period analysis 1994 - 1995;
title 'period analysis 1994 - 1995';
%rel_surv(infile = period, patientid=id,
	age = age, sex=sex, origin=dx, 
	entry = p_entry, 
	exit = exit, scale = 365.24,
	censor = status(0,4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat = sex, crude=1);
	
	
*	complete ;

*	complete analysis 1994 - 1995;
title 'complete analysis 1994 - 1995';
%rel_surv(infile = melanoma (where = (stage=1)), patientid=id,
	age = age, sex=sex, origin=dx, 
	exit = exit, scale = 365.24, 
	censor = status(0, 4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat = sex, crude=1);
	
