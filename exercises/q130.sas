/*
==================
 EXERCISE 130
 
 Understanding splines
 genmod, %lexis(), sgplot , means, %rcsgen

reviewed:  5 sep 2020
=====================
*/

options fmtsearch = (data.melanoma_formats);

proc format;
	value status
		1 = 'failed'
		0 = 'censored';

	value sex
		0 = 'Male'
		1 = 'Female';
run;
*	(a) truncate follow-up time to 10 years;
data melanoma10;
	set data.melanoma;
	
*	create censoring variable;
	d = 0;							*	alive (censored);
	if status in(1 2) then d = 1;	*	deceased;
	
*	truncate followup to 10 years (120 months);
	if surv_mm >=120 then do;
		surv_mm = 120;
		d = 0;
	end;
	
	entry = 0; 						*	need an entry variable.  set it to origin;
	female = 0;
	if sex = 2 then female = 1;
	
	label d = 'status at exit';
	format d status.
		female sex.;
run;

*	Split the data with narrow (1 month) time intervals;
%lexis(data = melanoma10, 
		out = melanoma10_monthly,
		origin = 0,
		entry = entry,
		exit = surv_mm,
		breaks = 0 to 10 by 1/12,
		scale = 12,
		fail = d,
		cens = 0,
		right = _t,
		nint = fu,
		left = _t0,
		risk = risktime,
		lrisk = ln_y);
	
proc means data = melanoma10_monthly 
	(where = (_st = 1))
	noprint nway;
	class fu female year8594 agegrp;
	output out = collapsed (drop = _freq_ _type_)
				sum(d risktime) = deaths risktime
				min(_t0) = start
				max(_t) = end;
run;

*	assign the center of the interval as followup time
	log of total person-time at risk needed for rate regression;
data collapsed;
	set collapsed;
	midtime = (end+start)/2;
	ln_y = log(risktime);
run;

title 'examine collapsed data';
proc freq data = collapsed;
	where fu in (1 5 10);
	table fu;
run;

title 'observations for tenth followup interval';
proc print data = collapsed noobs;
	where fu = 10;
	format deaths 4.;
run;
	
* 	(b) Fit a model with a parameter for each interval;
title 'fit followup intervals (baseline hazard estimates)';
proc genmod data = collapsed ;
	class fu;
	model deaths = fu/dist = poisson offset=ln_y noint;
	output out = haz_grp pred = pred_deaths;
run;

*	(c) predict the baseline, one observation for each interval,
	and scale the rate to represent deaths / 1,000 person-years;
data haz_grp ;
	set haz_grp;
	haz_grp = 1000* pred_deaths/risktime;
run;

title 'scatter plot of interval hazard estimates';
proc sgplot data = haz_grp;
	scatter x = midtime y = haz_grp;
	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

		
*	(d) 2 linear splines (1 knot at 1.5 years);
data collapsed;
	set collapsed;
	
*	define position of knot (an intercept at 1.5);
	lin_int 					= 0;
	if midtime>1.5 then lin_int = 1;

	lin_s1 		= midtime;
	lin_s2 		= (midtime - 1.5)*(lin_int);
run;

*	Fit two separate linear regression lines (4 parameters);
proc genmod data = collapsed;
	model deaths = lin_s1 lin_int lin_s2/dist = poisson offset=ln_y;
	output out = lin_haz pred = pred_deaths;
run;

data lin_haz;
	set lin_haz (in=a)
	haz_grp;

	lin_haz1 = .;
	lin_haz2 = .;
	if a then do;
		if midtime <= 1.5 then lin_haz1 = 1000* pred_deaths/risktime;
		else lin_haz2 = 1000* pred_deaths/risktime;
	end;
run;

title 'two linear splines plus separate intercept';
proc sgplot data = lin_haz;
	scatter x = midtime y = haz_grp;
	series x = midtime y = lin_haz1;
	series x = midtime y = lin_haz2;
	
	refline 1.5 /axis=x lineattrs=(pattern=dash);

	yaxis
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

*	(e) Force the functions to join at the knot (3 parameters);		
proc genmod data = collapsed;
	model deaths = lin_s1 lin_s2/dist = poisson offset=ln_y;
	output out = lin_haz2 pred = pred_deaths;
run;

data lin_haz2;
	set lin_haz2 (in=a)
	haz_grp;
	if a then lin_haz1 = 1000* pred_deaths/risktime;
run;

title 'two linear spines forced to join at knot';
title2 '(no separate intercept)';
proc sgplot data = lin_haz2;
	scatter x = midtime y = haz_grp;
	series x = midtime y = lin_haz1;
	refline 1.5 /axis=x lineattrs=(pattern=dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

*	(f) Now use cubic polynomials with 1 knot at 2 years;
data collapsed;
	set collapsed;
	
*	0 for times <= 2, 1 for times greater than 2;
	cubic_int 					= 0;		
	if midtime>2 then cubic_int = 1;

	cubic_s1 	= midtime;	
	cubic_s2 	= midtime**2;
	cubic_s3 	= midtime**3;

	cubic_lin 	= (midtime - 2)*(cubic_int);
	cubic_quad 	= ((midtime - 2)**2)*(cubic_int);
	cubic_s4 	= ((midtime - 2)**3)*(cubic_int);
run;

proc genmod data = collapsed;
	model deaths = cubic_s1-cubic_s4 cubic_int cubic_lin cubic_quad /dist=poisson offset=ln_y;
	output out = lin_haz3 pred = pred_deaths;
run;

data lin_haz3;
	set lin_haz3 (in=a)
	haz_grp;
	
	haz_cubic1 =.;
	haz_cubic2 = .;
	if a then do;
		if midtime <=2 then haz_cubic1 = 1000* pred_deaths/risktime;
		else haz_cubic2 = 1000* pred_deaths/risktime;
	end;
run;

title 'two cubic polynomials with knot at 2 years';
title2 'note discontinuty at knot'; 
proc sgplot data = lin_haz3;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic1;
	series x = midtime y = haz_cubic2;
	refline 2 /axis=x lineattrs=(pattern=dash);

	yaxis
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

*	(g) constrain to join at knots (drop separate intercept);	
proc genmod data = collapsed;
	model deaths = cubic_s1-cubic_s4 cubic_lin cubic_quad /dist=poisson offset=ln_y;
	output out = lin_haz4 pred = pred_deaths;
run;

data lin_haz4;
	set lin_haz4 (in=a)
	haz_grp;
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;
		
run;

title 'two cubic polynomials with knot at 2 years';
title2 'polynomials joined at knot'; 
proc sgplot data = lin_haz4;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic;
	refline 2 /axis=x lineattrs=(pattern=dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;
		
*	(h) continuous 1st derivative (drop second linear term);	
proc genmod data = collapsed;
	model deaths = cubic_s1-cubic_s4 cubic_quad /dist=poisson offset=ln_y;
	output out = lin_haz5 pred = pred_deaths;
run;

data lin_haz5;
	set lin_haz5 (in=a)
	haz_grp;
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;
		
run;

title 'two cubic polynomials with knot at 2 years';
title2 'polynomials joined at knot continuous in 1st derivative'; 
proc sgplot data = lin_haz5;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic;
	refline 2 /axis=x lineattrs=(pattern=dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;
	
*	(i) continuous 2nd derivative (drop second quadratic term);
proc genmod data = collapsed;
	model deaths = cubic_s1-cubic_s4  /dist=poisson offset=ln_y;
	output out = lin_haz6 pred = pred_deaths;
run;

data lin_haz6;
	set lin_haz6 (in=a)
	haz_grp;
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;
		
run;

title 'two cubic polynomials with knot at 2 years';
title2 'polynomials joined at knot continuous in 2nd derivative'; 
proc sgplot data = lin_haz6;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic;
	refline 2 /axis=x lineattrs=(pattern=dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;
		
*	(j) restricted cubic splines
	generate splines with 5 knots (4 df)
	default has boundary knots at min and max of 'midtime'
	with df-1 internal knots
	default variables created are rcs1 - rcs4;

*	default behaviour of %rcsgen is to orthogonalise the computed spline variables
	here we turn it off, to see the raw splines;
%rcsgen(midtime, df=4, fw=deaths, set = collapsed, orthog = 0);

*	keep the selected knot points in a new macro string;
%let new_knots = &save_knots.;

proc sort data = collapsed out = s nodupkey; by midtime;run;

title 'view spline values for first 12 intervals';
proc print data = s noobs;
	where fu <= 12;
	var fu start end rcs1-rcs4;
run;

proc genmod data = collapsed;
	model deaths = rcs1 /dist=poisson offset=ln_y;
	output out = lin_haz7 pred = pred_deaths;
run;

data lin_haz7;
	set lin_haz7 (in=a)
	haz_grp;
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;
		
run;

title 'plot with just the linear spline term';
title2 ;
proc sgplot data = lin_haz7;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic;
	refline &save_knots. /axis=x lineattrs = (pattern=dash);
	
	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

*	(k) now add remaining spline terms;
%put &=save_knots;
proc genmod data = collapsed;
	model deaths = rcs1-rcs4 /dist=poisson offset=ln_y;
	output out = lin_haz8 pred = pred_deaths;
run;

data lin_haz8;
	set lin_haz8 (in=a)
	haz_grp;
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;
		
run;

title 'with all 4 restricted cubic splines';
proc sgplot data = lin_haz8;
	scatter x = midtime y = haz_grp;
	series x = midtime y = haz_cubic;

*	knot positions from %rcsgen;
	refline &new_knots. /axis=x lineattrs=(pattern=dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;

*	likelihood ratio test all spline terms vs just the linear term;

*	enter LogLikelihood value;
data lrt_pval;
	lrt = 2*abs(3164.9928 -3133.7141	);
	df = 3;							*	in this case;
	p_value = 1 - probchi(LRT,df);
run;

proc print data=lrt_pval;
	title1 "LR test statistic and p-value";
run;


*	(l) look at impact of moving boundary knots to within data, or
	to moving around the internal knots;
*	difficult to get macro parameters passed to a program outside
	of the macro.  Simplest solution is to embed the whole
	process within a macro program
	in this macro, simply supply the knot positions in a single parameter,
	such as:
	%test_knots(0 1 2 3 10)   - alternative internal knots, original boundaries
	%test_knots( 1 2 3) 	- boundary knots within the data
	
	;
	
%macro test_knots(k);
%rcsgen(midtime, fw=deaths, gen=r_new, set = collapsed, knots = &k.);

*	need to know the number of parameters;
%let n_parm = %eval(%sysfunc(countw(&k.," ")) -1);

proc genmod data = collapsed (keep = deaths r_new1 - r_new&n_parm. ln_y midtime risktime);
	model deaths = r_new1 - r_new&n_parm. /dist=poisson offset=ln_y;
	output out = lin_haz_knots&n_parm. pred = pred_deaths;
run;

data lin_haz_knots&n_parm.;
	set lin_haz_knots&n_parm. (in=a)
	haz_grp (keep = midtime haz_grp);
	
	if a then  haz_cubic = 1000* pred_deaths/risktime;	
run;

title 'Evaluate different knot positions';
title2 "knots at &k";
proc sgplot data = lin_haz_knots&n_parm.;
	scatter x = midtime y = haz_grp/ legendlabel = 'interval-specific hazards';
	series x = midtime y = haz_cubic/ legendlabel = 'restricted cubic splines';;
	refline &k., / axis = x lineattrs=(pattern = dash);

	yaxis 
		values = (5 10 20 50 100 150)
		label="Baseline Hazard /1,000 pyrs";
	xaxis label = "Years from diagnosis";
run;
%mend;

%test_knots(1 2 3);
%test_knots(0 1  2  3 10);
%test_knots(0 1  2  3 8);
