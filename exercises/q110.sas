/*
==================
 EXERCISE 110
 procs in use
 	genmod (poisson regression)
 	univariate
 	freq 
 	sgplot
 	print
 	%mrate

reviewed:  4 May 2020
==================
*/

options fmtsearch = (data.diet_formats);

/*	scale survival time to years	*/
data diet;
	set data.diet;
	
	years = (dox-doe)/365.25;
	ln_yr = log(years);
	
	label energy = 'Total Energy (kcals/day)';
run;


*	incidence rate by age hieng;
%mrate(in = diet, out = table,  byvar = hieng, per = 1000,
	eventvar = chd, eventlist = 1 , timevar = years);
	
title ;
*	poisson regression (use estimate statement for Incidence Rate Ratio (exp of parameter estimate);
proc genmod data = diet;
	model chd = hieng/ dist = poisson offset = ln_yr;
    estimate 'Beta hieng' hieng 1 -1/ exp;
run;
        

*	create energy categorical variables;
proc univariate data = diet noprint;
	histogram energy /normal;
run;

*	create new variable with values of lower cut points
	these are close to the 25th and 75th percentile of the energy distribution;
data diet;
	set diet;
	if energy < 2500 then eng3 = 1500;
	else if 2500 <= energy < 3000 then  eng3 = 2500;
	else eng3 = 3000;
run;

proc freq data = diet;
	table eng3;
run;

*	incidence rate by age energy categorical variable;
%mrate(in = diet, out = table,  byvar = eng3, per = 1000,
	eventvar = chd, eventlist = 1 , timevar = years);
	
proc sgplot data = table;
	scatter x= eng3 y = rate/YERRORLOWER=cll YERRORupper=clh;
run;

*	create your own dummy variablesan 
	use boolen test to set flags;
data diet;
	set diet;
	x1 = energy < 2500;			*	1 if energy <2500, otherwise 0;
	x2 = 2500 <= energy < 3000;	*	1 if energy <= 2500, < 3000, otherwise 0;
	x3 = energy >= 3000;		*	1 if energey >= 3000, otherwise 0;
run;
	
*	poisson regression with dummy variables;
proc genmod data = diet;
	model chd = x2 x3/ dist = poisson offset = ln_yr;
    estimate 'Beta x2' x2 1 / exp;
    estimate 'Beta x3' x3 1 / exp;
run;

proc genmod data = diet;
	model chd = x1 x3/ dist = poisson offset = ln_yr;
    estimate 'Beta x1' x1 1 / exp;
    estimate 'Beta x3' x3 1 / exp;
run;

*	poisson regression with class statement;
proc genmod data = diet;
	class eng3 (ref='1500');
	model chd = eng3/ dist = poisson offset = ln_yr;
    estimate 'eng3 2500' eng3 1 0 -1 /  exp;	*	2500 vs reference;
   	estimate 'eng3 3000' eng3 0 1 -1 /  exp;	*	3000 vs reference;

run;

*	compute rate manually;
proc means data = diet noprint;
	output out = means mean( years chd )=;
run;

data rate;
	set means;
	deaths = chd*_freq_;
	py = years*_freq_;
	rate = deaths/py;
run;

proc print data = rate noobs;
	var deaths py rate;
run;

*	mortality rate with no stratification variable, to check on above calculation;
%mrate(in = diet, out = table, per = 1000,
	eventvar = chd, eventlist = 1 , timevar = years);
