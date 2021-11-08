/*
==================
 EXERCISE 112
 
 poisson regression, change timescale to attained age
 	genmod
 	
 	%hazard_late
 	%smooth
 	%lexis    ***  note use of where clause to exclude _st ^= 1   ****

reviewed: 4 May 2020
====================
*/

options fmtsearch = (data.diet_formats);

/*	version of dataset with age at diagnosis, age at death as time scale	*/
data diet1;
	set data.diet;
	
	energy = energy*1000;
	dead = (dox-dob)/365.24;	*	time scale is attained age (death-birth);
	entry = (doe-dob)/365.24;	*	entry is start of at risk period (age at diagnosis);
	time_at_risk = dead - entry;
	bmi=weight/(height/100*height/100);
	
	label energy = 'Total Energy (kcals/day)';
run;

*	smoothed hazard plots (overall and by hieng) ;

title 'Time scale is age at diagnosis';
%hazard_late(data=diet1,  entry= entry, exit = dead, fail =chd,
	 censor = 0);
	
%hazard_late(data=diet1, strat = hieng, entry= entry, exit = dead, fail =chd,
	 censor = 0);

	
* 	Timescale: Time-in-study;
data diet2;
	set data.diet;
	
	energy = energy*1000;
	dead = (dox-doe)/365.24;	*	time scale is years since diagnosis;
	entry = 0;
	time_at_risk = dead - entry;
	ln_time = log(time_at_risk);
	bmi=weight/(height/100*height/100);
	
	label energy = 'Total Energy (kcals/day)';
run;

* Hazard function (overall and by hieng);
title 'Time scale is years since diagnosis';
%hazard_late(data=diet2,  entry= entry, exit = dead, fail =chd,
	 censor = 0);
	 
%hazard_late(data=diet2, strat = hieng, entry= entry, exit = dead, fail =chd,
	 censor = 0);
	
* 	Modelling the rate, no adjustment for timescale;
proc genmod data = diet2;
	class hieng;
	model chd=hieng/dist=poisson offset=ln_time;
	estimate "hienergy" hieng -1 1/exp;
run;

*	Adjustment for confounders job and bmi;
proc genmod data = diet2;
	class hieng;
	model chd=hieng job bmi/dist=poisson offset=ln_time;
	estimate "hienergy" hieng -1 1/exp;
	estimate "job" job  1/ exp;
	estimate "bmi" bmi 1/exp;
run;

*	model the rate, adjust for timescale age;


* 	Adjustment for confounders job and bmi;
*	time scale here is age;
%lexis (data = diet1, 
		out = diet1_split ,
		breaks = %str(30 ,50 ,60, 72),      
		origin = dob,
		entry = doe,
		exit = dox,
		scale = 365.24,
		fail = chd,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,				/*	person-years survived in interval	*/
		lrisk = ln_y,
		nint = ageband,
		st = _st);
		
proc print data = diet1_split noobs ;
	where id <=10 and _st = 1;
	var id _t0 _t ageband y chd ;
run;


proc freq data = diet1_split;
	where _st = 1;
	table ageband*chd/missing norow nocol nopercent;
run;


* 	Poisson regression adjusted for ageband;

proc genmod data = diet1_split;
	where _st=1;
	class hieng ageband;
	model chd=hieng ageband/dist=poisson offset=ln_y;
	estimate "hienergy" hieng -1 1/exp;
	estimate "age 50-59" ageband -1 1 0/exp;
	estimate "age 60-72" ageband -1 0 1/exp;
run;
	
* 	Adjustment for confounders job, bmi;

proc genmod data = diet1_split;
	where _st=1;
	class hieng ageband job;
	model chd=hieng ageband job bmi/dist=poisson offset=ln_y;
	estimate "hienergy" hieng -1 1/exp;
	estimate "age 50-59" ageband -1 1 0/exp;
	estimate "age 60-72" ageband -1 0 1/exp;
	estimate "job 1" job  -1 1 0/ exp; 
	estimate "job 2" job  -1 0 1/ exp; 
	estimate "bmi" bmi 1/exp;
run;
	

* 	Modelling the rate, adjusting for timescale time-in-study;

data diet3;
	set data.diet;
	
	energy = energy*1000;
	dead = dox/365.24;	*	time scale is attained age;
	entry = doe/365.24;
	time_at_risk = dead - entry;
	ln_time = log(time_at_risk);
	bmi=weight/(height/100*height/100);
	
	label energy = 'Total Energy (kcals/day)';
run;


/* Split follow up time */

%lexis (data = diet3 , 
		out = diet3_split ,
		breaks = %str(0,5,10,15,22),      
		origin = doe,
		entry = doe,
		exit = dox,
		scale = 365.24,
		fail = chd,
		cens = 0,
		right = _t,
		left = _t0,
		risk = y,				/*	person-years survived in interval	*/
		lrisk = ln_y,
		nint = fuband,
		st = _st);

proc print data = diet3_split noobs ;
	where id <=10 and _st = 1;
	var id _t0 _t fuband y chd ;
run;	
	
proc freq data = diet3_split;
	where _st=1;
	table fuband*chd/missing nocol norow nopercent;
run;

title 'Time scale is time since diagnosis';
proc genmod data = diet3_split;
	where _st=1;
	class hieng fuband;
	model chd=hieng fuband/dist=poisson offset=ln_y;
	estimate "hienergy" hieng -1 1/exp;
	estimate "fuband 5" fuband -1 1 0 0/exp;
	estimate "fuband 10" fuband -1 0 1 0/exp;
	estimate "fuband 10" fuband -1 0 0 1/exp;
run;


* Adjustment for confounders job and bmi;
proc genmod data = diet3_split;
	where _st=1;
	class hieng fuband job;
	model chd=hieng fuband job bmi/dist=poisson offset=ln_y;
	estimate "hienergy" hieng -1 1/exp;
	estimate "fuband 5" fuband -1 1 0 0/exp;
	estimate "fuband 10" fuband -1 0 1 0/exp;
	estimate "fuband 10" fuband -1 0 0 1/exp;
	estimate "job 1" job  -1 1 0/ exp; 
	estimate "job 2" job  -1 0 1/ exp; 
	estimate "bmi" bmi 1/exp;
run;

