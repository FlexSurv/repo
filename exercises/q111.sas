/*
==================
 EXERCISE 111
 
 poisson regression with interactions
 
 	lifetest
  	genmod
  	sort
  macros
  	%stset
  	%lexis   *****  note use of where clause to exclude _st ^= 1  *****
  	%mrate

reviewed:  4 May 2020
==================
*/

options fmtsearch = (data.melanoma_formats);

data melanoma;
	set data.melanoma (where = (stage = 1));
	years = surv_mm/12;
	death = status = 1;
	ln_yr = log(years);
run;

*	survival and hazard plots;
proc lifetest data = melanoma notable plots= (hazard survival (nocensor));
	time years*status(0,2,4);
	strata year8594;
run;

* Hazard rates prior to Poisson regression;
%mrate(in = melanoma, out = table,  byvar = year8594, per = 1000,
	eventvar = status, eventlist = 1 , timevar = years);


*	introduction of %stset (not required for this analysis);
%stset(melanoma, status(1), years , id ) ;

*	use of %mrate after %stset;
%mrate(byvar = year8594, per = 1000);

	
title ;
*	poisson regression of cancer deaths using genmod with period of diagnosis
	as only parameter to be estimated.  ;
proc genmod data = melanoma;
	model death = year8594/ dist = p offset = ln_yr;
    estimate 'year of diagnosis' year8594 1  /  exp;
run;
	
*	Poisson regression restricted to 10 years follow up. 
	non-cancer deaths are censored;
data melanoma_10;
	set melanoma ;
*	truncate follow-up at 10 years;
	if surv_mm >120 then do;	
		surv_mm = 120;
		years = surv_mm/12;
		ln_yr = log(years);
		if status in(1) then do;
			status = 0;
			death = 0;
		end;
	end;
	zero = 0;
run;

%mrate(in = melanoma_10 , out = table,  byvar = year8594, per = 1000,
	eventvar = status, eventlist = 1 , timevar = years);
	
title ;
proc genmod data = melanoma_10 (where = (surv_mm <= 120));
	model death = year8594/ dist = p offset = ln_yr;
    estimate 'year of diagnosis' year8594 1  /  exp;
run;

*	split follow-up time into 1-year intervals;

*	Split the data to obtain one observation for each life table interval for each individual.
	The scale must be transformed to years;
%lexis (data = melanoma_10 , 
		out = melanoma_split ,
		breaks = 0 to 10 by 1,      
		origin = zero,
		entry = zero,
		exit = surv_mm,
		scale = 12,
		fail = death,
		right = right,
		risk = y,				/*	person-years survived in interval	*/
		lrisk = ln_y,
		lint = length,			/*	length of interval					*/
		nint = fu,
		st = _st);
	
%mrate(in = %str(melanoma_split where  _st = 1), out = table,  byvar = left, per = 1000,
	eventvar = death, eventlist = 1 , timevar = y);

title ;	
proc sgplot data = table;
	scatter x = left y = rate/yerrorlower=cll yerrorupper = clh;
	xaxis values = (0 to 10);
	yaxis values = (0 to 60 by 20);
run;

proc lifetest data = melanoma_10 (where = (_st=1)) 
	plots( only)= hazard notable;
	time years*status(0,2,4);
run;

*	Estimate incidence rate ratios as a function of follow-up;
proc genmod data = melanoma_split (where = (_st=1));
	class left (ref='0');
	model death =  left/dist = poisson offset = ln_y;
	estimate '1' left 1 0 0 0 0 0 0 0 0 -1/exp;
	estimate '2' left 0 1 0 0 0 0 0 0 0 -1/exp;
	estimate '3' left 0 0 1 0 0 0 0 0 0 -1/exp;
	estimate '4' left 0 0 0 1 0 0 0 0 0 -1/exp;
	estimate '5' left 0 0 0 0 1 0 0 0 0 -1/exp;
	estimate '6' left 0 0 0 0 0 1 0 0 0 -1/exp;
	estimate '7' left 0 0 0 0 0 0 1 0 0 -1/exp;
	estimate '8' left 0 0 0 0 0 0 0 1 0 -1/exp;
	estimate '9' left 0 0 0 0 0 0 0 0 1 -1/exp;
run;

* 	Poisson regression adjusting for time since diagnosis;
proc genmod data = melanoma_split (where = (_st=1));
	class left (ref='0') year8594;
	model death =  left year8594/dist = poisson offset = ln_y;
	estimate '1' left 1 0 0 0 0 0 0 0 0 -1/exp;
	estimate '2' left 0 1 0 0 0 0 0 0 0 -1/exp;
	estimate '3' left 0 0 1 0 0 0 0 0 0 -1/exp;
	estimate '4' left 0 0 0 1 0 0 0 0 0 -1/exp;
	estimate '5' left 0 0 0 0 1 0 0 0 0 -1/exp;
	estimate '6' left 0 0 0 0 0 1 0 0 0 -1/exp;
	estimate '7' left 0 0 0 0 0 0 1 0 0 -1/exp;
	estimate '8' left 0 0 0 0 0 0 0 1 0 -1/exp;
	estimate '9' left 0 0 0 0 0 0 0 0 1 -1/exp;
	estimate 'era' year8594 -1 1/exp;
run;

* 	Poisson regression adjusting for age, calendar period and sex;
proc genmod data = melanoma_split (where = (_st=1));
	class left year8594 sex agegrp /ref=first;
	model death =  left year8594 sex agegrp /
		dist = poisson offset = ln_y
		type3 wald;
	estimate '1' left  1 0 0 0 0 0 0 0 0 -1/exp;
	estimate '2' left  0 1 0 0 0 0 0 0 0 -1/exp;
	estimate '3' left  0 0 1 0 0 0 0 0 0 -1/exp;
	estimate '4' left  0 0 0 1 0 0 0 0 0 -1/exp;
	estimate '5' left  0 0 0 0 1 0 0 0 0 -1/exp;
	estimate '6' left  0 0 0 0 0 1 0 0 0 -1/exp;
	estimate '7' left  0 0 0 0 0 0 1 0 0 -1/exp;
	estimate '8' left  0 0 0 0 0 0 0 1 0 -1/exp;
	estimate '9' left  0 0 0 0 0 0 0 0 1 -1/exp;
	estimate 'era' year8594 1 -1/exp;
	estimate 'females' sex -1 1 /exp;
	estimate '45-59' agegrp  1 0 0 -1/exp;
	estimate '60-74' agegrp  0 1 0 -1/exp;
	estimate '75+'  agegrp   0 0 1 -1/exp;
				
run;


* 	Poisson regression ;
proc genmod data = melanoma_split (where = (_st=1));
	class left year8594 sex agegrp /ref=first;
	model death =  left agegrp year8594|sex  /
		dist = poisson offset = ln_y;
/* Effect of females, 1985-94 */
estimate 'females, 1985-94' sex -1 1 year8594*sex -1 1 0 0/exp;
run;
		
* 	Creating dummies;
data melanoma_split;
	set melanoma_split (where = (_st=1));
	
	sex_early=(sex=2)*(year8594=0);
	sex_later=(sex=2)*(year8594=1);
run;

* 	Poisson regression;
proc genmod data = melanoma_split ;
	class left year8594 sex agegrp sex_early sex_later/ref=first;
	model death =  left year8594 sex agegrp sex_early sex_latter/
		dist = poisson offset = ln_y;
estimate 'females, 1985-94' year8594 1 -1 sex -1 1 sex_early -1 1/exp;		
run;

*	poisson regression stratified on calendar period;
proc sort data = melanoma_split;
	by year8594;
run;

proc genmod data = melanoma_split;
	by year8594 ;
	class left year8594 sex agegrp sex_early sex_later/ref=first;
	model death =  left sex agegrp/
		dist = poisson offset = ln_y;	

run;

* 	Poisson regression with only interactions;
proc genmod data = melanoma_split;
	class left year8594 sex agegrp/ref=first;
	model death =  left|year8594 sex|year8594 agegrp|year8594/
		dist = poisson offset = ln_y type3;	
run;
