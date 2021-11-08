/*
	A 104 standardise regression estimates with meansurvwt

	17 June 2020

	based on 
	http://pauldickman.com/software/stata/model-based-standardisation/
	by Paul Dickman, March 2019


*/

options fmtsearch = (data.colon_formats);

*	(a) create age groups to the ICSS standard using C-SPAN formats and create dummy variables
	for age groups to be used in regression model;
data colon;
	set data.colon (where = (stage=1));

	surv_mm = surv_mm/12;
	id = _n_;

	agegrp2 = 45 <= age < 55;
	agegrp3 = 55 <= age < 65;
	agegrp4 = 65 <= age < 75;
	agegrp5 = age >= 75;

	agegroup = input(put(age,z4.), agegrp_o.);
	format agegroup agegrpn.;

	keep surv_mm status yydx age id agegroup agegrp2 agegrp3 agegrp4 agegrp5;
run;

*	(b) tabulate age groups and save proportions (saved as percents)
	and assign individual-level weights;
proc freq data= colon;
	table agegroup/out= freq;
run;

*	create individual weights, based on ;
proc sql;
	create table colon_w as select
		a.*,
		c.weight/b.percent*100 as w			/*	individual-level weights		*/
		from colon a left join freq b
			on a.agegroup = b.agegroup

		left join data.icss c
			on a.agegroup = c.agegroup
			and standard = 2;
quit;

*	display the individual-level weights that have been created;
proc means data = colon_w mean min max;
	class agegroup;
	var w;
run;

*	(c) fit regression model and compute weighted survival estimates;
%stset(colon_w, status(1 2), surv_mm, id);

%rcsgen(yydx, df=3, gen=yearspl);

%stpm2 (yearspl1 yearspl2 yearspl3 agegrp2 agegrp3 agegrp4 agegrp5, scale= hazard, df = 4, options = eform,
  tvc= yearspl1 yearspl2 yearspl3 agegrp2 agegrp3 agegrp4 agegrp5,  dftvc= 1);

%range(temptime, 0, 10, 101);

*  (population-averaged) survival;
%predict( s_unweighted, meansurv, timevar= temptime);

*  (standardised) survival;
%predict( s_weighted, meansurv, meansurvwt = w, timevar= temptime);

*	(d) compare unweighted and weighted regression estimtates;
title 'A104: Colon cancer age standardised survival curves';
proc sgplot data = _events_;
	where temptime ^=.;
	series x = temptime y = s_unweighted/
		legendlabel = "survival (unstandardised)" lineattrs = (pattern = solid) name = '2';

	series x = temptime y = s_weighted/
		legendlabel = "survival (standardised)" lineattrs = (pattern = dot) name = '1';

	yaxis values = (0 to 1 by .2) label = "All-cause survival";
	xaxis label = "Years since diagnosis";

	keylegend '1' '2' / position = bottomleft location = inside across = 1;
run;


*	(e) Now let's just estimate 5-year survival to get insight into what meansurv is doing;

data _events_; set _events_; t5=5; run;

* 5-year survival conditional on covariates;
%predict(s5, survival, timevar = t5);

* 5-year  (population-averaged) survival;
%predict(s5_unweighted, meansurv, timevar = t5);

* 5-year  (population-averaged) survival standardised to ICSS;
%predict(s5_weighted, meansurv, meansurvwt = w, timevar = t5);

proc sort data = _events_ out = sorted nodupkey; by agegroup yydx ; run;

title 'sample of estimates at 5 years'; 
proc print data = sorted (obs = 5) label;
	var yydx agegroup s5 s5_unweighted s5_weighted;
	label agegroup 	= 'age group'
		s5			= 'conditional on covariates'
		s5_unweighted = 'population-weighted'
		s5_weighted	  = 'ICSS standardised';
run;

data _events_;
	set _events_;
	s5w = s5*w;

	label s5 = 'Unweighted estimates'
		s5w	 = 'Weighted estimate';

run;

title 'Comparison of survival estimates at 5 years';
proc means data = _events_ mean min max maxdec = 4;
	var s5 s5w;
run;


