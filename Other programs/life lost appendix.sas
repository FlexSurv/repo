/*
	life lost Appendix
	
	29 Apr 2023
	
	Appendix A
	years life lost example using melanoma data

*/

*	folder structure
	rsrc
		- data		location of datasets used in exercises
		- exercises	exercises from Albuquerque workshop
		- examples	some additional sas material
		- macros	sas macros for relative survival and flexible parametric regression;
	
/*%let rsrc = <path to local copy of the GitHub survival repository>;*/
%let rsrc = S:\S & E Unit - System files\system files\Survival resources;

*	macros required for fitting and prediction;
%include "&rsrc.\programs\survival macros\regression methods.sas";
%include "&rsrc.\programs\survival macros\formats.sas";

libname data "&rsrc.\data";

options fmtsearch = (data.melanoma_formats);

proc format ;
	value sex
	1 = 'Males'
	2 = 'Females';
run;

*	load data, restrict to patients 50+;
data melanoma;
	set data.melanoma (where = (age >= 50));
			
	if surv_mm > 120.5 then do;		*	truncate followup at 10 years;
		surv_mm = 120.5;
		status = 0;
	end;

	surv = surv_mm/12;		*	scale survival time to years;
	format sex sex.;
run;

*	append population mortality rate at exit time;
proc sql;
	create table melanoma_bh as select a.*,
		b.rate
		from melanoma a left join data.popmort b
		on min(int(a.age + a.surv),99) = b._age
		and a.sex = b.sex
		and min(int(a.yydx + a.surv), 2000) = b._year
	order by a.id;
quit;

*	create analytic dataset;
%stset(melanoma_bh, status(1 2), surv, id);

*	spline variables for age and year of diagnosis;
%rcsgen (age, df = 3, gen = a);			*	spline variables for age;
%rcsgen (yydx, df = 3, gen = y);		*	spline variables for year of diagnosis;

*	Fit a flexible parametric model including year, age and sex
	allow for non-proportionality in year;
%stpm2 (a1 a2 a3 y1 y2 y3  sex, scale = hazard, df = 4, bhazard =rate,
	tvc =  y1 y2 y3, dftvc = 2);	
		
*	Predict ll, loss in expectation of life for each subject with confidence limits
	also create these variables for each patient in the dataset
	-------		------------------
	survobs		expected life years for this cancer patient
	survobs_lci	lower 95% confidence limit of survobs
	survobs_uci	upper 95% confidence limit of survobs

	survexp		expected life years for population member with same age, sex
				identified in the year of diagnosis of the patient;

%predict (ll, lifelost, mergeby = _age sex _year, diagage = age, diagyear = yydx,
	 nodes = 20, tinf = 80, using = data.popmort, maxyear = 2000, options = ci);


*	trend in expectation of life for selected ages;
proc sort data = _events_ (where = (age in (55 75)))
	out = for_plot1 nodupkey;
	by yydx age sex;
run;

title 'trends in expectation of life for selected ages';
title2 'melanoma';
proc sgpanel data = for_plot1;
	panelby sex age/ novarname ; 

	band x = yydx upper = survobs_uci lower = survobs_lci/
		legendlabel = '95% confidence band';

	series x = yydx y = survexp/ legendlabel = 'general population'
		lineattrs = (pattern = dash color = red);

	series x = yydx y = survobs/ legendlabel = 'melanoma cases';

	colaxis label="Year of diagnosis";
	rowaxis values=(0 to 40 by 10) label="Expectation of life (years)"; 
run;


*	expected life years lost for a patient diagnosed in 1994,
	by age at diagnosis;
proc sort data = _events_ (where = (yydx = 1994))
	out = for_plot2 nodupkey;
	by  age sex;
run;

title 'life expectation years lost, by age';
title2 'melanoma diagnosed in 1994';
proc sgpanel data = for_plot2;
	panelby sex / novarname ; 

	series x = age y = ll/ legendlabel = 'melanoma cases';

	colaxis label="age at diagnosis";
	rowaxis values=(0 to 10 by 2) label="loss in expectation of life (years)"; 
run;
