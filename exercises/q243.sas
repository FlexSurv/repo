/*==================
 EXERCISE 243		
 REVISED MAY 2015	

%rel_surv, using ICSS weights and internal weights	

reviewed: 7 May 2020
====================

*/

options fmtsearch = (data.melanoma_formats);

proc format ;
	value female
	0 = 'Males'
	1 = 'Females';
	
	value sex
	1 = 'Males'
	2 = 'Females';
	
	invalue icss
	0-44 = 1
	45-54 = 2
	55-64 = 3
	65-74 = 4
	75-99 = 5;
	
	value icsscat
	0-44 = '<45'
	45-54 = '45-54'
	55-64 = '55-64'
	65-74 = '65-74'
	75-99 = '75+';
	
	value icsslab
	1 = '<45'
	2 = '45-54'
	3 = '55-64'
	4 = '65-74'
	5 = '75+';
run;


* (	a) Load melanoma , retain stage 1 cases;
* 	generate an age group variable for the 5 groupings;
data melanoma;
	set data.melanoma (where = (stage=1));

	agegrpICSS = input(age,icss.);

	label 	agegrpICSS = 'Age groups for ICSS';
run;

*	Generate the internal weights based on the age distribution of the data;
proc freq data = melanoma noprint;
	table agegrpICSS/ out = iweights ;
run;

*	append iweights to case data;
proc sql;
	create table melanoma_w as select 
		a.*,
		b.percent/100 as standwei label  = 'Internal age group weights'
		from melanoma a left join iweights b
		on a.agegrpICSS = b.agegrpICSS;
quit;

* 	Age-standardised using iweights;

*	using finer intervals early in followup because of the hazard transform
	assumption of constant hazard within the interval;
title 'melanoma relative survival using internal weighting';
%rel_surv(infile = melanoma_w, 
	weight_var = standwei, standstrata = agegrpICSS,
	patientid=id,
	age = age, 
	exit = surv_mm,
	yydx = yydx,
	censor = status(0 4),
	scale = 12,
	popmort = data.popmort,
	intervals = %str(0 to 1 by 1/12, 1 to 3 by .25, 4 to 10 by 1),
	crude = (mod(right,1) = 0),
	age_adj = (mod(right,1) = 0),
	std_estimates = internal);
	
/* (b) */ 

/* We use ICSS 2 weights for melanoma 

15-44 year: 28%, 45-54 years: 17%, 
55-64years: 21%, 65-74 years: 20%,
75+ years: 14%

*/

*	Generate a variable with the external weights;
data melanoma_w;
	set melanoma_w;
	
	select (agegrpICSS);
		when (1) icss2wei = .28;
		when (2) icss2wei = .17;
		when (3) icss2wei = .21;
		when (4) icss2wei = .20;
		when (5) icss2wei = .14;
		otherwise;
	end;		
run;	

* 	Age-standardised using external weights implemented with iweights;
title 'melanoma relative survival using ICSS (2) weighting';
%rel_surv(infile = melanoma_w, 
	weight_var = icss2wei, standstrata = agegrpICSS,
	patientid=id,
	age = age, 
	exit = surv_mm,
	yydx = yydx,
	censor = status(0 4),
	scale = 12,
	popmort = data.popmort,
	intervals = %str(0 to 1 by 1/12, 1 to 3 by .25, 4 to 10 by 1),
	crude = 1,
	std_estimates = external1); 
	
*	or using built-in ICSS weight table and associated age formats;

%include "&fpsaus./macros/C-SPAN Formats.sas";
title 'melanoma relative survival using ICSS (2) weighting (from reference library)';
%rel_surv(infile = melanoma_w, 
	patientid=id,
	age = age, 
	weight_lib = data,
	stnd = icss(2),		/*	ICSS weights for melanoma	*/
	exit = surv_mm,
	yydx = yydx,
	censor = status(0 4),
	scale = 12,
	popmort = data.popmort,
	intervals = %str(0 to 1 by 1/12, 1 to 3 by .25, 4 to 10 by 1),
	crude = (mod(right,1) = 0),
	age_adj = (mod(right,1) = 0),
	std_estimates = external2); 	
	
*	(c)	compare internal and ICSS weights for each age group;
proc sort data = melanoma_w out = comp nodupkey;
	by agegrpicss;
run;

title 'weights assigned';
proc print data = comp noobs;
	var agegrpicss standwei icss2wei;
	format agegrpicss icsslab.;
run;

* 	It is also possible to compare the standardised estimates ;
data comp_est;
	set internal (in=a) external1 (in=b) external2 (in=c);
	by interval;

	length source $16;
	if a then source = 'internal';
	else if b then source = 'ICSS (2)';
	else if c then source = 'ICSS lib (2)';
	if right  = '5.00';
run;

title 'age standardised estimates with different weighting schemes';
proc print data =  comp_est noobs;
	var source right as_rel ci_rel;
run;

*	(d) Let's now use ICSS 1 weights for melanoma; 
title 'melanoma relative survival using ICSS (1) weighting (from reference library)';
%rel_surv(infile = melanoma_w, 
	age = age, 
	patientid=id,
	weight_lib = data,
	stnd = icss(1),
	exit = surv_mm,
	yydx = yydx,
	censor = status(0 4),
	scale = 12,
	popmort = data.popmort,
	intervals = %str(0 to 1 by 1/12, 1 to 3 by .25, 4 to 10 by 1),
	crude = (mod(right,1) = 0),
	age_adj = (mod(right,1) = 0),
	std_estimates = external3); 

data comp_est2;
	set internal (in=a) external1 (in=b) 
	external2 (in=c) external3 (in=d);
	by interval;

	length source $16;
	if a then source = 'internal';
	else if b then source = 'ICSS (2)';
	else if c then source = 'ICSS lib (2)';
	else if d then source = 'ICSS lib (1)';
	if right  = '5.00';
run;

title 'age standardised estimates with different weighting schemes';
proc print data =  comp_est2 noobs;
	var source right as_rel ci_rel;
run;
	
