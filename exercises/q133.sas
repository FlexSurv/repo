/*
==================
 EXERCISE 133
 
 Extending Q131 with alternative link functions

reviewed:  4 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;
*	truncate follow-up time to 60.5 months, scale to years;
data melanoma5;
	set data.melanoma (where = (stage = 1));
		
	ca = status in( 1);
	
	if surv_mm >60.5 then do;
		surv_mm = 60.5;
		if ca = 1 then ca = 0;
	end;
	
*	agegroup dummies;
	a0 = 0;
	a1 = 0;
	a2 = 0;
	a3 = 0;

	select (agegrp);
		when (1) a1 = 1;
		when (2) a2 = 1;
		when (3) a3 = 1;
		otherwise;
	end;

	label a0 = '0-44'
		a1 = '45-59'
		a2 = '60-74'
		a3 = '75+';
	
	surv_mm = surv_mm/12;	*	change scale from months to years;
	female = sex = 2;		*	recode female to be 0 = males, 1 = females;

	keep id ca surv_mm female a0-a3 year8594;
run;


%stset( melanoma5, ca(1), surv_mm, id );

*	(a) PH Model;
%stpm2(female a1 a2 a3 year8594, scale = hazard, df = 4, options = eform);
data ph_fit;set _fit_;run;

*	macro for predicting survival and hazard functions for each age group
	first paramter names data file for saving resulting hazard estimates
		file of survival estimates from PH model will be named PH_<s_out>
		file of survival estimates from PO model will be named PO_<s_out>
	first paramter names data file for saving resulting hazard estimates
		similarly for hazard estimates (PH_<h_out>>,  PO_<h_out>);
%macro H_S_byage(h_out, s_out);
*	alternate time variable to speed up processing;
%range(tt, 0, 5, 501);

* 	get scale parameter from _model_ database;
data _temp_;
	set _model_ (where = (comp = 'scale'));
	call symput('scale', desc);
run;

%if &scale. = hazard %then %let s = PH;			*	model is on the proportional hazards scale;
%if &scale. = odds %then %let s = PO;		*	model is on the proportional odds scale;

*	reference age group;
	%predict(s_age0_&s., survival, at=zero, timevar = tt);
	
	%predict(h_age0_&s., hazard, at=zero, timevar = tt);
	
*	other age groups;	
	%do ia = 1 %to 3;
		%predict(s_age&ia._&s., survival, at = a&ia.:1 zero, timevar = tt);
		%predict(h_age&ia._&s., hazard, at = a&ia.:1 zero, timevar = tt);
	%end;
%mend;

%h_s_byage(ph_haz, ph_surv);

*	(b) Proportional Odds Model;
%stpm2(female a1 a2 a3 year8594, scale = odds, df = 4, options = eform);
data po_fit;set _fit_;run;

%h_s_byage(po_haz, po_surv);


*	(c) Compare survival and hazard functions;

title 'survival plot from different link functions';
title2 'solid lines are proportional hazards, dashed lines are proportional odds';
proc sgplot data = _events_;
	series x = tt y = s_age0_ph/legendlabel = '<40 PH';
	series x = tt y = s_age0_po/legendlabel = '<40 PO' lineattrs = (pattern = dash);
	series x = tt y = s_age3_ph/legendlabel = '75+ PH';;
	series x = tt y = s_age3_po/legendlabel = '75+ PO' lineattrs = (pattern = dash);
	keylegend / down=2;
	yaxis label = 'Survival proportion';
	xaxis label = 'Time since diagnosis (years)';
run;


title 'hazard plot';
proc sgplot data = _events_;
	series x = tt y = h_age0_ph/legendlabel = '<40 PH';;
	series x = tt y = h_age0_po/legendlabel = '<40 PO';
	series x = tt y = h_age3_ph/legendlabel = '75+ PH';
	series x = tt y = h_age3_po/legendlabel = '75+ PO';

	keylegend / down=2;
	yaxis label = 'Mortality rate';
	xaxis label = 'Time since diagnosis (years)';
run;

*	(d) Compare AIC and BIC;
data fits; merge ph_fit po_fit (rename=(value = po_value)); by descr;run;
proc print data = fits noobs label;
	var descr value po_value;
	label value = 'value from PH model'
		po_value = 'value from PO model';
run;

*	(e) Hazard ratio for female;
%predict(or_female_age0_7584, hrnum = female: 1 zero, hrdenom = female:0 zero, timevar = tt, options=ci);

title "Odds Ratio for sex (age<45, diagnosed 1975-1984)";
proc sgplot data = _events_;
	band x = tt upper = or_female_age0_7584_uci lower = or_female_age0_7584_lci / legendlabel = '95% CI';
	series x = tt y = or_female_age0_7584/legendlabel = 'Odds ratio Females vs Males';
	xaxis label = 'Years since diagnosis';
	yaxis label = 'Odds ratio ';
run;

*	(f) Compare ods ratios for different covariate patterns;		
%predict(or_female_age3_7584, hrnum = female:1 a3:1 zero, hrdenom = female:0 a3:1 zero, timevar = tt); 

title 'odds ratios for sex, diagnosed 1975-84';
title2 'different age groups';
proc sgplot data = _events_;
	series x = tt y = or_female_age0_7584/legendlabel = '<45';
	series x = tt y = or_female_age3_7584/legendlabel = '75+';

	xaxis label = 'Years since diagnosis';
	yaxis label = 'Odds ratio ' values = (.5 to .7 by .05);
run;

*	(g) Fit Aranda-Ordaz link function;
%stpm2(female a1 a2 a3 year8594, scale = theta, df = 4);
data fits; merge fits _fit_ (rename=(value = ao_value)); by descr;run;

title 'fit statistics from differing link functions';
proc print data = fits noobs label;
	var descr value po_value ao_value;
	label value = 'value from PH model'
		po_value = 'value from PO model'
		ao_value = 'value from AO model';
run;

*	(h) Show estimate of theta with 95% CI;
data theta_est; 
	set _parms_;  
	where parameter = 'ln_theta'; 

	theta = exp(estimate);
	se = theta*standarderror;		*	delta method gives variance of exp(ln_theta);
	lower = exp(lower);
	upper = exp(upper);
run;
title 'theta value from A-O model';
proc print data	= theta_est noobs;
	var alpha theta se tvalue probt lower upper;
run;
