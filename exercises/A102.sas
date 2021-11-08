/*
	A 102 life table and regression estimates

	how do life table estimates differ from regression estimates
	
1 June 2020

*/

options fmtsearch = (data.melanoma_formats);

*	(a) relative survival using life table methods
	quarterly intervals to one year, anually afterwards;
%rel_surv(infile = data.melanoma, 
	age 		= age, 
	yydx 		= yydx,
	exit 		= surv_mm,
	scale 		= 12,
	censor 		= status(0,4), 
	intervals 	= %str(0 to 1 by .25, 2 to 10 by 1),
	popmort 	= data.popmort, 
	patientid 	= id);
	
*	(b) append population probabilites for regression modelling of relative survival
	add survival time as time from diagnosis to death/censoring
	adjust censoring variable and event date if death is greater than 10 years
	attained age is constrained to a maximum of 99
	attained year is constrained to a  maximum of 2000;
proc sql;
	create table melanoma_times as select 
		a.*,
		ifn(a.surv_mm/12 < 10.05, a.status, 0) as new_status,
		min(a.surv_mm/12, 10.05) as surv,
		b.rate
		from data.melanoma a left join data.popmort b
		on floor(min(a.age + min(a.surv_mm/12, 10.05),99)) = b._age
		and a.sex = b.sex
		and floor(min(yydx + min(a.surv_mm/12, 10.05), 2000)) = b._year
		order by a.id;
quit;

*	(c) build standard dataset, accept all causes of death as events
	age as continuous with 3 spine variables;
%stset(melanoma_times, new_status(1 2), surv, id, options = noprint);
%rcsgen(age, df=3, gen = agercs);

*	save output in excel spreadsheet;

/*ods excel file = "&fpsaus.\reports\A102 life table and regression.xlsx"*/
/*	options (index = 'yes' embedded_titles='yes' sheet_interval='none');*/

/*ods excel options(sheet_label="base");*/

*	(d) fit base model (no covariates);
%stpm2( df=5,  bhazard=rate);

*	predict survival from base model;
%range(tt, 0, 10, 1001);
%predict(surv0,  survival, timevar = tt, options = ci );

*	combine and plot;
proc sort data = _events_  (where = (surv0 ^=.)) out = sorts0 ;by tt;run;
data sorts0;
	set sorts0 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

* 	Plot the predicted survival function;
/*ods excel options(sheet_label="base plot" sheet_interval='now');*/

title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'base model with no covariates';
proc sgplot data = sorts0 noborder;

	/*	confidence limits		*/
	band upper = surv0_uci lower = surv0_lci x = tt / transparency = .5;

	/*	the modeled estimates	*/
	series y = surv0 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	
	/* the life table estimates (these are the Pohar estimates)		*/
	scatter x = right y=cr_p/ legendlabel = 'Life Table' name = 'b'
	yerrorupper = hi_cr_p yerrorlower = lo_cr_p;	
	
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);

	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;
	
*	compare excess hazard plots;
*	smoothed hazard plot on the log scale;
*	note that the width is reduced from the default (1/5 of the observed range, ie  2).  The smoothing
	algorithm should have fine intervals to work with, be sure to use %range to create 1000 intervals;
%smooth (data=sorts0 (where = (surv0 ^=.)), out = smooth0, width = 1,
	time=tt, survival=surv0, option = noplot);

*	predict excess hazard from base model and scale to per 1,000 p-y ;

*	scale excess hazard to 1,000 person-years;
data combine_haz;
	set smooth0 (in=a)
	grouped;

	lambda = lambda*1000;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'base model with no covariates';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;


*	(e) add age into the model;
/*ods excel options(sheet_label="age" sheet_interval='now');*/
%stpm2(agercs1 agercs2 agercs3, df=5,  bhazard=rate);

*	model predictions using meansurv
	equivalent to an age-standardised analysis - as in Pohar life table approach;
	
*	predict survival averaged over all covariate combinations;
%predict(surv1,  meansurv, timevar= tt, options = ci );

proc sort data = _events_ nodupkey out = sorts1;by tt;run;
data sorts1;
	set sorts1 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

/*ods excel options(sheet_label="age plot" sheet_interval='now');*/
* 	Plot the predicted survival function;
title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age';
proc sgplot data = sorts1 noborder;
	band upper = surv1_uci lower = surv1_lci x = tt / transparency = .5;
	series y = surv1 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	scatter x = right y=cr_p/legendlabel = 'life table estimates' name = 'b'
		yerrorupper = hi_cr_p yerrorlower = lo_cr_p;	
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
	keylegend  'a' 'b' / position = topright across = 1 location = inside;
run;

*	compare excess hazard plots;
%smooth (data=sorts1 (where = (surv1 ^=.)), time=tt, width = 1, 
	out = smooth1, survival=surv1, option = noplot);

*	scale hazard to 1,000 person-years;
data combine_haz;
	set smooth1
	grouped;

	LAMBDA = 1000*LAMBDA;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;

*	sex and age into the model;
/*ods excel options(sheet_label="age sex" sheet_interval='now');*/
%stpm2(sex agercs1 agercs2 agercs3, df=5,  bhazard=rate);

%predict(surv2,  meansurv,  timevar = tt, options = ci );

proc sort data = _events_ (where = (surv2 ^=.))  out = sorts2;by tt;run;
data sorts2;
	set sorts2 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

/*ods excel options(sheet_label="age sex plot" sheet_interval='now');*/
* 	Plot the predicted survival function;
title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age and sex';
proc sgplot data = sorts2 noborder;
	band upper = surv2_uci lower = surv2_lci x = tt / transparency = .5;
	series y = surv2 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	scatter x = right y=cr_p/ legendlabel = 'life table estimates' name = 'b'
		yerrorupper = hi_cr_p yerrorlower = lo_cr_p;	
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
	keylegend  'a' 'b' / position = topright across = 1 location = inside;
run;

*	compare excess hazard plots;
%smooth (data=sorts2 (where = (surv2 ^=.)), time=tt, width = 1, out = smooth2, 
	survival=surv2, option = noplot);

*	scale hazard to 1,000 person-years;
data combine_haz;
	set smooth2
	grouped;

	LAMBDA = 1000*LAMBDA;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age and sex';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;

run;

*	add age TVC into the model;

/*ods excel options(sheet_label="age(tvc) sex" sheet_interval='now');*/
%stpm2(sex agercs1 agercs2 agercs3, df=5, tvc = agercs1 agercs2 agercs3, dftvc = 3, bhazard=rate);

%predict(surv3,  meansurv,  timevar = tt, options = ci );

proc sort data = _events_ (where = (surv3 ^=.)) out = sorts3;by tt;run;
data sorts3;
	set sorts3 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

* 	Plot the predicted survival function;
/*ods excel options(sheet_label="age(tvc) sex plot" sheet_interval='now');*/
title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with sex age (age TVC)';
proc sgplot data = sorts3 noborder;
	band upper = surv3_uci lower = surv3_lci x = tt / transparency = .5;
	series y = surv3 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	scatter x = right y=cr_p/ legendlabel = 'life table estimates' name = 'b'
		yerrorupper = hi_cr_p yerrorlower = lo_cr_p;		
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
	keylegend  'a' 'b' / position = topright across = 1 location = inside;
run;

*	compare excess hazard plots;
%smooth (data=sorts3 (where = (surv3 ^=.)), time=tt, width = 1, out = smooth3,
	survival=surv3, option = noplot);

*	scale hazard to 1,000 person-years;
data combine_haz;
	set smooth3
	grouped;

	LAMBDA = 1000*LAMBDA;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with sex age (age TVC)';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;

run;

*	model with sex, age, sex as TVC;
ods excel options(sheet_label="age sex(tvc)" sheet_interval='now');
%stpm2(sex agercs1 agercs2 agercs3, df=5, tvc = sex, dftvc = 3, bhazard=rate);

*	have to use all observations because of meansurv;
%predict(surv4,  meansurv,  timevar = tt, options = ci );

proc sort data = _events_ (where = (surv4 ^=.)) out = sort4;by tt;run;
data sort4;
	set sort4 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

* 	Plot the predicted survival function;
/*ods excel options(sheet_label="age sex(tvc) plot" sheet_interval='now');*/
title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age sex (TVC sex)';
proc sgplot data = sort4 noborder;
	band upper = surv4_uci lower = surv4_lci x = tt / transparency = .5;
	series y = surv4 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	scatter x = right y=cr_p/ legendlabel = 'life table estimates' name = 'b'
		yerrorupper = hi_cr_p yerrorlower = lo_cr_p;		
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;

*	compare excess hazard plots;
%smooth (data=sort4 (where = (surv4 ^=.)), time=tt, width = 1, out = smooth4, 
	survival=surv4, option = noplot);

*	scale hazard to 1,000 person-years;
data combine_haz;
	set smooth4
	grouped;

	LAMBDA = 1000*LAMBDA;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age sex (TVC sex)';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;

run;


*	full model:  age, sex, both as TVC;
/*ods excel options(sheet_label="both tvc" sheet_interval='now');*/
%stpm2(sex agercs1 agercs2 agercs3, df=5, 
	tvc = sex agercs1 agercs2 agercs3, dftvc = 3, bhazard=rate);

%predict(surv5,  meansurv,  timevar = tt, options = ci );

proc sort data = _events_ (where = (surv5 ^= .)) out = sorts5;by tt;run;
data sorts5;
	set sorts5 
	grouped(keep = right cr_p hi_cr_p  lo_cr_p in = a);
run;

* 	Plot the predicted survival function;
/*ods excel options(sheet_label="both tvc plot" sheet_interval='now');*/
title 'A102:  Survival:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age and sex, both as TVC';
proc sgplot data = sorts5 noborder;
	band upper = surv5_uci lower = surv5_lci x = tt / transparency = .5;
	series y = surv5 x = tt/legendlabel = 'Overall Survival' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	scatter x = right y=cr_p/ legendlabel = 'life table estimates' name = 'b'
		yerrorupper = hi_cr_p yerrorlower = lo_cr_p;	
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;

*	compare excess hazard plots;
%smooth (data=sorts5 (where = (surv5 ^=.)), time=tt, width = 1, out = smooth5,
	survival=surv5, option = noplot);

*	scale hazard to 1,000 person-years;
data combine_haz;
	set smooth5
	grouped;

	LAMBDA = 1000*LAMBDA;
	excess = 1000*excess;
run;

title 'A102:  Excess hazard:  all stages melanoma';
title2 'Net survival vs regression estimates from excess hazard model';
title3 'model with age and sex, both as TVC';
proc sgplot data = combine_haz;
	series x = s y = lambda/legendlabel = 'regression model' name = 'a';
	scatter x = right y = excess/legendlabel = 'life table estimates' name = 'b';
	yaxis label ='excess mortality /1,000 p-y (log scale)' type = log;
	xaxis label ='time since diagnosis (years)';
	keylegend  'a' 'b'/ position = topright across = 1 location = inside;
run;
/*ods excel close;*/
