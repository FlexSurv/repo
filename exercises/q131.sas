/*
==================
 EXERCISE 131

 Model cause-specific mortality using flexible parametric models

	lifetest
 	%stpm2 %predict

reviewed:  4 May 2020
(some random knot locations create inestimable situations)
==================
*/

options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;

*	truncate follow-up time to 120.5 months, scale to years;
data melanoma10;
	set data.melanoma (where = (stage = 1));

	d = status in( 1);

	if surv_mm >120.5 then do;
		surv_mm = 120.5;
		if d = 1 then d = 0;
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

	surv_mm = surv_mm/12;	*	stata example has a scale parameter of 12;
	female = sex = 2;
run;


proc freq data = melanoma10;
	table d ;
run;

proc means data = melanoma10 max min ;var surv_mm;run;

*	(a) Kaplan-Meier curve;
proc lifetest data = melanoma10
	notable
	plots(only)= survival(nocensor)
	outsurv=km_surv;;
	time surv_mm*d(0);
run;

*	get smothed version of K-M hazard function
	default output is _plt_;
%smooth(data=km_surv, time=surv_mm,  survival=survival, width=1, option=noplot);


*	(b) Weibull model (using stpm2);
*	first, set up the standard datset for what follows;
%stset(melanoma10, d(1), surv_mm, id);

*	flexible parametric model with 1 df is equivalent to a Weibull model;
%stpm2( scale=hazard, df=1);		*	no covariates (null model);

*	predict survival curve from this model;
%predict(s1, survival);

*	combine K-M and Weibull survival estimates;
*	select single observation to plot;
proc sort data = _events_  out = surv; by _t_; run;

data plots;
	set km_surv surv;
run;

title 'comparison of non-parametric survival curve with Weibull fit';
proc sgplot data = plots;
	series y = survival x = surv_mm/legendlabel = 'K-M estimate';
	series y = s1 x = _t_ / legendlabel = 'Weibull fit';
	yaxis values = (.6 to 1 by .2);
	xaxis label='Analysis time';
run;

*	(c) Plot hazard function
	predict hazard function;
%predict(h1, hazard);

*	combine with smoothed K-M hazard function;
proc sort data = _events_ out = hz; by _t_;run;
data ploth;
	set _plt_ hz;
run;

title 'comparison of non-parametric hazard curve with hazard from Weibull fit';
proc sgplot data = ploth;
	series y = lambda x = s/ legendlabel = 'Smoothed K-M hazard function';
	series y = h1 x = _t_/legendlabel = 'Weibull fit' ; 
	xaxis label='Analysis time';
run;

*	(d) try stpm2 with 4 df;
%stpm2(scale=hazard, df=4);

%predict(s4, survival);

proc sort data = _events_ out = surv; by _t_; run;
data plots; set km_surv surv; run;

title 'comparison of survival estimates';
proc sgplot data = plots;
	series y = survival x = surv_mm/ legendlabel = 'K-M survival';
	series y = s4 x = _t_/ legendlabel = 'Flexible parametric (4df)';
	yaxis label = 'Survival proportion';
	xaxis label='Analysis time';
run;

%predict(h4, hazard);
proc sort data = _events_ out = hz; by _t_; run;
data ploth; set _plt_ hz; run;

title 'comparison of hazard function estimates';
proc sgplot data = ploth;
	series y = lambda x = s/ legendlabel = 'Smoothed K-M hazard function';
	series y = h4 x = _t_/ legendlabel = 'Flexible parametric (4df)';
	yaxis label = 'hazard (mortality) rate';
	xaxis label='Analysis time';
run;

*	(e)	fit Cox model with one covariate;
proc phreg data =  melanoma10;
	model surv_mm*d(0) = year8594 / 
		risklimits;
run;

*	(f) Flexible parametric model with one covariate and 4 df for splines,
		note option to report exponentiated estimates as in Cox model;
%stpm2(year8594, scale=hazard, df=4, options = eform);

*	(g) predicted values;
*	predict survival at the covariate value for each subject;
%predict(s1ph, survival);

proc sort data = _events_
	(keep = _t_ s1ph  year8594)
	out = s4_plot;
	by _t_;
run;

title 'fitted survival curves from proportional hazards model';
proc sgplot data = s4_plot;
	series y = s1ph x = _t_/group=year8594;
	yaxis label='proportion surviving';
	xaxis label = 'Time since diagnosis';
run;

*	similarly for the hazard functions.  scale to per 1,000 person-years;
%predict(h1ph, hazard, per=1000);

proc sort data = _events_
	(keep = _t_ h1ph  year8594)
	out = h4_plot;
	by _t_;
run;

title 'fitted hazard functions from proportional hazards model';
proc sgplot data = h4_plot;
	series y = h1ph x = _t_/group=year8594;
	yaxis label='Cause-specific mortality rate (per 1,000 pyrs)';
	xaxis label = 'Time since diagnosis';
run;

title 'fitted hazard functions from proportional hazards model';
proc sgplot data = h4_plot;
	series y = h1ph x = _t_/group=year8594;
	yaxis label='Cause-specific mortality rate (per 1,000 pyrs)'
		type = log;
	xaxis label = 'Time since diagnosis';
run;

*	(h) hazard on log scale;
title 'fitted hazard functions from proportional hazards model';
proc sgplot data = h4_plot;
	series y = h1ph x = _t_/group=year8594;
	yaxis label='Cause-specific mortality rate (per 1,000 pyrs)'
		type = log;
	xaxis label = 'Time since diagnosis';
run;


* (i) sensitivity to knots;
*	build macro to cycle through 6 models;
%macro t_knots;
	data est; run;
	data fit;run;
	%do k = 1 %to 6;
		%stpm2(year8594, scale=hazard, df=&k., options = eform);
		%predict(s&k., survival, at=year8594:0);
		%predict(h&k., hazard, at=year8594:0, per = 1000);

		data est;
			set est _parms_ (in=a
				keep = estimate standarderror parameter
				where = (parameter = 'year8594'));
			if a then knots = &k.;
		run;

		data fit;
			set fit _fit_ (in=a
				keep = value  descr
				where = (substr(descr,1,4) in( 'AIC ', 'BIC ')));
			if a then knots = &k.;

		run;
	%end;
%mend;

*	invoke the macro and print table of estimates, SE,
	AIC, BIC;
%t_knots;

*	assemble the fit statistics and print;
data fit2;
	set fit;
	by knots;

	if first.knots then do;
		desc1 = descr;
		val1 = value;
	end;
	retain desc1 val1;

	else output;
run;
data fit3;
	merge fit2 est;
	by knots;
run;

title 'fit statistics from 6 different flexible parametric fits';
proc print data = fit3 noobs label;
	where knots ^=.;
	var knots estimate standarderror val1 value;
	label knots = '# of knots'
		estimate = 'Estimate'
		standarderror = 'SE'
		val1 = 'AIC'
		value = 'BIC';
run;


*	(j) compare hazard and survival estimates;
proc sort data = _events_;
	by _t_;
run;

title 'survival estimates from df = 1 - 6';
proc sgplot data=_events_;
	series x =_t_ y = s1/ legendlabel ='1 knot';
	series x =_t_ y = s2/ legendlabel ='2 knots';
	series x =_t_ y = s3/ legendlabel ='3 knots';
	series x =_t_ y = s4/ legendlabel ='4 knots';
	series x =_t_ y = s5/ legendlabel ='5 knots';
	series x =_t_ y = s6/ legendlabel ='6 knots';
	yaxis label = 'survival proportion';
	xaxis label = 'Time since diagnosis (years)';
run;

proc sort data = _events_;
	by _t_;
run;

title 'hazard estimates from df = 1 - 6';
proc sgplot data=_events_;
	series x =_t_ y = h1/ legendlabel ='1 knot';
	series x =_t_ y = h2/ legendlabel ='2 knots';
	series x =_t_ y = h3/ legendlabel ='3 knots';
	series x =_t_ y = h4/ legendlabel ='4 knots';
	series x =_t_ y = h5/ legendlabel ='5 knots';
	series x =_t_ y = h6/ legendlabel ='6 knots';
	yaxis label = 'mortality rate / 1,000 p-y';
	xaxis label = 'Time since diagnosis (years)';
run;

* (k) knot locations;
/* run the following code to fit 10 models with 5df (6 knots) where
	 the 4 internal knots are selected as random centiles of the
	 distribution of event times.

	As there are many ties in this data we add a small random number to the survival times
	 (otherwise we risk having knots in the same location)
*/

*	re-sort the dataset so predict will work;
proc sort data = _events_; by _study_id_;run;
data _events_;
	set _events_;
	_t_ = _t_+rand("uniform")*0.001;
run;

%macro rand_knots;

	data est; length knots $24;call streaminit(123456);knots = '';run;
	data fit;length knots $24;knots = '';run;

	%do irand = 1 %to 10;
		data _events_;set _events_; sp&irand. = .; hp&irand. = .;run;
		%let plist= ;
		%do j = 1 %to 4;
			%let next = %sysevalf(100*%sysfunc(rand(uniform),6.3));
			%let plist = &plist. &next.;
		%end;
		%let plist = &plist;
		%put random centiles:  &plist.;
		%stpm2(year8594, scale=hazard, knots=&plist., options = noprint,
			knscale=centile);
		%if %sysfunc(exist(_fit_)) ne 0 %then %do;
			%predict(sp&irand.,survival, at=year8594:0);
			%predict(hp&irand.,hazard, at=year8594:0);

			data est;
				set est _parms_ (in=a
				keep = estimate standarderror parameter
				where = (parameter = 'year8594'));
				if a then knots = "&plist.";
			run;

			data fit; set fit _fit_ (in=a
				keep = value  descr
				where = (substr(descr,1,4) in( 'AIC ', 'BIC ')));
				if a then knots = "&plist.";
			run;

			proc datasets lib=work nolist force;
				delete _fit_ ;
				quit;
			run;
		%end;
	%end;
%mend;

%rand_knots;

title 'effect of knot positions on estimate';
proc print data = est noobs split = '*';
	var  estimate standarderror knots;
	label knots ='random knot*positions*(centiles)'
		estimate = 'log(HR)*estimate'
		standarderror = 'SE of*estimate';
run;

proc sort data = _events_ out = surv; by _t_; run;
proc sgplot data = surv;
	series x = _t_ y = sp1/legendlabel = '1';
	series x = _t_ y = sp2/legendlabel = '2';
	series x = _t_ y = sp3/legendlabel = '3';
	series x = _t_ y = sp4/legendlabel = '4';
	series x = _t_ y = sp5/legendlabel = '5';
	series x = _t_ y = sp6/legendlabel = '6';
	series x = _t_ y = sp7/legendlabel = '7';
	series x = _t_ y = sp8/legendlabel = '8';
	series x = _t_ y = sp9/legendlabel = '9';
	series x = _t_ y = sp10/legendlabel = '10';
run;

proc sgplot data = surv;
	series x = _t_ y = hp1/legendlabel = '1';
	series x = _t_ y = hp2/legendlabel = '2';
	series x = _t_ y = hp3/legendlabel = '3';
	series x = _t_ y = hp4/legendlabel = '4';
	series x = _t_ y = hp5/legendlabel = '5';
	series x = _t_ y = hp6/legendlabel = '6';
	series x = _t_ y = hp7/legendlabel = '7';
	series x = _t_ y = hp8/legendlabel = '8';
	series x = _t_ y = hp9/legendlabel = '9';
	series x = _t_ y = hp10/legendlabel = '10';
run;


* (l) Include effect of age group and sex;
proc phreg data=_events_;
	class female year8594 agegrp/ref=first;
	model _t_*d(0) = female year8594 agegrp;
	ods output parameterestimates = cox_est;
run;

%stpm2 (female year8594 a1 a2 a3,  df=4, scale=hazard, options = eform);

*	(m)  compare to Cox model (drop the values related to only the flexible model);
data compare;
	merge cox_est 
		_parms_ (rename = (
				estimate = estimate_flex
				standarderror = stderr_flex)
		where = (substr(parameter,1,3) not in('con' 'rcs')));
run;

title 'comparison of parameter estimates';
title2 'Results from Cox PH model and corresponding flexible parametric model';
title3 'constant and spline parameters not printed';
proc print data = compare noobs split = '*' label;
	var parameter  estimate stderr estimate_flex stderr_flex;
	label estimate 		= 'Cox*estimate' 
		stderr 			= 'Cox SE'
		estimate_flex 	= 'flex*estimate'
		stderr_flex		= 'flex SE';
run;


* (n) obtaining predictions;
* predict using at option;

*	define a new time variable (201 values equally spaced from 0 to 10);
%range(temptime, 0, 10, 201);

*	get survival prediction with all covariates held at reference (zero) level;
%predict (S0, survival, at=zero, timevar=temptime);

title 'survival of males, diagnosed 1975-84, youngest age group';

*	plot functions require data to be sorted by x axis variable;
proc sort data = _events_ out = surv; by temptime;
proc sgplot data = surv;
	series y = s0 x = temptime;
	yaxis values = (0.6 to 1 by .2)  label = 'Survival proportion';
	xaxis label = 'Time since diagnosis (years)';
run;

*	predict survival for specific values of covariates
	request confidence limits;
%predict (S_F_8594_age75, survival, options = ci,
	at=female:1 year8594:1 agegrp:3 zero, timevar=temptime);

title 'survival of females, diagnosed 1985-94, oldest age group';
proc sort data = _events_ out = surv; by temptime;
proc sgplot data = surv;
	band x = temptime lower = 	S_F_8594_age75_lci 
		upper = S_F_8594_age75_uci
		/legendlabel = '95% CI' fillattrs = (transparency = .6);

	series x = temptime y = S_F_8594_age75 /
		legendlabel = 'ages 75+';

	yaxis values = (0.6 to 1 by .2)  label = 'Survival proportion';
	xaxis label = 'Time since diagnosis (years)';
run;
