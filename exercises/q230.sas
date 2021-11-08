/*
==================
 EXERCISE 230
 
 Modelling excess mortality using flexible parametric models I
 
reviewed:  8 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);

*	(a) Load Data, merge expected mortality and fit inital model ;
data melanoma;
	set data.melanoma;
	
	if surv_mm > 120.5 then do;
		surv_mm = 120.5;
		status = 0;
	end;
	
	female = sex = 2;

	surv = surv_mm/12;
	agegrp2 = agegrp = 1;
	agegrp3 = agegrp = 2;
	agegrp4 = agegrp = 3;
run;

proc sql;
	create table melanoma_mm as select a.*,
		b.rate
		from melanoma a left join data.popmort b
		on floor(a.age + a.surv) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv) = b._year
	order by id;
quit;

proc freq data = melanoma_mm;
	table status;
run;

*	fit initial model;
%stset(melanoma_mm, status(1 2), surv, id);
%stpm2( ,df=3, bhazard=rate);

* 	(b) Plot the predicted hazard and survival functions;

%predict(h1, hazard, per=1000, options = ci );
%predict(s1, survival, options = ci );

proc sort data = _events_ out = sorted; by _t_;run;
title 'relative survival plot, melanoma data';
proc sgplot data = sorted;
	band x = _t_ upper = s1_uci lower = s1_lci/legendlabel = '95% CI';
	series y = s1 x = _t_/legendlabel = 'Melanoma';
	xaxis label='time since diagnosis';
	yaxis label='relative survival';
run;

title 'excess hazard plot, melanoma data';
proc sgplot data = sorted;
	band x = _t_ upper = h1_uci lower = h1_lci/legendlabel = '95% CI';
	series y = h1 x = _t_/legendlabel = 'Melanoma';
	xaxis label='time since diagnosis';
	yaxis label='excess mortality rate (per 1,000 py)';
run;


*	(c) Compare fitted values for various df;
%macro comp;
	data fit_comp;run;
	%do _df_ = 2 %to 6 %by 2;
		%stpm2( , df=&_df_., scale=hazard, bhazard=rate );
		%predict(h&_df_., hazard, per=1000);
		%predict(s&_df_., survival);
		data fit_comp;set fit_comp _fit_ (in=a); if a then df = &_df_.;run;
	%end;
%mend;

%comp;


proc sort data = fit_comp; by df descr;run;

*	Plot the excess hazard functions;

proc sort data = _events_ out = sorted; by _t_; run;
title 'excess hazard plot for different number of splines used in fitting';
proc sgplot data = sorted;
	series x = _t_ y = h2/legendlabel = 'df = 2';	
	series x = _t_ y = h4/legendlabel = 'df = 4';	
	series x = _t_ y = h6/legendlabel = 'df = 6';
	xaxis label='time since diagnosis';
	yaxis label='excess mortality rate (per 1,000 py)';
run;

*	Plot the survival functions;
title 'relative survival plot for different number of splines used in fitting';
proc sgplot data = sorted;
	series x = _t_ y = s2/legendlabel = 'df = 2';	
	series x = _t_ y = s4/legendlabel = 'df = 4';	
	series x = _t_ y = s6/legendlabel = 'df = 6';
	xaxis label='time since diagnosis';
	yaxis label='relative survival';
run;

*	combine AIC, BIC and ll values for display;
data fits;set fit_comp; by df descr;

	if first.df then do;
		aic = .;
		bic = .;
		ll = .;
	end;
	retain aic bic ll;
	select (descr);
		when ('AIC (smaller is better)') aic = value;
		when ('BIC (smaller is better)') bic = value;
		when ('-2 Log Likelihood') ll = value/(-2);
		otherwise;
	end;
	
	if last.df;
run;

title 'compare fit statistics for different models';
proc print data = fits noobs label;
	where df ^=.;
	var df aic bic ll;
	label df = 'number of splines'
		aic = 'AIC'
		bic = 'BIC'
		ll = 'log likelihood';
run;

*	(d) Fit a proportional excess hazards model (note eform option);
%stpm2( agegrp2 agegrp3 agegrp4 female year8594, 
	df=3, scale=hazard, bhazard=rate, options = eform);

%predict (h2, hazard,  per=1000 );
%predict (s2, survival);

proc sort data = _events_ out = sorted;by _t_; run;

*	(e)	Plot predicted excess hazard rate;

title 'proportional excess hazard model, males diagnosed 1985-94';
proc sgplot data = sorted;
	where agegrp in(0 3) and female = 0 and year8594=1;
	series x = _t_ y = h2/group=agegrp;
	xaxis label="Time since diagnosis";
	yaxis label="Excess hazard rate";
run;

title 'proportional excess hazard model, males diagnosed 1985-94';
title2 'note log scale for excess hazard';
proc sgplot data = sorted;
	where agegrp in(0 3) and female = 0 and year8594=1;
	series x = _t_ y = h2/group=agegrp;
	
	xaxis label="Time since diagnosis";
	yaxis type=log label="Excess hazard rate (log)";
run;

*	(f) Time-dependent effects for age group;
%stpm2( agegrp2 agegrp3 agegrp4 female year8594, 
	df=3, scale=hazard, bhazard=rate,
	tvc = agegrp2 agegrp3 agegrp4, dftvc = 2 );

%predict (h3, hazard,  per=1000);

%predict (s3, survival);

proc sort data = _events_ out = sorted; by _t_; run;

title 'non-proportional excess hazard model, males diagnosed 1985-94';
proc sgplot data = sorted;
	where agegrp in(0 3) and female = 0 and year8594=1;
	series x = _t_ y = h3/group=agegrp;
	xaxis label="Time since diagnosis";
	yaxis label="Excess hazard rate" ;
run;

proc sgplot data = sorted;
	where agegrp in(0 3) and female = 0 and year8594=1;
	series x = _t_ y = h3/group=agegrp;
	
	xaxis label="Time since diagnosis";
	yaxis type=log label="Excess hazard rate (log)";
run;


*	(g) predicted (time-dependent) excess hazard ratios
		hr denominator all covariates held at reference level;
%range(tt, 0, 10, 401);
%predict (hr2, hrnum=agegrp2: 1 zero, hrdenom =zero , timevar = tt, options = ci);

%predict (hr3, hrnum=agegrp3: 1 zero, hrdenom =zero, timevar = tt, options = ci );

%predict (hr4, hrnum=agegrp4: 1 zero, hrdenom =zero, timevar = tt, options = ci );

proc sort data = _events_ out = sorted; by tt; run;

title 'time-dependent excess hazard ratios, males diagnosed 1975-84';
proc sgplot data = sorted;
	series x = tt y = hr2/legendlabel = 'agegrp2';	
	series x = tt y = hr3/legendlabel = 'agegrp3';	
	series x = tt y = hr4/legendlabel = 'agegrp4';
	xaxis label = 'Time since diagnosis';
	yaxis label = 'Excess mortality rate ratios';
run;

proc sgplot data = sorted;
	band x = tt upper=hr4_uci lower=hr4_lci;
	series x = tt y = hr4/legendlabel = 'agegrp4';

	refline 1;
	xaxis label = 'Time since diagnosis';
	yaxis label = 'Excess mortality rate';
run;

*	(h) difference in relative survival;
%predict( sdiff4, sdiff1= agegrp4:1 female:0 year8594:1 zero,
                sdiff2 = agegrp4:0 female:0 year8594:1 zero, 
                per = 1000, timevar = tt,
                options = ci);
     
proc sort data = _events_ out = sorted; by tt; run;
title 'time-dependent relative survival differences, males diagnosed 1985-94';
proc sgplot data = sorted;
	band x = tt upper = sdiff4_uci lower=sdiff4_lci;
	series x = tt y = sdiff4;
	xaxis label = 'Time since diagnosis';
	yaxis label = 'Difference in Relative Survival';
run;	
       
*	(i) difference in excess mortality rate;					
%predict (hdiff4, hdiff1=agegrp4: 1 female: 0 year8594: 1 zero, timevar = tt,
                hdiff2=agegrp4: 0 female: 0 year8594: 1 zero, per=1000,
                options = ci);
                
proc sort data = _events_ out = sorted; by tt; run;
title 'time-dependent excess hazard differences, males diagnosed 1985-94';
proc sgplot data = sorted;
	band x = tt upper = hdiff4_uci lower=hdiff4_lci;
	series x = tt y = hdiff4;
	xaxis label = 'Time since diagnosis';
	yaxis label = 'Difference in Excess Hazard (per 1000 py)';
run;	
          
