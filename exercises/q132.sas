/*
==================
 EXERCISE 132
 
 Extending Q131  Time-dependent effects in flexible parametric models
 %stpm2(); %predict(); in-line macros

reviewed:  4 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;

*	truncate follow-up time to 60.5 months, scale to years;
data melanoma5;
	set data.melanoma  (where = (stage = 1));;
		
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
	
	surv_mm = surv_mm/12;	*	convert time variable to years survival;
	female = sex = 2;		*	0 = male, 1 = female;
run;

*	(a) fit a Cox model 
	Assees the PH assumption for age group;
proc phreg data = melanoma5;
	class agegrp year8594 female /ref=first;
	model surv_mm*ca(0) = year8594 female agegrp;
	hazardratio year8594;
	hazardratio female;
	hazardratio agegrp;
	output out = ph_test ressch = year8594 female agegrp_1 agegrp_2 agegrp_3;
run;

title 'schoenfeld residual plots for age groups';
title2 'ages 45-59';
proc sgplot data = ph_test;
	where ca=1;
	loess x = surv_mm y = agegrp_1;
	refline 0 ;
run;

title2 'ages 60-74';
proc sgplot data = ph_test;
	where ca=1;
	loess x = surv_mm y = agegrp_2;
	refline 0 ;
run;

title2 'ages 75+';
proc sgplot data = ph_test;
	where ca=1;
	loess x = surv_mm y = agegrp_3;
	refline 0 ;
run;

*	(b) fit flexible parametric model;
*	create standard file for %stpm2 and %predict;
%stset( melanoma5, ca(1), surv_mm, id );

*	fit model with main effects of age, sex and period of diagnosis.  request HR
	cumulative log hazard scale is the default scale;
%stpm2(female year8594 a1 a2 a3, df=4, options = eform) ;

*  	define alternate time variable to speed up prediction and reduce then
	number of points the plot requires;
%range(temptime, 0, 5, 501);

%predict(h0_ag1, hazard, at=zero, per=1000, timevar = temptime);

%predict(h0_ag2, hazard, at=a1:1 zero, per=1000, timevar = temptime);

%predict(h0_ag3, hazard, at=a2:1 zero, per=1000, timevar = temptime);

%predict(h0_ag4, hazard, at=a3:1 zero, per=1000, timevar = temptime);

title 'Age specific hazard plots';

proc sgplot data = _events_;
	series x = temptime y = h0_ag1/lineattrs = (color = red)
		legendlabel = '< 45' name = 'a';

	series x = temptime y = h0_ag2/lineattrs = (color = blue)
		legendlabel = '45 - 59' name = 'b';

	series x = temptime y = h0_ag3/lineattrs = (color = magenta)	
		legendlabel = '60 - 74' name = 'c';

	series x = temptime y = h0_ag4/lineattrs = (color = green)
		legendlabel = '75 +' name = 'd';
	xaxis label="Time since diagnosis (years)";
	yaxis label="Cause specific mortality rate (per 1000 py's)" ;

	keylegend 'd' 'c'  'b' 'a'  / location = inside position = topright down = 4;
run; 

title 'Age specific hazard plots (log hazard scale)';
proc sgplot data = _events_;
	series x = temptime y = h0_ag1/lineattrs = (color = red)
		legendlabel = '< 45' name = 'a';

	series x = temptime y = h0_ag2/lineattrs = (color = blue)
		legendlabel = '45 - 59' name = 'b';

	series x = temptime y = h0_ag3/lineattrs = (color = magenta)	
		legendlabel = '60 - 74' name = 'c';

	series x = temptime y = h0_ag4/lineattrs = (color = green)
		legendlabel = '75 +' name = 'd';
	xaxis label="Time since diagnosis (years)";
	yaxis label="Cause specific mortality rate (per 1000 py's) log scale" type = log;

	keylegend 'd' 'c'  'b' 'a'  / location = inside position = topright down = 4;
run; 

*	(c) add Time-dependent effects for age group.  suppress printed output;
%stpm2(female year8594 a1 a2 a3, df=4, options = noprint,
	tvc = a1 a2 a3, dftvc=2) ;
	
*	lrt is difference in log-likelihood between this model and the previous one
	log-likelihood is read from the _fit_ data file in the work directory;
data lrt_pval;
	lrt = abs(5017.2-4997.1);
	df = 6;							*	6 new spline terms;
	p_value = 1 - probchi(LRT,df);
run;

title1 "Likelihood Ratio test statistic (p-value) for incusion of TVC variables";

proc print data=lrt_pval noobs split = '*';
	var lrt df p_value;
	label lrt = 'Likelihood*ratio test*(chi-square)'
		df = 'chi-square*degrees*of freedom'
		p_value = 'p-value';
run;


*	(d) predict the hazard for each age group (with all other covariates
		held at the reference level (males, diagnosed 1975-84);
%predict(h1_ag1, hazard, at=zero, per=1000, timevar = temptime);

%predict(h1_ag2, hazard, at=a1:1 zero,  per=1000, timevar = temptime);

%predict(h1_ag3, hazard, at=a2:1 zero,  per=1000, timevar = temptime);

%predict(h1_ag4, hazard, at=a3:1 zero,  per=1000, timevar = temptime);

title 'Age specific hazard plots from main effects and from TVC model';
title2 'dashed lines are from time-dependent model';

proc sgplot data = _events_;
	series x = temptime y = h0_ag1/lineattrs = (color = red)
		legendlabel = '< 45' name = 'a';	
	series x = temptime y = h1_ag1/lineattrs = (color = red pattern=dash);	

	series x = temptime y = h0_ag2/lineattrs = (color = blue)
		legendlabel = '45 - 59' name = 'b';	
	series x = temptime y = h1_ag2/lineattrs = (color = blue pattern=dash);	

	series x = temptime y = h0_ag3/lineattrs = (color = magenta)	
		legendlabel = '60 - 74' name = 'c';	
	series x = temptime y = h1_ag3/lineattrs = (color = magenta pattern=dash);	

	series x = temptime y = h0_ag4/lineattrs = (color = green)	
		legendlabel = '75 +' name = 'd';	
	series x = temptime y = h1_ag4/lineattrs = (color = green pattern=dash);	

	xaxis label="Time since diagnosis (years)";
	yaxis label="Cause specific mortality rate (per 1000 py's)";

	keylegend 'd' 'c'  'b' 'a'  / location = inside position = topright down = 4;
run; 

*	(e) time-dependent hazard ratios;
%predict(hr2, hrnum = a1:1 zero, timevar = temptime);

%predict(hr3, hrnum = a2:1 zero, timevar = temptime) ;

%predict(hr4, hrnum = a3:1 zero, timevar = temptime, options=ci); 

title 'Hazard ratio plots (versus <45 age group)';
proc sgplot data = _events_;
	series x = temptime y = hr2/legendlabel = '45 - 64';
	series x = temptime y = hr3/legendlabel = '65 - 74';
	series x = temptime y = hr4/legendlabel = '75 +';
	xaxis label = "Time since diagnosis (years)";
	yaxis type = log label = "Hazard ratio";
run;

title 'Hazard ratio plots (oldest versus <45 age group)';
title2 'showing 95% confidence band';
proc sgplot data = _events_;
	band x = temptime upper = hr4_uci lower = hr4_lci/legendlabel = '95% CI';
	series x = temptime y = hr4/legendlabel = '75 +';
	xaxis label = "Time since diagnosis (years)";
	yaxis type = log label = "Hazard ratio";
run;

* (f) Difference in hazard rates;

%predict(hdiff4, hdiff1 = a3: 1 zero, per = 1000, timevar = temptime, options=ci);

title 'Difference in hazard rates (age 75+ vs <45)';
title2 'showing 95% confidence band';
proc sgplot data = _events_;
	band x = temptime upper = hdiff4_uci lower = hdiff4_lci
		/legendlabel = '95% CI';
	series x = temptime y = hdiff4 /legendlabel = 'diff in hazard';
	xaxis label = "Time since diagnosis (years)";
	yaxis label = "Difference in mortality rate (per 1,000pyr)";
run;


*	(g) Difference in survival functions;
%predict(s1, survival, at = female:1 year8594:1 zero, timevar = temptime);

%predict(s2, survival, at = a3:1 female:1 year8594:1  zero, timevar = temptime);

title 'two age-specific survival functions';
title2 'estimates for females diagnosed 1985-94'; 
proc sgplot data = _events_;
	series x = temptime y = s1/legendlabel='<45';
	series x = temptime y = s2/legendlabel='75+';
	xaxis label = "Time since diagnosis (years)";
	yaxis label = "S(t)";
run;
		
*	now predict the difference (with confidence limits);
%predict(sdiff4, sdiff1=a3:1 female :1 year8594:1 zero,
				sdiff2 = female :1 year8594:1 zero, options = ci, timevar = temptime);

title 'difference between two age-specific survival functions';
title2 'estimates for females diagnosed 1985-94'; 
proc sgplot data = _events_;
	band x = temptime upper = sdiff4_uci lower = sdiff4_lci/legendlabel = '95% CI';
	series x = temptime y = sdiff4/legendlabel = '75+ vs <45';
	xaxis label = "Time since diagnosis (years)";
	yaxis label = "Difference in survival functions";
run;
	

* (h) varying df for time-dependent effects;
%macro tvc3;
	%do i_tvc = 1 %to 3;
		%stpm2(female year8594 a1 a2 a3, df=4, scale=hazard,
		tvc=a1 a2 a3, dftvc=&i_tvc.);
			
		%predict(hr_df&i_tvc., hrnum=a3:1 zero, timevar = temptime);		/*	, options = ci);		*/

	%end;
%mend;
	
%tvc3;
title 'hazard ratio plots for varying TVC models of age';
proc sgplot data = _events_;
	series x = temptime y = hr_df1 /legendlabel = 'linear (1 df)';
	series x = temptime y = hr_df2/legendlabel = '2 df';
	series x = temptime y = hr_df3/legendlabel = '3 df';
	yaxis type=log;
run;

* (	i) two time-dependent effects;
%stpm2( female  a1 a2 a3, df=4, scale=hazard,
	tvc = a1 a2 a3 female, dftvc=3);
	
*	HR for females vs males, estimated with all other covariates at reference level;
%predict (hr_f_ag1, hrnum = female:1 zero, options=ci, timevar = temptime);

*	HR for females vs males, 
	estimated for oldest age group, all other covariates at reference level;
%predict(hr_f_ag4, hrnum = female:1 a3:1 zero, hrdenom = a3:1 zero, options=ci, timevar = temptime);

title 'Hazard ratio for females vs males';
title 'Showing effect of level of other TVC covariates';
proc sgplot data = _events_;
	series x = temptime y = hr_f_ag1/lineattrs = (color = red)
	 legendlabel  ='with age at youngest';
	series x = temptime y = hr_f_ag1_lci/lineattrs = (color = red pattern = dash);
	series x = temptime y = hr_f_ag1_uci/lineattrs = (color = red pattern = dash);

	series x = temptime y = hr_f_ag4/lineattrs = (color = blue)
	 legendlabel  ='with age at oldest';
	series x = temptime y = hr_f_ag4_lci/lineattrs = (color = blue pattern = dash);
	series x = temptime y = hr_f_ag4_uci/lineattrs = (color = blue pattern = dash);
	yaxis type = log;
run;

