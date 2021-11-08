/*
==================
 EXERCISE 261
 
 cure models 

%stpm2 (with cure option)
%predict (cure, centile options)

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.colon_formats);

*	(a) load the data and fit a cure model with period of diagnosis as a covariate;
data colon;
	set data.colon;
	
	if surv_mm > 120.5 then do;
		surv_mm = 120.5;
		if status in(1 2) then status = 0;
	end;
	id = _n_;
	surv_mm = surv_mm/12;	*	stata example has a scale parameter of 12;	
run;

*	append population values;
proc sql;
	create table colon_bh as select a.*,
		b.rate
		from colon a left join data.popmort b
		on min(floor(a.age + a.surv_mm), 99) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv_mm) = b._year;
quit;

*	create standard dataset;
%stset(colon_bh, status(1 2), surv_mm, id);

*	Cure model using stpm2 (hazard scale is the default);
%stpm2( year8594, df = 6, bhazard = rate, options = cure);

*	(b)	predict cure proportion and median survival in the uncured; 

%predict(cure1, cure);
%predict(median1, survival, uncured = 1, cent = 50);

title 'cure proportion by year of diagnosis';

*	each row has a copy of the cure proportion.  Just select one
	for each covariate pattern;

proc sort data = _events_ out = sorted nodupkey; by year8594;run;
proc print data = sorted noobs label split = '*';
	var year8594 cure1;
	label cure1 = 'cure*proportion*from*PEH model';
run;

title 'median survival in uncured patients, by year of diagnosis';
proc print data = sorted noobs label split = '*';
	var year8594 median1;
	label median1 = 'median*survival';
run;

*	(c) fit cure model with a tvc variable;
%stpm2( year8594, df = 6, bhazard = rate, tvc = year8594, 
	dftvc = 4, options = cure );

%predict(cure2, cure);
%predict(median2, survival, uncured = 1, cent = 50 );

proc sort data = _events_ out = sorted nodupkey; by year8594;run;
title 'cure proportion, by year of diagnosis';
title2 'model with period of diagnosis as TVC variable';
proc print data = sorted noobs label split = '*';
	var year8594 cure2;
	label cure2 = 'cure*proportion*from*PEH model';
run;
title 'median survival of uncured patients';
title2 'model with period of diagnosis as TVC variable';

proc print data = sorted noobs label split = '*';
	var year8594 median2;
	label median2 = 'median*survival';
run;


*	(d) estimate overall survival function and the survival function for the uncured;

*	survival estimation for all subjects;
%predict(all_surv, survival);

*	survival curve for the uncured population;
%predict(uncured, survival, uncured = 1);

*	get median survivals and cure proportions as macro strings;
proc sql noprint;
	select distinct  put(median2, 4.2) into :med separated by " "
		from _events_;

	select distinct put(cure2, 4.2) into :cure separated by " "
		from _events_;
quit;

%let c = %scan(&cure., 1, " ");
%let m = %scan(&med., 1, " ");
%put &=c.;
%put &=m.;

title 'survival plots for all patients vs uncured population';
title2 'diagnosed 1975 - 1984';

proc sort data = _events_ out = sorted; by _t_;run;
proc sgplot data = sorted;
	where year8594 = 0;
	series x = _t_ y = all_surv/legendlabel = 'Survival Overall';
	series x = _t_ y = uncured/legendlabel = 'Survival for uncured';
	xaxis values = (0 to 10 by 2) label = 'Years since diagnosis';
	yaxis label='proportion surviving';

	refline .5;
	refline &m. / 
		labelloc = inside
		axis = x label = "Median in uncured (*ESC*){unicode '000a'x}  = &m.";
	refline &c. / 
		label = "Cure Proportion (*ESC*){unicode '000a'x}    = &c." 
		lineattrs = (pattern = dot thickness = 2);
run;

title 'survival plots for all patients vs uncured population';
title2 'diagnosed 1975 - 1984';
%let c = %scan(&cure.,2," ");
%let m = %scan(&med.,2," ");
proc sgplot data = sorted;
	where year8594 = 1;
	series x = _t_ y = all_surv/legendlabel = 'Survival Overall';
	series x = _t_ y = uncured/legendlabel = 'Survival for uncured';
	xaxis values = (0 to 10 by 2) label = 'Years since diagnosis';
	yaxis label='proportion surviving';

	refline .5;
	refline &m. / 
		labelloc = inside
		axis = x label = "Median in uncured (*ESC*){unicode '000a'x}  = &m. ";
	refline &c./ 
		label = "Cure Proportion (*ESC*){unicode '000a'x}    = &c." 
		lineattrs = (pattern = dot thickness = 2);
run;

 
