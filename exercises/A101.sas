/*
	A101 life table and regression estimates of conditional survival

	methods of A100 extended to include conditional estimates

	added June 2020
	
*/

options fmtsearch = (data.melanoma_formats);

*	(a) conditional relative survival using life table methods
	quarterly intervals to one year, anually afterwards

 	compute life table estimates based on patients who have survived the first year
	this is equivalent to a 'late entry' analysis.  All patients enter the risk
	interval at 1 year after their date of diagnosis.
	the first interval starts at year 1;

data cond;
	set data.melanoma;
	enter = dx+365.24;
run;

*	compute new estimates. all results are saved in the default dataset 'grouped';
%rel_surv(infile = cond, 
	age 		= age, 
	origin		= dx,
	entry		= enter,
	exit 		= exit,
	censor 		= status(0,4), 
	intervals 	= %str(0 to 10 by 1),
	popmort 	= data.popmort, 
	list 		= pohar,
	patientid 	= id,
	strat 		= year8594	);
	
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

*	(c) build standard dataset, accept all causes of death as events;
%stset(melanoma_times, new_status(1 2), surv, id, options = noprint);

*	age as continuous with 3 spine variables;
%rcsgen(age, df=3, gen = agercs);

*	fit model with age and time period as non-proportional;
%stpm2(year8594 agercs1 agercs2 agercs3, scale=hazard, df=5, options = noprint,
	tvc= year8594 agercs1 agercs2 agercs3, dftvc=3,  bhazard=rate);


*	(d) conditional relative survival comparing life table and regression estimates ;

*	'meansurv' measure with a 1-year condition ('tcond = 1');
%range(tt, 0 , 10, 201);
%predict(ms15a,  meansurv,  ifp = year8594 = 0, tcond = 1, timevar = tt, options = ci);
%predict(ms15b,  meansurv,  ifp = year8594 = 1, tcond = 1, timevar = tt, options = ci);


*	(e) prepare life table estimates for plotting
	estimates are already sorted by time period and right end of interval;
data conditional;
	set grouped (keep = year8594 right  cr_p lo_cr_p hi_cr_p);
	by year8594 right;

		if year8594 = 0 then do;
			cr_p0 = cr_p;
			lo_cr_p0 = lo_cr_p;
			hi_cr_p0 = hi_cr_p;
		end;
		if year8594 = 1 then do;
			cr_p1 = cr_p;
			lo_cr_p1 = lo_cr_p;
			hi_cr_p1 = hi_cr_p;
		end;
run;

*	combine estimates for plotting;
data _events_;
	set _events_ conditional;
run;

title 'A101:  1 year conditional survival:  all stages melanoma';
title2 'Net survival and regression estimates from excess hazard model';

proc sgplot data = _events_;
	where tt ^ = . or right ^=.;
	
	band x = tt upper = ms15a_uci lower = ms15a_lci/fillattrs = (transparency = .5);
	band x = tt upper = ms15b_uci lower = ms15b_lci/fillattrs = (transparency = .5);

	series x = tt y = ms15a/legendlabel = '1975-84 regression' name = 'a'
		lineattrs = (thickness = 2 color = blue);
	series x = tt y = ms15b/legendlabel = '1985-94 regression' name = 'b'
		lineattrs = (thickness = 2 color = red);

	scatter x = right y=cr_p0/ yerrorupper = hi_cr_p0 
		yerrorlower = lo_cr_p0 markerattrs = (color = blue)
		legendlabel = 'life table' name = 'c';	

	scatter x = right y=cr_p1/ yerrorupper = hi_cr_p1 
		yerrorlower = lo_cr_p1 markerattrs = (color = red)
		legendlabel = 'life table' name = 'd';	
	
	keylegend 'b' 'd' 'a' 'c'  / position = topright across = 2 location = inside;

	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);
run;

