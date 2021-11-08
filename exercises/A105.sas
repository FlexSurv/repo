/*
	A105   life table and regression estimates of conditional survival

	methods of A100 extended to include conditional estimates and age standardisation

	added June 2020

	age standardised to ICSS(2) standard
	
*/

options fmtsearch = (data.melanoma_formats);

*	(a) relative survival using life table methods
	compute life table estimates based on patients who have survived the first year
	this amounts to a 'late entry' design, with patients only coming into
	the risk period after the first year of followup;
data melanoma;
	set data.melanoma;

	enter_months = 12 ;		*	entry time (in months);

*	create age groups according to ICSS standard for link to ICSS weights below;
	agegroup = input(put(age,z4.), agegrp_o.);
	format agegroup Agegrpn.;
run;

title 'age standardised net survival, conditional on 1 year survival';
%rel_surv(infile = melanoma, 
	age 		= age, 
	yydx 		= yydx,
	exit 		= surv_mm,
	entry		= enter_months,
	scale 		= 12,
	censor 		= status(0,4), 
	intervals 	= %str(0 to 1 by .25, 2 to 10 by 1),
	popmort 	= data.popmort, 
	list 		= pohar,
	stnd		= icss(2),
	weight_lib	= data,
	patientid 	= id,
	strat 		= year8594,

/*	uncomment the following line to view all the life table estimates	*/
/*	age_adj		= 1,*/

	std_estimates = cond_out);		*	save estimates for plotting later;
	
*	(b) append population probabilites for regression modelling of relative survival
	add survival time as time from diagnosis to death/censoring
	adjust censoring variable and event date if death is greater than 10 years
	attained age is constrained to a maximum of 99
	attained year is constrained to a  maximum of 2000;

*	weight individual cases according ratio of ICSS standard and observed proportions; 

*	tabulate age groups by year strata;
proc freq data = melanoma ;
	table agegroup*year8594/out = frout norow nopercent outpct;
run;

*	now we can append population death rates and assign the individual weights for standardisation;
proc sql;
	create table melanoma_times as select 
		a.*,
		ifn(a.surv_mm/12 < 10.05, a.status, 0) as new_status,
		min(a.surv_mm/12, 10.05) as surv,
		b.rate,
		100*c.weight/d.pct_col as weight

		from melanoma a left join data.popmort b
			on floor(min(a.age + min(a.surv_mm/12, 10.05),99)) = b._age
			and a.sex = b.sex
			and floor(min(yydx + min(a.surv_mm/12, 10.05), 2000)) = b._year

		left join data.icss c
			on a.agegroup = c.agegroup
			and c.standard = 2

		left join frout d
			on a.agegroup = d.agegroup
			and a.year8594 = d.year8594

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

*	'meansurv' measure with a 1-year condition ('tcond = 1')
	uses individual weights as above;
%range(tt, 0, 10, 401);
%predict(ms15a,  meansurv,  ifp = year8594 = 0, meansurvwt = weight, timevar = tt, tcond = 1, options = ci);
%predict(ms15b,  meansurv,  ifp = year8594 = 1, meansurvwt = weight, timevar = tt, tcond = 1, options = ci);

*	(e) prepare life table estimates for plotting
	estimates are already sorted by time period and right end of interval;
data conditional;
	set cond_out (keep = year8594 right  as_rel lo_rel hi_rel);
	by year8594 right;

		if year8594 = 0 then do;
			cr_p0 = as_rel;
			lo_cr_p0 = lo_rel;
			hi_cr_p0 = hi_rel;
		end;
		if year8594 = 1 then do;
			cr_p1 = as_rel;
			lo_cr_p1 = lo_rel;
			hi_cr_p1 = hi_rel;
		end;
run;

*	(f) combine estimates for plotting;
data combine;
	set _events_ (where = (tt ^=.)) 
		conditional (where = (right ^=.));
run;

title 'A105:  1 year conditional Survival:  all stages melanoma';
title2 'Age standardised Net survival from life table and excess hazards modelling';

proc sgplot data = combine;

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

