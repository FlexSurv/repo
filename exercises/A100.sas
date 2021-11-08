/*
	A100 life table and regression estimates

	model net survival in two time periods with continuous age in model
	comparison with Pohar estimate for the same data

	added June 2020
	
*/

options fmtsearch = (data.melanoma_formats);

*	relative survival using life table methods
	quarterly intervals to one year, anually afterwards;
%rel_surv(infile = data.melanoma, 
	age 		= age, 
	yydx 		= yydx,
	exit 		= surv_mm,
	scale 		= 12,
	list		= pohar,
	censor 		= status(0,4), 
	intervals 	= %str(0 to 1 by .25, 2 to 10 by 1),
	popmort 	= data.popmort, 
	patientid 	= id,
	strat 		= year8594,
	crude		= right = 10);
	
*	append population probabilites for regression modelling of relative survival
	add survival time as time from diagnosis to death/censoring
	adjust censoring variable and event date if death is greater than 10 years
	attained age is constrained to a maximum of 99
	attained year is constrained to a  maximum of 2000;
proc sql;
	create table melanoma_times as select 
		a.*,
		ifn(a.surv_mm/12 < 10.01, a.status, 0) as new_status,
		min(a.surv_mm/12, 10.01) as surv,
		b.rate
		from data.melanoma a left join data.popmort b
		on floor(min(a.age + min(a.surv_mm/12, 10.01),99)) = b._age
		and a.sex = b.sex
		and floor(min(yydx + min(a.surv_mm/12, 10.01), 2000)) = b._year
		order by a.id;
quit;

*	build standard dataset, accept all causes of death as events;
%stset(melanoma_times, new_status(1 2), surv, id, options = noprint);

*	age as continuous with 3 spine variables;
%rcsgen(age, df=3, gen = agercs);

*	fit model with age and time period as non-proportional;
%stpm2(year8594 agercs1 agercs2 agercs3, scale=hazard, df=5, options = noprint,
	tvc= year8594 agercs1 agercs2 agercs3, dftvc=3,  bhazard=rate);
	
*	model predictions using meansurv
	equivalent to an age-standardised analysis where the age distribution
	is different in each time period - as in Pohar life table approach;
	
%range(tt, 0 , 10, 201);
*	first call uses only cases from earlier time period (with confidence interval);
%predict(rs1,  meansurv, ifp = year8594 = 0, timevar = tt, options = ci debug3 );

*	now cases only from the second time period (with confidence interval);
%predict(rs2, meansurv , ifp = year8594 = 1, timevar = tt, options = ci);


*	combine modelled and life table estimates for plotting;
data _events_;
	set _events_ 
	grouped(keep = year8594 right cr_p hi_cr_p  lo_cr_p in = a);

	if a then do;
		if year8594 = 0 then do;
			cr_p1 = cr_p;
			lo_cr_p1 = lo_cr_p;
			hi_cr_p1 = hi_cr_p;
		end;
		if year8594 = 1 then do;
			cr_p2 = cr_p;
			lo_cr_p2 = lo_cr_p;
			hi_cr_p2 = hi_cr_p;
		end;
	end;
run;


* 	Plot the predicted survival functions;
title 'A100:  Survival:  all stages melanoma';
title2 'Net survival and regression estimates from excess hazard model';

proc sgplot data = _events_ noborder;
	where tt ^=. or right ^=.;

	/*	confidence limits		*/
	band upper = rs1_uci lower = rs1_lci x = tt / transparency = .5;
	band upper = rs2_uci lower = rs2_lci x = tt /transparency = .5;

	/*	the modeled estimates	*/
	series y = rs1 x = tt/legendlabel = '1975-84 regression' name = 'a'
		lineattrs = (thickness = 2 color = blue);		
	series y = rs2 x = tt/legendlabel = '1985-94 regression' name = 'b'
		lineattrs = (thickness = 2 color = red);		
	
	/* the life table estimates (these are the Pohar estimates)		*/
	scatter x = right y=cr_p1/ yerrorupper = hi_cr_p1 
		yerrorlower = lo_cr_p1 markerattrs = (color = blue)
		legendlabel = 'life table' name = 'c';	

	scatter x = right y=cr_p2/ yerrorupper = hi_cr_p2 
		yerrorlower = lo_cr_p2 markerattrs = (color = red)
		legendlabel = 'life table' name = 'd';	
	
	xaxis label='time since diagnosis (years)';
	yaxis label='relative survival'
		values = (.6 to 1 by .1);

	keylegend 'b' 'd' 'a' 'c'  / position = topright across = 2 location = inside;
run;
