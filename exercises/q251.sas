/*
==================
 EXERCISE 251
 
 Probability of death in a competing risks framework (flexible parametric models) 
  
 %predict
 %stpm2CM

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

*	(a) load data, assign age group binary variables and assign population
	mortality rates;
data melanoma10;
	set data.melanoma;
	
	if surv_mm > 60.5 then do;
		surv_mm = 60.5;
		if status in(1 2) then status = 0;
	end;
	
	surv_mm = surv_mm/12;	*	stata example has a scale parameter of 12;

	agegrp2 = agegrp = 1;
	agegrp3 = agegrp = 2;
	agegrp4 = agegrp = 3;
	
run;

proc sql;
	create table melanoma_bh as select a.*,
		b.rate
		from melanoma10 a left join data.popmort b
		on floor(a.age + a.surv_mm) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv_mm) = b._year;
quit;

*	(b)  Excess hazard model with age groups as TVC;
%stset(melanoma_bh, status(1 2), surv_mm, id);
%stpm2( agegrp2 agegrp3 agegrp4, df=5, scale=hazard, bhazard=rate,
	tvc =  agegrp2 agegrp3 agegrp4, dftvc=3);
	
*	(c) predict failure functions for each age group and plot;
%range(temptime, 0 , 10, 1001);
%predict(nm1, failure, at=zero, timevar =temptime);

%predict(nm2, failure, at=agegrp2:1 zero, timevar=temptime);

%predict(nm3, failure, at=agegrp3:1 zero, timevar=temptime);

%predict(nm4, failure, at=agegrp4:1 zero, timevar=temptime);

title ;
proc sgplot data = _events_;
	series x = temptime y = nm1/legendlabel='0-44';	
	series x = temptime y = nm2/legendlabel='45-59';	
	series x = temptime y = nm3/legendlabel='60-74';	
	series x = temptime y = nm4/legendlabel='75+';
	yaxis label='Net probability of death';
	xaxis label='Time since diagnosis';
	keylegend /position=topleft location=inside down=4;
run;

*	(d)  estimate crude mortality from excess hazard model and plot ;
%stpm2CM(cm1, at= zero, using = data.popmort, sex =1, diagyear = 1985, diagage = 40, tgen = temptime);

%stpm2CM(cm2, at= agegrp2:1 zero, using = data.popmort,  sex =1, diagyear = 1985,  diagage = 55, tgen = temptime);

%stpm2CM(cm3, at= agegrp3:1 zero, using = data.popmort,  sex =1, diagyear = 1985, diagage = 70, tgen = temptime);

%stpm2CM(cm4, at= agegrp4:1 zero, using = data.popmort,  sex =1, diagyear = 1985, diagage = 80, tgen = temptime);

title ;
proc sgplot data = _events_;
	series x = temptime y = cm1_d/legendlabel='0-44';	
	series x = temptime y = cm2_d/legendlabel='45-59';	
	series x = temptime y = cm3_d/legendlabel='60-74';	
	series x = temptime y = cm4_d/legendlabel='75+';
	yaxis label='crude probability of death (cancer)';
	xaxis label='Time since diagnosis';
	keylegend /position=topleft location=inside down=4;
run;
proc sgplot data = _events_;
	series x = temptime y = cm1_o/legendlabel='0-44';	
	series x = temptime y = cm2_o/legendlabel='45-59';	
	series x = temptime y = cm3_o/legendlabel='60-74';	
	series x = temptime y = cm4_o/legendlabel='75+';
	yaxis label='crude probability of death (other)';
	xaxis label='Time since diagnosis';
	keylegend /position=topleft location=inside down=4;
run;

*	(e) stacked graphs ;
data for_plot;
	set _events_;
	age = 40;
	lo = cm1_d;
	hi = cm1_all;
	output;	
	
	age = 50;
	lo = cm2_d;
	hi = cm2_all;
	output;	
	
	age = 70;
	lo = cm3_d;
	hi = cm3_all;
	output;	
	
	age = 80;
	lo = cm4_d;
	hi = cm4_all;
	output;
	keep age temptime lo hi;
run;

proc sgpanel data = for_plot;
	panelby age;
	band x = temptime lower = 0 upper = lo/legendlabel='P(dead cancer)';
	band x = temptime lower = lo upper = hi/legendlabel='P(dead other causes)';
	band x = temptime lower = hi upper = 1/legendlabel='P(Alive)';
	colaxis label='Time since diagnosis';
	rowaxis label='Crude probability of death';
run;

	
*	(f) splines for age instead of age categories;
%rcsgen(age, gen=rcsage, df=4, orthog=1);
%let age_knots = &save_knots.;
proc datasets lib = work nolist force;
	delete agemat;
	change Tmat = agemat;
quit;
run;

%stpm2(rcsage1 rcsage2 rcsage3 rcsage4, scale = hazard, df=5, bhazard=rate,
	tvc=rcsage1 rcsage2 rcsage3 rcsage4, dftvc=2);

*	macro to assign spline values for a specific age;
data single; age = 1; run;
%macro sp_set (ag, _v, _deg);
*	ag 		= scalar age to generate spline values for
	_v		= variable stub to hold spline values generated
	_deg	= number of variables to generate;
	
	%global &_v.1 &_v.2 &_v.3 &_v.4 &_v.5;

	%rcsgen(age, gen = &_v. ,  knots = &age_knots., tmatrix = agemat, 
		set = single, scalar = &ag., orthog=1);
	proc sql noprint;
		%do _ind = 1 %to &_deg;
			select &_v.&_ind. into :&_v.&_ind. from single ;
		%end;
	quit;
%mend;

%sp_set(40, a, 4);
%stpm2CM(cm1_rcs, at= rcsage1:&a1.  rcsage2:&a2. rcsage3:&a3. rcsage4:&a4.,sex =1, 
	using = data.popmort, diagyear = 1985, diagage = 40, tgen = temptime);
	
%sp_set(55, a, 4);
%stpm2CM(cm2_rcs, at= rcsage1:&a1.  rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., sex =1, 
	using = data.popmort, diagyear = 1985, diagage = 55, tgen = temptime);
	
%sp_set(70, a, 4);
%stpm2CM(cm3_rcs, at= rcsage1:&a1.  rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., sex =1, 
	using = data.popmort, diagyear = 1985, diagage = 70, tgen = temptime);
	
%sp_set(80, a, 4);
%stpm2CM(cm4_rcs, at= rcsage1:&a1.  rcsage2:&a2. rcsage3:&a3. rcsage4:&a4., sex =1, 
	using = data.popmort, diagyear = 1985, diagage = 80, tgen = temptime);
	
title 'dashed lines estimates from model with age continuous';
proc sgplot data = _events_;
	series x = temptime y = cm1_rcs_d/legendlabel='40' lineattrs=(pattern=dash);	
	series x = temptime y = cm2_rcs_d/legendlabel='55' lineattrs=(pattern=dash);	
	series x = temptime y = cm3_rcs_d/legendlabel='70' lineattrs=(pattern=dash);	
	series x = temptime y = cm4_rcs_d/legendlabel='80' lineattrs=(pattern=dash);

	series x = temptime y = cm1_d/legendlabel='0-44' ;	
	series x = temptime y = cm2_d/legendlabel='45-59';	
	series x = temptime y = cm3_d/legendlabel='60-74';	
	series x = temptime y = cm4_d/legendlabel='75+';

	yaxis label='crude probability of death (cancer)';
	xaxis label='Time since diagnosis';
	keylegend /position=topleft location=inside down=4;
run;	
