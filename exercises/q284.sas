/*
==================
 EXERCISE 284	
 REVISED MAY 2015	

	estimate years life lost

%rcsgen
%stpm2
%predict

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);


*	(a) load data and create append expected rates;
proc format ;
	value female
	0 = 'Males'
	1 = 'Females';
	
	value sex
	1 = 'Males'
	2 = 'Females';
run;

data melanoma;
	set data.melanoma;
	
	patid = _n_;
	
	if surv_mm > 120.5 then do;
		surv_mm = 120.5;
		status = 0;
	end;
	
	surv = surv_mm/12;
		
	fem = sex=2;
	format fem female.;
run;

*	merge on expected rates at exittime;
proc sql;
	create table melanoma_bh as select a.*,
		b.rate
		from melanoma a left join data.popmort b
		on min(int(a.age + a.surv),99) = b._age
		and a.sex = b.sex
		and int(a.yydx + a.surv) = b._year
	order by a.patid;
quit;

*	create standard dataset for analysis;
%stset(melanoma_bh, status(1 2), surv, patid);

*	(b) Fit a flexible parametric model including year, age and sex;
%rcsgen (age, df = 4, gen = sag, orthog = 1);		*	spline variables for age;
%rcsgen (yydx, df = 4, gen = syr, orthog = 1);		*	spline variables for year;


%stpm2 (sag1 sag2 sag3 sag4 syr1 syr2 syr3 syr4 fem, 
	scale = hazard, df = 5, bhazard =rate, options = orthog,
		tvc = sag1 sag2 sag3 sag4 syr1 syr2 syr3 syr4 fem, dftvc = 3);	
		
*	(c)  Predict loss in expectation of life for each subject;
%predict (ll,lifelost, mergeby = _age  sex _year, diagage = age, diagyear =yydx,
	 nodes = 20, tinf = 80, using = data.popmort, maxyear = 2000, options = ci);


*	(d) Create a graph that shows how the loss in expectation of life 
	varies over age, for males diagnosed in 1994.;

proc sort data = _events_ (where = (sex = 1 and yydx = 1994))
	out = for_plot;
	by age;
run;

title 'loss in expectation of life by age at diagnosis';
title2 'males, diagnosed 1994';
proc sgplot data = for_plot;
	series x = age y = ll;
	xaxis label="Age at diagnosis";
	yaxis values=(0 to 45 by 5) label="Future Years of Life Lost"; 
run;

*	(e) List life expectancy and loss in expectation of life for selected ages, 
	as well as total number of life years lost, for dx 1994;
proc sort data = _events_ 
	(where = (age in(50 60 70 80) and sex in(1 2) and yydx = 1994))
	nodupkeys
	out = list1;
	by age sex yydx;
run;

title 'life expectancy and loss in life expectancy for selected ages';
proc print data = list1 noobs;
	var age sex yydx survexp survobs ll;
run;
	
proc means data = _events_ sum;
	where yydx = 1994;
  	var ll ;
run;

*	(f) Predict loss in expectation of life if males had the same 
	relative survival as females;
data _events_;
	set _events_;
	
	fem = 1;
run;

%predict (ll_alt, lifelost, mergeby = _age sex _year, diagage = age,
	 diagyear = yydx, nodes = 40, tinf = 80, using = data.popmort, 
	  maxyear = 2000); 

*	(g) Calculate number of life years that could potentially be saved 
	if males diagnosed in 1994 had the same relative survival as females;
	
data _events_;
	set _events_;
	
	ll_diff = ll - ll_alt;
run;

title 'total life years lost in males diagnosed in 1994 due to male survival deficit';
proc means data = _events_ sum n;
	where yydx = 1994;
  	var ll_diff ;
run;

proc sort data = _events_ 
	(where = (age in(50 60 70 80) and sex = 1 and yydx = 1994))
	nodupkeys
	out = list2;
	by age sex yydx;
run;

title 'life years lost by a male patient due to male survival deficit';
title2 'by age at diagnosis';
proc print data = list2 noobs;
	by sex yydx;
	id sex yydx;
	var age  ll ll_alt ll_diff;
run;
	
