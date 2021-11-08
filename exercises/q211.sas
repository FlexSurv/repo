/*
==================
 EXERCISE 211
 
 Modelling excess mortality using Poisson regression 
 
 %rel_surv
 %rcsgen
 genmod with user link
 sgplot
 iml

reviewed:  17May 2020
==================

*/

*	load the melanoma data, truncate followup at 10 years;
data melanoma_mm;
	set data.melanoma;
	
	if surv_mm > 120.5 then do;
		surv_mm = 120.5;
		status = 0;
	end;

	zero = 0;
run;

*	(a) fine split (monthly to 6 years);
%rel_surv(infile = melanoma_mm ,
	yydx = yydx, patientid=id,
	age = age, sex=sex, origin=zero, exit = surv_mm, scale = 12,
	censor = status(0,4), intervals = 0 to 6 by 0.0833,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat=sex agegrp year8594);
	
*	(b) Create variables and fit PEH spline model;
data grouped;
	set grouped;
		
	female = sex = 2;
	midtime = (right+left)/2;
run;

%rcsgen(midtime, knots=.05 0.25 0.75 1.5 2.5 4,  set = grouped);

proc genmod data = grouped ;
	class agegrp female year8594 /ref = first;
	
	fwdlink link = log(_MEAN_-d_star);		
	invlink ilink = exp(_XBETA_)+d_star;	

	model d = rcs1 - rcs5 agegrp female year8594 / 
		dist=poisson offset = ln_y;

	estimate "female" female 1  -1/exp;
	estimate "year8594" year8594 1  -1/exp;
	estimate "age 45-59" agegrp 1 0 0 -1/exp;
	estimate "age 60-74" agegrp 0 1 0 -1/exp;
	estimate "age 75 +" agegrp  0 0 1 -1/exp;

	output out = obs xbeta = pred;
run;


*	(c) Predict excess hazard rate and plot for oldest and youngest groups;
data obs;
	set obs;
	
	haz = exp(pred-ln_y);
run;

title 'Excess hazard plot';
proc sgplot data = obs;
	where agegrp in(0 3) and female = 0 and year8594 = 1;
	series x = midtime y = haz/group = agegrp;
	xaxis label = 'Time since diagnosis (years)';
	yaxis label = 'Excess Hazard';
run;

title 'Excess hazard plot (log scale)';
proc sgplot data = obs;
	where agegrp in(0 3) and female = 0 and year8594 = 1;
	series x = midtime y = haz/group = agegrp;
	xaxis label = 'Time since diagnosis (years)';
	yaxis type = log label = 'Excess Hazard';
	inset "The lines are parallel as the model assumes proportional excess hazards";

run;

*	(d) Time dependent effects;
*	Interactions between spline variables and age group;
*	defining new variables for interactions makes the iml code
	easier later;

*	use some macro code to simplify the data step code
	note that the result of this macro is returned as 
	text to be insterted into the data step where it is called
	this is allowed because the macro itself only contains macro
	statements and fuctions;
%macro inter(lev1, lev2, gen);
%let code = ;
%do i = 1 %to &lev1.;
	%do j = 1 %to &lev2.;
		%let code = &code. 
			%str(&gen.&i.rcs&j. = (agegrp = (&i.-1))*rcs&j.;);
	%end;
%end;
 		&code.
%mend inter;

	
data grouped_int;
	set grouped;
	%inter(4, 5, age);	*	generate 20 age-rcs variables from age1rcs1 to age4rcs5;

*	prototype of code that is generated from the above call;
*	age<i>rcs<j> = (agegrp = <i>-1)*(rcs<j>);
*	where <i> runs from 1 - 4 and <j> runs from 1 - 5;
*	agegroup has values 0 - 3, so <i>-1 runs from 0 - 3;
*	(agegrp = <i>-1) has a value of 1 if the record is for agegroup i-1 and is zero otherwise;
	
run;

proc genmod data = grouped_int ;
	class agegrp  /ref = first;
	
	fwdlink link = log(_MEAN_-d_star);		
	invlink ilink = exp(_XBETA_)+d_star;	

	model d = rcs1-rcs5 agegrp female year8594
	age2rcs1-age2rcs5 age3rcs1-age3rcs5 age4rcs1-age4rcs5/ covb
		dist=poisson offset = ln_y;

	ods output parameterestimates = M_sp_peh;
	ods output covb = cov_peh;
	output out = rcs_obs xbeta = pred stdxbeta = stderr;
run;

*	(e) Obtain excess hazard ratios for each age group with 95% CI 
	with youngest age as the reference category;
proc sql;
	create table ehratio as select
		a.agegrp, a.female, a.year8594, a.midtime,
		exp((a.pred - a.ln_y) - (x.pred - x.ln_y)) as ehr,
		exp((a.pred - a.ln_y) - (x.pred - x.ln_y))*(a.stderr+x.stderr) as stderr
		from rcs_obs a left join rcs_obs x
			on a.female = x.female
			and a.year8594 = x.year8594
			and a.fu = x.fu
			
			where x.agegrp = 0
		order by a.midtime;
		
		alter table ehratio add lo_ehr num add hi_ehr num;
		
		update ehratio set lo_ehr = ehr-1.96*stderr ,
			hi_ehr = ehr+1.96*stderr;
quit;
		
title 'Excess hazard ratio plot (females, diagnosed 1985-94)';
proc sgplot data = ehratio;
	where  female = 1 and year8594 = 1 and agegrp >0;
	series x = midtime y = ehr/group = agegrp;
	xaxis label = 'Time since diagnosis (years)';
	yaxis type = log label='Excess Hazard Ratio';

run;

* (f) Plot excess hazard ratio for oldest group with 95% CI;
title 'Excess hazard ratio plot (females, diagnosed 1985-94)';
title2 'EHR for oldest vs youngest age group';
proc sgplot data = ehratio;
	where  female = 1 and year8594 = 1 and agegrp = 3;
	band x = midtime lower=lo_ehr upper=hi_ehr /legendlabel='95% conidence interval';
	series x = midtime y = ehr;
	yaxis type = log label='Excess Hazard Ratio';
	refline 2.9 /label='EHR from base model' labelloc=inside;

run;

*	(g) Estimated Relative Survival Curves with 95% CIs ;

*	generate spline variables for 200 equally spaced intervals;
*	that knots specified earlier are available in local macro string;
data rs_est;
	do midtime = 0.025 to 6 by .025;
		output;
	end;
run;

%rcsgen(midtime, knots = &save_knots., gen = rcs, set = rs_est); 

*	iml to compute relative survival from estimated cumulative excess hazard;
*	the matrix calculations are not complicated, but it is tricky to get
	the right parameters and corresponding variances and covariances for
	a specific covariate pattern;
proc iml;
*	parameter estimates;
	use m_sp_peh;
	read all var {estimate} into parm ;
*	print parm;

*	variance/covariance matrix of estimated parameters;	
	use cov_peh;
	read all var _num_ into cov;
*	print cov;

*	need a lower triangluar matrix for summing cumulative hazard;
	cummat = j(240,240,0);
	do i = 1 to 240;
		do j = 1 to i;
			cummat[i,j] = 0.025;
		end;
	end;

*	parameter estimates and covariance matrix for youngest age group (males):
	intercept rcs1-rcs5 year8495 ;
	b1 = parm[{1 2 3 4 5 6 12}];
	v1 = cov[{1 2 3 4 5 6 11},{1 2 3 4 5 6 11}];
	
*	now for the oldest age group:	
	intercept rcs1-rcs5 agegrp4 year8495 age4rcs1 - age4rcs5;
	b4 = parm[{1 2 3 4 5 6 9 12 23 24 25 26 27}];
	v4 = cov[{1 2 3 4 5 6 9 11 22 23 24 25 26},{1 2 3 4 5 6 9 11 22 23 24 25 26}];
	
*	time points for estimation (as spline variables);
	use rs_est;
	read all into time;

*	covariate pattern for youngest age group (with splines);	
	des1 = j(nrow(time),1,1)||time[,{2 3 4 5 6}]||J(nrow(time),1,1);
	des4 = j(nrow(time),1,1)||time[,{2 3 4 5 6}]||J(nrow(time),2,1)||time[,{2 3 4 5 6}];
	
*	predicted excess hazard function;
	predeh1 = exp(des1*b1);

*	covariance matrix of predicted excess hazard;
	predeh1_vcov = diag(predeh1)*des1*v1*des1`*diag(predeh1);
	
*	cumulative excess hazard and variance estimate;	
	predceh1 = cummat*predeh1;		
	predceh1_var = vecdiag(cummat*predeh1_vcov*cummat`);
	
*	similar computations for oldest age group;
	predeh4 = exp(des4*b4);
	predeh4_vcov = diag(predeh4)*des4*v4*des4`*diag(predeh4);
	predceh4 = cummat*predeh4;		
	predceh4_var = vecdiag(cummat*predeh4_vcov*cummat`);

*	estimates and time variable for plotting;
	rs = j(nrow(time),7,.);
	rs[,1] = time[,1];
	rs[,2] = exp(-predceh1);
	rs[,3] = exp(-(predceh1 + 1.96*sqrt(predceh1_var)));
	rs[,4] = exp(-(predceh1 - 1.96*sqrt(predceh1_var)));
	
	rs[,5] = exp(-predceh4);
	rs[,6] = exp(-(predceh4 + 1.96*sqrt(predceh4_var)));
	rs[,7] = exp(-(predceh4 - 1.96*sqrt(predceh4_var)));
	
	create rs from rs[colname = {"midtime" "rs1" "rs1_lci" "rs1_uci" "rs4" "rs4_lci" "rs4_uci" }];
	append from rs ;
quit;

*	get life table estimates of relative survival from rel_surv ;
data obs;
	set grouped (keep = right left cr female agegrp year8594);
	where female = 0 and year8594 =1 and agegrp in(0 3);
	midtime = (right+left)/2;
	
	drop left right;
run;

proc sort data = obs;
	by agegrp female year8594 midtime;
run;

*	combine life table and regression estimates;
data rs;
	set rs obs;
run;

title 'Modelled and life table Relative survival estimates, males diagnosed 1985 - 1994';
proc sgplot data = rs;
	band x = midtime upper = rs1_uci lower = rs1_lci;
	band x = midtime upper = rs4_uci lower = rs4_lci;
	series y = rs1 x = midtime;
	series y = rs4 x = midtime;
	yaxis values= (0.4 to 1 by .2) label='Relative Survival';
	xaxis values = (0 to 6 by 1) label='Time since diagnosis (years)';
	scatter y = cr x = midtime/group=agegrp;
run;
