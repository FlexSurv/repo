/*
==================
 EXERCISE 210
 
 Modelling excess mortality using Poisson regression 
 
 nlmixed genmod 
 %rel_surv()

reviewed:  5 May 2020
=====================

*/

options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;

data melanoma;
	set data.melanoma (where = (stage = 1));
run;

*	Estimate relative survival for each combination of of the covariates;
%rel_surv(infile = melanoma , patientid=id,
	age = age, sex=sex, origin=dx, exit = exit, scale = 365.24,
	censor = status(0,4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat=sex agegrp year8594);
	
*	Restrict to first 5 years of follow-up;
data grouped;
	set grouped; 
	where right <=5;
run;

data _cases_;
	set _cases_; 
	if right <= 5;
run;

*	grouped poisson regression of excess risk (needs user-defined link functions);
proc genmod data = grouped;
	class right agegrp sex year8594/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right agegrp sex year8594 /dist=poisson offset = ln_y;
	ods output parameterestimates = grouped_parms;
run;

*	d) poisson regression with interaction;
proc genmod data = grouped;
	class right agegrp sex year8594/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right|agegrp sex year8594 /dist=poisson offset = ln_y;
 run;

*	enter -2LogLikelihood value;
data lrt_pval;
	lrt = 2*abs(2088.7339-2091.4340);
	df = 12;							*	in this case;
	p_value = 1 - probchi(LRT,df);
run;

proc print data=lrt_pval;
	title1 "LR test statistic and p-value";
run;

*	e) individual survival data containing individual subject-band observations;
proc genmod data = _cases_ ;
	where right < 6;
	class right agegrp sex year8594/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right agegrp sex year8594 /dist=poisson offset = ln_y;
   	ods output parameterestimates = individ_parms;
 run;
 
*	esteve - define custom likelihood function, maximise with nlmixed;
*	need dummy variables for interval, age group;
data esteve;
	set _cases_;
	where right <=5;
	
	int2 = right = 2;
	int3 = right = 3;
	int4 = right = 4;
	int5 = right = 5;
	
	a2 = agegrp = 1;
	a3 = agegrp = 2;
	a4 = agegrp = 3;
	
	cons = 1;
	sex = 2-sex;	*	recode sex to 0/1;
	
	len = 1/length;
run;

*	check we have created the right age dummies;
proc freq data = esteve;
	table agegrp*(a2 a3 a4)/norow nopercent nocol;
run;

	
*	(f) using proc nlp to fit esteve model
	this proc is not available in sas UE;
*	run one or the other (or both),  should be the same;
 proc nlp data = esteve; 
 	parms int i2-i5 ab2-ab4 s yr; 
 	 
 	theta = int + i2*int2 + i3*int3 + i4*int4 + i5*int5 + 
 			ab2*a2 + ab3*a3 +ab4*a4+ 
 			s*sex + yr*year8594; 
 			 
 	lnf = -exp(theta)*y + (log(-log(p_star)/length + exp(theta)))*d;
 		 
 	max lnf; 
 run; 

*	Sas UE makes proc nlmixed available;
proc nlmixed data = esteve;
	theta = int + i2*int2 + i3*int3 + i4*int4 + i5*int5 +
			ab2*a2 + ab3*a3 +ab4*a4+
			s*sex + yr*year8594;
			
	lnf = -exp(theta)*y + (log(-log(p_star)/length + exp(theta)))*d;
		
	model d ~ general(lnf); 
	ods output parameterestimates = esteve_parms;
run;

data RER;
	set esteve_parms;
	
	rer = exp(estimate);
run;

proc print data = rer;
	var parameter rer;
run;

*	(g) Hakulinen-Tenkanen;
/*************************************************************************
Grouped survival times
Binomial error structure (Hakulinen-Tenkanen)
This approach should not be used for modelling period estimates
**************************************************************************/
ods output parameterestimates=Hakulinen_parms /* parameter estimates */
           modelinfo=modelinfo        /* Model information */
           modelfit=modelfit          /* Model fit information */
           convergencestatus=converge /* Whether the model converged */
           ModelANOVA=type3estimates;      /* Type III estimates */

proc genmod data=grouped(where=(right le 5)) order=formatted;
title2 'Binomial error model fitted to grouped data [model 1 in Dickman et al. (2004)]';
title3 'Main effects model (follow-up, sex, age, and dgnyear)';
variance var = _mean_*(1-_mean_/l_prime);
fwdlink link = log(-log(_mean_/p_star));
invlink ilink = exp(-exp(_xbeta_))*p_star;
class right sex agegrp year8594/ref=first;
model ns/l_prime = right agegrp sex year8594 / error=bin type3;

run;

ods output close;

data all_parmest;
merge grouped_parms (keep = parameter level1 estimate stderr 
			where = (grp_est ^= 0)
			rename = (estimate = grp_est stderr = grp_stderr ))
	individ_parms (keep = parameter estimate stderr 
			where = (ind_est ^= 0)
			rename = (estimate = ind_est stderr = ind_stderr ))
	esteve_parms (keep = parameter estimate standarderror 
			rename = (estimate = est_est standarderror = est_stderr ))
	Hakulinen_parms (keep = parameter estimate stderr 
			where = (hak_est ^= 0)
			rename = (estimate = hak_est stderr = hak_stderr ));
			
	array rer grp_est ind_est est_est hak_est;
	do over rer;
		rer = exp(rer);
	end;
	
	if parameter = 'Intercept' then 
		do over rer;
			rer = 100*rer;
		end;
	length parm $32;
	parm = parameter||level1;
	if parameter = 'Intercept' then parm = 'Intercept (*100)';
run;

title 'Excess hazard ratios for various models';
proc print data = all_parmest noobs label;
	where parameter ^= 'Scale';
	var parm grp_est ind_est est_est hak_est;
	label
		parm = 'Parameter' 
		grp_est = 'Grouped'
		ind_est = 'Individual'
		est_est = 'Esteve'
		hak_est = 'Hakulinen';
	format grp_est ind_est est_est hak_est 6.3;
run;

*	h) all stages included;
data melanoma; set data.melanoma; run;

%rel_surv(infile = melanoma ,  patientid=id,
	age = age, sex=sex, yydx = yydx, exit = surv_mm, scale = 12,
	censor = status(0,4), intervals = 0 to 10 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	strat=sex agegrp year8594 stage,
	crude = %str(scan(interval,-1,' ') in('3.00','10.0')));

*	Restrict to first 5 years of follow-up;
data grouped;set grouped; 
	right = input(scan(interval,-1,' '),5.2);
	if right <= 5;
run;

*	main effects model;
proc genmod data = grouped;
	class stage right agegrp sex year8594 stage/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right stage agegrp sex year8594 /dist=poisson offset = ln_y;
run; 

*	all main effects with interaction of period by stage ;
proc genmod data = grouped;
	class stage right agegrp sex year8594 stage/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right|stage agegrp sex year8594 /dist=poisson offset = ln_y;
run;

*	interaction but no main effect of stage; 
proc genmod data = grouped;
	class stage right agegrp sex year8594 stage/ref=first;
	
	fwdlink link = log(_MEAN_-d_star);		*	link function;
	invlink ilink = exp(_XBETA_)+d_star;	*	inverse of link function;
	
   	model d = right right*stage agegrp sex year8594 /dist=poisson offset = ln_y;
run; 
