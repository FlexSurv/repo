/*
==================
 EXERCISE 250
 
 Probability of death in a competing risks framework (life table relative survival) 
 
 %rel_surv

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

*	(a) Load melanoma data and produce life tables by age group and sex ;

title 'crude probabilities of death, melanoma data, stage 1';
%rel_surv(infile = data.melanoma (where = (year8594 = 1)), 
	patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), 
	intervals = %str(0 to 2.75 by .25, 3 to 5 by 1),
	popmort = data.popmort, mergeby = _age sex _year, 
	list = l d w cp_p cp_star_w cr_p cgc cgo,
	strat = agegrp sex, 
	crude= mod(right,1) = 0);
	
*	(b) similarity between F and cgc (death due to cancer) in youngest group;
data grouped;
	set grouped;
	f = 1-cp_p;			/*	obverse of observed (weighted) survival	*/
	f2 = cgc+cgo;
	prob_c = cgc/f;
	prob_o = cgo/f;
	net = 1-cr_p;		/*	obverse of PPE (net) survival	*/
	label f 	= 'cumulative*failure'
		f2 		= 'all cause prob'
		prob_c 	= 'prop. deaths*due to cancer'
		prob_o 	= 'prop. death*due to other causes'
		net 	= 'net prob*death*due to cancer';
run;

title 'melanoma crude probability of death for youngest age group';
title2 'all cause mortality and deaths due to cancer';
proc print data = grouped noobs label split = '*';
	where agegrp= 0 and sex = 1;
	by agegrp;
	id agegrp;
	var right sex f cgc;
run;

*	(c) Relationship between crude probability and all-cause probability;
title 'melanoma crude probabilities and all-cause mortality';
proc print data = grouped noobs label split = '*';
	where agegrp=2 and right = 5;
	by agegrp; id agegrp;
	var right sex f cgc cgo f2;
run;

*	(d) Proportion of deaths due to cancer and other causes;
title 'melanoma crude probabilities and all-cause mortality';
proc print data = grouped noobs label split = '*';
	where agegrp=2 and right = 5;
	by agegrp; id agegrp;
	var right  sex f cgc cgo prob_c prob_o;
run;

*	(e) Plot overall, net and crude probability of death by age group;
title 'melanoma death due to cancer, compared to net and all-cause mortality';
proc sgpanel data = grouped;
	where sex = 1;
	panelby agegrp/novarname;
	series x = right y = F/legendlabel='Overall';
	series x = right y = net/legendlabel='Net';
	series x = right y = cgc/legendlabel='Cancer';
	rowaxis label = 'end of interval';
	colaxis label= 'probability of death';
run;
	



