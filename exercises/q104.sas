/*
==================
 EXERCISE 104
 
 K-M survival, statistical tests
 	lifetest
 	%mrate()

reviewed:  4 May 2020
==================
*/

options fmtsearch = (data.melanoma_formats);

/* Data set used */
data melanoma;
	set data.melanoma (where = (stage=1));
run;

*	K-M by calendar time (censor deaths not due to cancer)
	log-rank, wilcoxon by default;
proc lifetest data = melanoma notable
	plots = (survival (nocensor), hazard);
	time surv_mm*status(0, 2 ,4);
	strata year8594;
run;

*	estimate mortality rates per (1,000) months by age group;
%mrate(in = melanoma, out = table,  byvar = agegrp, per = 1000,
	eventvar = status, eventlist = 1 , timevar = surv_mm);

*	use of scale parameter to show rates per year;
%mrate(in = melanoma, out = table,  byvar = agegrp, per = 1000,
	scale = 12,
	eventvar = status, eventlist = 1 , timevar = surv_mm);

title ;
*	K-M by age group;
proc lifetest data = melanoma notable
	plots = (survival (nocensor), hazard);
	time surv_mm*status(0, 2 ,4);
	strata agegrp;
run;

*	estimate mortality rates by sex;
%mrate(in = melanoma, out = table,  byvar = sex, per = 1000, scale = 12,
	eventvar = status, eventlist = 1 , timevar = surv_mm);

title ;
*	log-rank and Wilcoxon tests for variable sex;
proc lifetest data = melanoma notable 
	plots = none;
	time surv_mm*status(0, 2 ,4);
	test sex ;
run;
