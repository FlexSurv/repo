/*
==================
 EXERCISE 103
 
 survival and hazard functions
 procs in use
 	lifetest
 	print
 	freq
 	%mrate

reviewed:  May the fourth be with you, 2020
==================
*/

options fmtsearch = (data.melanoma_formats);

proc freq data = data.melanoma;
	table stage ;
run;

* Survival and hazard function by stage - Kaplan Meir typically used for clinical data analyses ;
proc lifetest data = data.melanoma
	notable
	plots(only) = (survival (nocensor) hazard);
	time surv_mm*status(0 2 4);
	strata stage;
run;

*	average mortaiity rates over that period  ;
%mrate(in = data.melanoma, byvar = stage, 
	eventvar = status, eventlist = 1, timevar = surv_mm);

*	convert months to years of person-time, scale to 1,000 p-y;
data melanoma;
	set data.melanoma;
	surv_mm = surv_mm/12;
	label surv_mm = 'Survival time in years';
run;

*	mortaity rates per 1,000 persons years;
%mrate(in = melanoma, out = table,  byvar = stage, per = 1000,
	eventvar = status, eventlist = 1 , timevar = surv_mm);

title ;
* Survival and hazard function by sex;
proc lifetest data = melanoma 
	notable
	plots = (survival (nocensor) hazard);
	time surv_mm*status(0 2 4);
	strata sex;
run;

*	mortaity rates per 1,000 persons years;
%mrate(in=melanoma, byvar=sex, per=1000, eventvar=status, 
	eventlist=1, timevar=surv_mm);

*	look at status variable;
proc format library = data.melanoma_formats 
	cntlout=list;run;

title ;
proc print data = list;
	where lowcase(fmtname) in( 'statusa','agegrpa');
	var fmtname end label;
run;

proc freq data = melanoma;
	table status agegrp;
run;

*	all-cause survival;
proc lifetest data = melanoma	
	notable 
	plots = survival(nocensor);
	time surv_mm*status(0,4);* alive and lost to follow-up and call those 'censored';
	strata stage;
run;

*	compare cause-specific and all cause srurvival;
title 'cancer-specific deaths, oldest age group';
proc lifetest data = data.melanoma (where = (agegrp=3))	
	notable 
	plots = survival(nocensor);
	time surv_mm*status(0, 2, 4); * and now you are accepting, alive (0), died not of cancer(2) and lost to follow-up (4) .. as 'cencesor;
	strata stage;
run;

title 'all causes deaths, oldest age group';
proc lifetest data = data.melanoma (where = (agegrp=3))	
	notable 
	plots = survival(nocensor);
	time surv_mm*status(0 );
	strata stage;
run;

*	compare cause-specific and all cause srurvival by age;
proc lifetest data = data.melanoma 	
	notable 
	outsurv = survcanc
	plots = survival(nocensor);
	time surv_mm*status(0 2 4);
	strata agegrp;
run;

proc lifetest data = data.melanoma 
	notable 
	outsurv = survall
	plots = survival(nocensor);
	time surv_mm*status(0 );
	strata agegrp;
run;

data plot;
	set survcanc (in = a)
		survall;

	length cause $16;

	if a then cause = 'Cancer';
	else cause = 'All Causes';
run;
title 'cancer and all causes survival by age group';
proc sgpanel data = plot;
	panelby cause/novarname;

	series x = surv_mm y = survival/group = agegrp;
run;
