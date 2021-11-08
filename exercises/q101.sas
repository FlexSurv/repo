/*
==================
 EXERCISE 101
 
 life table, K-M survival
 procs used:
 	import
 	lifetest

reviewed:	3 May2020
==================
*/

options fmtsearch = (data.colon_formats);

* 	Create 0/1  variable for cause-specific (cancer, in this case) outcome;
data colon;
	set  data.colon_sample;
	if status = 1 then csr_fail = 1; *cause-specific survival failure - here they have failed because they have died of cancer, so =1;
	else if status ^= . then csr_fail = 0;
run;

* 	Life tables(lt) (1 year intervals) - using survival in years;
proc lifetest data = colon 
	method = lt 
	plots = none
	width = 1; * i.e. width of interval which in this case is 1 year;
	time surv_yy*csr_fail(0); * you are telling SAS that these are censured events;
run;

* 	Life tables (1 yeara intervals) - using survival in months;
proc lifetest data = colon
	method = lt
	plots = none
	width = 12;
	time surv_mm*csr_fail(0);
run;

* Kaplan-Meier (survival time in months);
proc lifetest data = colon
	plots = survival(atrisk nocensor);
	time surv_mm*csr_fail(0);
run;

* Kaplan-Meier (survival time in months);
title 'Kaplan-Meier estimate of cause-specific survival';
proc lifetest data = colon notable
	plots = survival(atrisk nocensor);
	time surv_mm*csr_fail(0);
*	strata sex;
	run;
