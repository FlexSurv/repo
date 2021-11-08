/*
==================
 EXERCISE 102
 
 life table, K-M survival
 procs in use
 	lifetest: life table method and K-M method

reviewed:  3 May 2020
==================
*/

options fmtsearch = (data.melanoma_formats);


* 	Create 0/1 outcome variable;
data melanoma;
	set data.melanoma (where = (stage = 1)) ;
	
	if status = 1 then csr_fail = 1;
	else if status ^= . then csr_fail = 0;
run;

* 	Life table;
proc lifetest data = melanoma 
	method = lt
	plots = none
	width = 1;
	time surv_yy*csr_fail(0);
run;

proc lifetest data = melanoma
	method = lt
	plots = none;
	time surv_mm*csr_fail(0);
run;

* Kaplan-Meier (survival time in years);
proc lifetest data = melanoma
	plots = survival(atrisk nocensor);
	time surv_yy*csr_fail(0);
run;

* Kaplan-Meier (survival time in months);
proc lifetest data = melanoma notable
	plots = survival(atrisk nocensor);
	time surv_mm*csr_fail(0);
run;
