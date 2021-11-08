/*
==================
 EXERCISE 202
 
 Life table estimates of cause-specific survival using strs  
 
 lifetest sgplot
 %rel_surv

reviewed:  4 May 2020
=====================

*/
options fmtsearch = (data.melanoma_formats);

*	Load the Melanoma data, keep those with localized stage;


data melanoma;
	set data.melanoma (where = (stage=1)); 
	
*	Create failure indicator variable;
	csr_fail = (status = 1);
	rsr_fail = (status in(1 2));
	surv_mm_12 = surv_mm/12;
	label csr_fail = 'deaths due to cancer'
		rsr_fail = 'all deaths';
run;

proc freq data = melanoma;
	table status csr_fail rsr_fail;
run;

*	(a) Life-table estiamtes of Cause specific mortality, annual intervals 0-20 years
	saving (cause specific) estimates;
title 'Cause specific mortality: rel_surv macro';
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex, exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0 2 4), intervals = 0 to 20 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude_estimates = csr,
	list = l d w p cp, crude=1);
	
*	(b) Life table estimates using lifetest with actuarial method;
title 'Cause specific mortality: lifetest';
proc lifetest data = melanoma 
	plots = none
	outsurv= actuarial
	method=act
	width=12
	maxtime=240;
	time surv_mm*status(0 2 4);
run;

*	(c) estimate relative survival rates;
title 'Relative survival:  rel_surv macro';
%rel_surv(infile = melanoma,  patientid=id,
	age = age, sex=sex,  exit = surv_mm, scale = 12,
	yydx = yydx, censor = status(0,4), intervals = 0 to 20 by 1,
	popmort = data.popmort, mergeby = _age sex _year, 
	crude_estimates = rsr,
	list = l d w p cp cr cr_p, crude=1);
	

*	(d) combine cause-specific and relative survival estimates and display;
*	relative survival estimates;
data rsr;
	set rsr;
	
	SE_RSR=se_cr/cr;
	RSR=cr;
	net=cr_p;

run;

*	cause-specific estimates;
data csr;
	set csr;
	
	SE_cSR=se_cp/cp;
	cSR=cp;
run;

*	compare cause-specific and relative survival;
data comp;
	merge csr rsr;
	by fu;
run;

title 'comparison of relative, net (PPE) survival, and cause-specific survival estimates';
proc print data = comp noobs label;
	var interval csr rsr net;
	label csr = 'Cancer-specifc survival'
		rsr = 'Relative survival (E2)'
		net = 'Net survival (PPE)';
	format csr rsr net 5.3;
run;

title 'comparison of relative, net (PPE) survival, and cause-specific survival estimates';
proc sgplot data = comp;
	series x = right y = csr/legendlabel = 'Cancer-specific';
	series x = right y = rsr/legendlabel = 'Relative survival (E2)';
	series x = right y = net/legendlabel = 'Net survival';
	yaxis label = 'Cumulative survival' values = (0 to 1 by .2);
run;
