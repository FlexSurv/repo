/*
==================
 EXERCISE 122
 
 cox regression for cause-specific mortality
 phreg

reviewed: 4 May 2020
====================
*/


options fmtsearch = (data.melanoma_formats);

/* Data set used */

*	truncate follow-up time to 10 years;
data melanoma10;
	set data.melanoma (where = (stage=1));
	
	if surv_mm >=120 then do;
		surv_mm = 120;
		if status in(1,2) then status = 0;
	end;
run;

proc means data = melanoma10 min max;
	var surv_mm;
run;

proc freq data = melanoma10;
	table status;
run;

*	all-cause survival;
proc phreg data = melanoma10;
	class agegrp sex (ref = last) year8594 /ref=first;
	model surv_mm*status(0 4) = sex year8594 agegrp /risklimits;
run;

*	cause-specific survival (using follow-up truncated at 10 years);
proc phreg data = melanoma10;
	class agegrp sex (ref = last) year8594 /ref=first;
	model surv_mm*status(0 2 4) = sex year8594 agegrp /risklimits;
run;

*	cause-specific survival (using all follow-up events);
proc phreg data = data.melanoma;
	class agegrp sex (ref = last) year8594 /ref=first;
	model surv_mm*status(0 2 4) = sex year8594 agegrp /risklimits;
run;

