/*
==================
 EXERCISE 242
 
 Age standardization using flexible parametric models.
 
 %stset
 %stpm2
 %predict

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

*	(a)	select early stage and truncate follow-up time to 10 years;
*	fit null model and predict relative survival;
data melanoma10;
	set data.melanoma (where = (stage=1));
	
	if surv_mm > 120.5 then do;
		surv_mm = 120.5;
		if status in(1 2) then status = 0;
	end;
	
	surv_mm = surv_mm/12;	

	agegrp2 = agegrp = 1;
	agegrp3 = agegrp = 2;
	agegrp4 = agegrp = 3;
run;

proc sql;
	create table melanoma_bh as select a.*,
		b.rate
		from melanoma10 a left join data.popmort b
		on floor(a.age + a.surv_mm) = b._age
		and a.sex = b.sex
		and floor(yydx + a.surv_mm) = b._year;
quit;


*	fit base model (no covariates);
%stset(melanoma_bh, status(1 2), surv_mm, id);
%stpm2( , df=5, scale=hazard, bhazard=rate );

*	create temporary time variable for estimation;
%range(temptime, 0, 10, 101);
%predict(rs0, survival, timevar = temptime, options = ci);


*	show estimated survival at 10 years;
title 'overall relative survival at 10 years';
proc print data = _events_ noobs split = '*';
	where temptime = 10;
	var rs0 rs0_lci rs0_uci;
	label rs0 = 'Overall*survival'
		rs0_lci = 'lower 95% CI'
		rs0_uci = 'upper 95% CI';

	format rs0 rs0_lci rs0_uci 5.3;
run;


*	(b) Fit PEH model and predict relative survival in each age group;
%stpm2(agegrp2 agegrp3 agegrp4 , df=5, scale=hazard, bhazard=rate, options = eform);

*	predictions for each age group;
%predict(rs1, survival, at = zero, timevar = temptime);

%predict(rs2, survival, at = agegrp2:1 zero, timevar = temptime);

%predict(rs3, survival, at = agegrp3:1 zero, timevar = temptime);

%predict(rs4, survival, at = agegrp4:1 zero, timevar = temptime);

title 'age-specific and all ages (crude) survival estimates';
proc sgplot data = _events_;
	series x = temptime y = rs0/legendlabel='All ages' lineattrs=(color=black pattern = dash);
	series x = temptime y = rs1/legendlabel='0 - 44' lineattrs=(color=red);
	series x = temptime y = rs2/legendlabel='45 - 59' lineattrs=(color=blue);
	series x = temptime y = rs3/legendlabel='60 - 74' lineattrs=(color=green);
	series x = temptime y = rs4/legendlabel='75 +' lineattrs=(color=orange);
	xaxis label='Years since diagnosis';
	yaxis label='Relative Survival';
run;

*	c) Tabulate agegrp and obtain weighted average of 4 relative survival curves;
proc freq data = _events_;
	table agegrp;
run;

data _events_;
	set _events_;
	rs_stand1 = 0.2751*rs1 + 0.2962*rs2 + 0.2888*rs3 + 0.1399*rs4;
run;

title 'age-specific, all ages and standardised survival estimates';
proc sgplot data = _events_;
	series x = temptime y = rs0/legendlabel='All ages' lineattrs=(color=black pattern = dash);
	series x = temptime y = rs_stand1/legendlabel='Age standardised' lineattrs=(color=black pattern = dot);
	series x = temptime y = rs1/legendlabel='0 - 44' lineattrs=(color=red);
	series x = temptime y = rs2/legendlabel='45 - 59' lineattrs=(color=blue);
	series x = temptime y = rs3/legendlabel='60 - 74' lineattrs=(color=green);
	series x = temptime y = rs4/legendlabel='75 +' lineattrs=(color=orange);
	xaxis label='Years since diagnosis';
	yaxis label='Relative Survival';
run;

*	(d) Age standardized relative survival at 10 years;
title 'age standardised estimate at 10 years';
proc print data = _events_ noobs;
	where temptime = 10;
	var rs_stand1;
run;

*	 (e) Use the meansurv option;
%predict(rs_stand2, meansurv, timevar = temptime, options=ci);

title 'alternate versions of age standardisation';
proc sgplot data = _events_;
	series x = temptime y = rs_stand1/legendlabel='Age standardised (using weights)' lineattrs=(color=black pattern = dash);
	series x = temptime y = rs_stand2/legendlabel='Age standardised (meansurv option)' lineattrs=(color=red pattern = dot);
	xaxis label='Years since diagnosis';
	yaxis label='Relative Survival';
run;

*	(f) Obtaining confidence intervals;
title "relative survival at 10 years with 95% CI from 'meansurv' option";
proc print data = _events_ noobs;
	where temptime = 10;
	var rs_stand2 rs_stand2_lci rs_stand2_uci;
run;

*	(g) Fit a PEH model for age group and calendar period;
%stpm2( agegrp2 agegrp3 agegrp4 year8594, scale=hazard, df=5, bhazard=rate);
%predict(rs, survival);

*	compute average rel survival across age groups and time periods;
title 'average survival by age and time period';
proc means data = _events_ mean nway;
	class agegrp year8594;
	var rs;
run;

*	(h) Change in age distribution between calendar periods;
title 'age distribution in dataset by time period';
proc freq data = _events_;
	table agegrp*year8594/norow nopercent;
run;

*	(i) Predict relative survival for each calendar period
	Standardize to age distribution in first calendar period
	by restricting estimation dataset to those in earlier period;
%predict(rs_7584, meansurv, at=year8594:0, timevar = temptime, ifp = year8594=0);

%predict(rs_8594, meansurv, at=year8594:1, timevar = temptime, ifp = year8594=0);

title 'age standardised survival estimates (earlier period provides internal weights)';
proc sgplot data = _events_;
	series x = temptime y = rs_7584/legendlabel='1975 - 1984' lineattrs=(color=black);
	series x = temptime y = rs_8594/legendlabel='1985 - 1994' lineattrs=(color=red);
	xaxis label='Years since diagnosis';
	yaxis label='Relative Survival';
run;

proc print data = _events_ noobs;
	where temptime = 10;
	var rs_7584 rs_8594;
run;

*	(j) Standardize to age distribution in second calendar period ;
%predict(rs_8594b, meansurv, at=year8594:1, timevar = temptime, ifp = year8594=1);

title 'age standardised survival estimates (different reference periods)';
proc sgplot data = _events_;
	series x = temptime y = rs_7584/legendlabel='1975 - 1984 (Reference 1975-1984)' lineattrs=(color=black);
	series x = temptime y = rs_8594/legendlabel='1985 - 1994 (Reference 1975-1984)' lineattrs=(color=red);
	series x = temptime y = rs_8594b/legendlabel='1985 - 1994 (Reference 1985-1994)' lineattrs=(color=blue);
	xaxis label='Years since diagnosis';
	yaxis label='Relative Survival';
	keylegend /down=3  location=inside;
run;

title 'age standardised survival estimates at 10 years (different reference periods)';
proc print data = _events_ noobs split = '*';
	where temptime = 10;
	var rs_7584 rs_8594 rs_8594b;
	label rs_7584 = 'earlier period*(ref=1975-84)'
		rs_8594 = 'later period*(ref=1975-84)'
		rs_8594b = 'later period*(ref=1985-94)';
run;
