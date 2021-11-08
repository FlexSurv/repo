/*
==================
 EXERCISE 282		
 REVISED MAY 2015	
	avoidable deaths	
	%rel_surv	

reviewed: 7 May 2020
====================
*/

options fmtsearch = (data.melanoma_formats);

* 	(a) Five year relative survival by age group and sex for melanoma data	;
%rel_surv(infile = data.melanoma (where = (year8594=1)),
	patientid=id,
	age = age, 
	exit = surv_mm,
	yydx = yydx,
	censor = status(0 4),
	scale = 12,
	strat =  agegrp sex ,
	popmort = data.popmort,
	intervals = %str(0 to 1 by 1/12, 2 to 3 by .25, 4 to 5 by 1),
	crude = right =5);

*	sort and retain just the elements we need 
	sort so females come before males;
proc sort data = grouped 
		(keep = fu right agegrp sex l cp_star cp cr);
	by agegrp fu  sex  ;
run;


* 	(b) estimated Number of incident cases per year 
	= total cases in each age-sex group at diagnosis / 10;
data grouped_wide;
	set grouped;
	by agegrp fu  sex  ;
	
	if fu = 1 then do;
		if sex = 1 then m_Nrisk = l/10;
		else f_Nrisk = l/10;
	end;
		
	if first.fu then do;
		m_l = l;
		m_cp_star = cp_star;
		m_cp = cp;
		m_cr = cr;
	end;
	
	retain m_Nrisk f_Nrisk m_l m_cp_star m_cp m_cr;
	
	if last.fu then do;
		f_l = l;
		f_cp_star = cp_star;
		f_cp = cp;
		f_cr = cr;
		
*	now do calculations of observed and expected deaths;		
		
		m_p_dead 	= 1-m_cp_star*m_cr;	
		m_Nd 		= m_Nrisk*m_p_dead;  
		m_NExp_d 	= m_Nrisk*(1-m_cp_star);
		m_ED 		= m_Nd - m_Nexp_d;
		
		f_p_dead 	= 1-f_cp_star*f_cr;	
		f_Nd 		= f_Nrisk*f_p_dead;  
		f_NExp_d 	= f_Nrisk*(1-f_cp_star);
		f_ED 		= f_Nd - f_Nexp_d;
		
*	calculation of avoidable deaths for this interval;		
		Nd_m_f 		= m_Nrisk*(1 - m_cp_star*f_cr);
		AD_m 		= m_Nd - Nd_m_f;
		output ;
	end;

	keep agegrp right m_Nrisk f_Nrisk m_p_dead m_Nd m_NExp_d m_ED f_p_dead f_Nd f_NExp_d f_ED
		AD_m;
run;

*	print results:  first, observed and expected deaths;
proc print data = grouped_wide noobs ;
	where right = 5;
	by agegrp ;
	id agegrp ;
	var  m_Nrisk  m_p_dead m_Nd m_NExp_d m_ED f_Nrisk  f_p_dead f_Nd f_NExp_d f_ED;
	format m_: f_: 5.1;
run;

proc print data = grouped_wide noobs;
	where  mod(right,1) = 0 and agegrp = 3;
	by agegrp; id agegrp;
	var  right AD_m ;
	format ad_m 5.1;
run;
