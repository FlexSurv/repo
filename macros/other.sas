/*
	survival macros for exercises
*/


*	failure time with left truncation;
%global knots;

%macro hazard_late(data= , strat = 1, entry= , exit = , fail =, censor =);

%if &strat. = 1 %then %do;
	%let by_strat = ;
	%let ph_strat =;
	%let group = ;
	%let merge = %str(set _h_ ;
	if _n_  = 1 then do;);

%end;
%else %do;
	%let by_strat = by &strat.;
	%let ph_strat = strata &strat.;
	%let group = group=&strat.;
	%let merge = %str(merge _rep_ (keep = &strat. &entry.)
	_h_;
	&by_strat.;
	if first.&strat. then do;);

	proc sort data = &data.;
		&by_strat.;
	run;
%end;

data _temp_;
	set &data.;
	failed = &fail. not in(&censor.);
run;

proc means data = _temp_ noprint;
	&by_strat.;
	output out = _rep_  sum(failed) = failed min(&entry)= &entry. max(&exit.) = &exit.;
run;

proc print data = _rep_ noobs label split = '*';
	&by_strat.;
	var _freq_   failed &entry. &exit.;
	label failed = 'total*events'
		_freq_ = 'total*subjects'
		&exit. = 'last observed*death'
		&entry. = 'earliest entry*time';
run;
	
proc phreg data = &data. noprint;
	model &exit.*&fail.(&censor.) =/entry = &entry.;
	baseline out = _surv_ (where = (&exit. ^= 0))
		survival = survival;
	&ph_strat.;
run;

proc datasets force noprint;
	delete _temp_;
quit;
run;

%smooth(time = &exit., yscale = linear, option=plot, strat = &strat.);

%mend;

%macro smooth (data=_last_, out = _plt_, time=, width=, survival=survival,
yscale = linear, option=plot, strat = 1);

/*********************************************************************
MACRO SMOOTH produces graphs of smoothed hazard functions using output
from either PROC LIFETEST or PROC PHREG. With PROC LIFETEST, it uses the
data set produced by the OUTSURV option in the PROC statement. With PROC
PHREG, it uses the data set produced by the BASELINE statement. SMOOTH
employs a kernel smoothing method described by H. Ramlau-Hansen (1983),
"Smoothing Counting Process Intensities by Means of Kernel Functions,"
The Annals of Statistics 11, 453-466. If there is more than one survival
curve in the input data set, SMOOTH will produce multiple smoothed
hazard curves on the same axes.

There are four parameters:

DATA     is the name of the data set containing survivor function
         estimates. The default is the most recently created data set.

OUT		 name of output dataset.  defalt is _plt_

TIME     is name of the variable containing event times.

SURVIVAL is the name of a variable containing survivor function
         estimates (the default is SURVIVAL, which is the automatic name in
         PROC LIFETEST).

WIDTH    is bandwidth of smoothing function. The default is 1/5 of the range
         of event times.
         
YSCALE	linear (the default)
		log
		
option	plot/noplot

Example of usage:

%smooth(data=my.data,time=duration,width=8,survival=s)

Author:  Paul D. Allison, University of Pennsylvania
         allison@ssc.upenn.edu
         
modified (RAD 6 Jan 2017) to account for late entry
added a parameter for stratification variable
30 sept 2019	sort stratified survival dataset
*************************************************************************/
%if &strat. = 1 %then %do;
	%let by_strat = ;
	%let first = %str( _n_ = 1);
%end;
%else %do;
	%let by_strat = by &strat.;
	%let first = %str( first.&strat.);
	%let group = group=&strat.;
	proc sort data = &data.;
		by &strat.;
	run;
%end;

data _inset_;
 set &data end=final;
 	&by_strat.;
 retain _grp_ _censor_ 0;
 t=&time;
 survival=&survival;
 if &first. then _grp_=_grp_+1;
 keep _grp_ t survival;
 if final and _grp_ > 1 then call symput('nset','yes');
   else if final then call symput('nset','no');
 if _censor_ = 1 then delete;
 if survival in (0,1) then delete;
run;

proc iml;
use _inset_;
read all var {t _grp_};
%if &width ne %then %let w2=&width;
  %else %let w2=(max(t)-min(t))/5;
w=&w2;
z=char(w,8,2);
call symput('width',z);
numset=max(_grp_);
create &out. var{ lambda s group};
setin _inset_ ;
do m=1 to numset;
  read all var {t survival _grp_} where (_grp_=m);
  n=nrow(survival);
  lo=t[1];			* rad changed 17 jan, 2017;
  * lo=w;			* mgh changed 1 feb 2010;
/*  hi=t[n] - w;*/
  hi=t[n] - lo;		*	rad changed 1 june 2020;
  npt=50;
  inc=(hi-lo)/npt;
  s=lo+(1:npt)`*inc;
  group=j(npt,1,m);
  slag=1//survival[1:n-1];
  h=1-survival/slag;
  x = (j(npt,1,1)*t` - s*j(1,n,1))/w;
  k=.75*(1-x#x)#(abs(x)<=1);
  lambda=k*h/w;
  append;
end;
quit;

%if &nset = yes %then %let c=/group=group;
  %else %let c=;
  
 %let scale =;
 %let type =;
 %if %upcase(&yscale) = LOG %then %do;
 	%let type = type=log;
 	%let scale = Log;
 %end;
%if &option = plot %then %do;
proc sgplot data=_plt_;
  series y=lambda x=s &c ;
  yaxis &type. label="&scale. Hazard Function" ;
  xaxis  label= "Analysis Time (bandwidth=&width)" ;
  
run;
%end;

%mend smooth;




%macro mrate(in = _events_, out = _rates_, byvar = Total, scale = 1,
	eventvar = _death_, eventlist = 1, options = ,
	per = 1, timevar = _t_, level = 95);
	
*	calculate event rates from survival time in censored survival data;

*	default values of in =, eventvar = and eventlist = depend on having 
	used macro %stset to build the standard dataset;
	
*	required parameters:
	in		dataset with observations
	eventvar	defines events to count
	eventlist	values of eventvar to consider events
	timevar	time variable (years)
	
	optional parameters:
	byvar	classification variable (can be blank)
	out		dataset to hold computed rates, etc
	options	noprint is the only option currently
	per		divisor for person-years computation. default is 1
	scale	to change months to years (scale = 12)
			or days to years (scale = 365.24)
			default = 1
	level	level for confidence intervals.  default is
			95% confidence intervals;

	%if &in. eq 
		or &byvar. eq 
		or &eventlist. eq 
		or &timevar. eq %then %do;
			%put All required parameters have not been specified;
			%goto fin;
		%end;
			
*	set option defaults;
	%let print = 1;

*	evaluate options parameter;
	%if "&options." ne "" %then %do;
		%let options = %sysfunc(lowcase(&options.));
		%if %sysfunc(findw(noprint, &options.)) ^= 0 %then %let print = 0;
	%end;

	%if &byvar. eq Total %then %do;
		%let sel = "Total" as total,;
	%end;
	%if &byvar. ne Total %then %do;
		%let sel = &byvar. ,;
	%end;
	
	proc sql; 
	create table &out. as
		select 
			&sel.
			sum(ifn(&eventvar. in (&eventlist.),1,0)) as d,
			sum(&timevar.)/(&per.*&scale.) as y
		from &in.
		group by &byvar.;
		
	alter table &out.
		add rate num,
		cll num,
		clh num;
		
	update &out.
		set rate = d/y,
		 cll = ifn(d=0, ., exp(log(d/y) + probit(0.5 - &level./200)*sqrt(1/d))),
		 clh = ifn(d=0, ., exp(log(d/y) - probit(0.5 - &level./200)*sqrt(1/d)));
	
	quit;
		
	%if &print. %then %do;	

		title "event rates per &per. person-time";
	%if &scale. ne 1 %then %do; title2 "event time = &timevar./&scale."; %end;
		proc print data = &out. noobs label;
			var &byvar. d y rate cll clh;
			label d = 'Events'
				y = 'person time'
				rate = "Rate"
				cll = "lower &level.% CI"
				clh = "upper &level.% CI";
		run;
	%end;

*	clean up;
	%if &out. eq _rates_ %then %do;
		proc datasets library = work nolist ;
		delete _rates_ ;
		quit;
		run;
	%end;
	
*	error conditions branch to here;
%fin:
	
%mend;




