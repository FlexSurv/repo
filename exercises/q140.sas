/*
==================
 EXERCISE 140
 Probability of death in a competing risks framework (cause-specific survival)
 lifetest phreg (fine & gray)

 sgplot sgpanel
 %stpm2 %stpm2cif

reviewed:  4 may 2020
=====================


*/

options fmtsearch = (data.colon_formats);


proc format ;
    value female
    0 = 'Males'
    1 = 'Females';
    
    value sex
    1 = 'Males'
    2 = 'Females';
run;


*   Load the colon data, exclude those with missing stages
	convert months survival to years survival;
data colon;
    set data.colon (where = (stage ^=0));
    
*   no id variable on the dataset provided;
    id = _n_;
    
    if surv_mm > 120.5 then do;
        surv_mm = 120.5;
        status = 0;
    end;
    
    surv = surv_mm/12;
        
    female = sex=2;
    format female female.;
run;

*   (a) Obtain the Kaplan-Meier survival estimate for each of the two causes 
	cancer (status = 1) and non-cancer causes (status = 2).  other values (0, 4)
	are censored events

    plot results for men;

*   cancer cause of death;
title 'cancer cause of death';
proc lifetest data = colon notable outsurv=cancer
    plots=none maxtime=10;
    time surv*status(0 2 4);
    strata female;
run;

*   other causes of death;
title 'other (on-cancer) cause of death';
proc lifetest data = colon notable outsurv=other
    plots=none maxtime=10;
    time surv*status(0 1 4);
    strata female;
run;

data failure;
    set cancer (rename = (survival = s_cancer))
        other (rename = (survival = s_other));
        
    f_cancer = 1-s_cancer;	*	failure is 1-survival;
    f_other = 1-s_other;
run;

title 'independent causes of death from K-M analysis, males only';
proc sgplot data = failure;
    where female=0;
    series x = surv y = f_cancer/legendlabel='Cancer';
    series x = surv y = f_other/legendlabel='Other';
    
    yaxis label= 'Probability of death';
    xaxis label= 'time since diagnosis (years)';
run;

*   (b) Obtain non-parametric cumulative incidence function 
    for each cause (kinda slow...) ;
title ;
proc lifetest data = colon notable
    outcif = cif_sex
    plots=none maxtime=10;
    time surv*status(0 4)/failcode;
    strata female;
run;

data cif_plot1;
    set cif_sex (keep = surv cif female failcode)
        failure;
    if failcode = 1 then cancer_cif = cif;
    else if failcode = 2 then other_cif = cif;  
run;

*	plotting as smooth series, as step plot requires more than default resources;
title 'non-parametric CIF plot (dashed lines) from competing risks (Fine and Gray)';
title2 'compared to independent failure plots from K-M  (males only)';
proc sgplot data = cif_plot1;
    where female = 1;
    series x = surv y = cancer_cif/legendlabel='Cancer cif'  lineattrs=(pattern = dash) ;
    series x = surv y = f_cancer/legendlabel='Cancer K-M';
    series x = surv y = other_cif/lineattrs=(pattern = dash) legendlabel='Other cif';
    series x = surv y = f_other/legendlabel='Other K-M';
    
    yaxis label= 'Probability of death';
    xaxis label= 'time since diagnosis (years)';
    
    keylegend /down=2;
run;


*   (c) Estimate CIFs by age group;
title 'lifetest by age';
proc lifetest data = colon notable
    outcif = cif_age
    plots=none maxtime=10;
    time surv*status(0 4)/failcode;
    strata agegrp;
run;

title 'non-parametric CIF by age group';
proc sgpanel data = cif_age;
    panelby failcode/novarname;
    step x = surv y = cif/group=agegrp;
    format failcode status.;
    rowaxis label='Probability of death';
    colaxis label= 'time since diagnosis (years)';
run;

*   (d) CIF by Stage group;
title 'lifetest by stage';
proc lifetest data = colon notable
    outcif = cif_stage
    plots=none maxtime=10;
    time surv*status(0 4)/failcode;
    strata stage;
run;

title 'non-parametric  CIF by stage group';
proc sgpanel data = cif_stage;
    panelby failcode / novarname;
    step x = surv y = cif/group=stage ;
    format failcode status.;
    rowaxis values= (0 to 1 by .2) label='Probability of death';
    colaxis label= 'time since diagnosis (years)';
run;

*   (e) Fit a competing risks model using Fine and Gray's method;
*   dataset to request prediction from PHREG;
data cov;
    sex = 1; output;
    sex = 2; output;
    format sex sex.;
run;

*	CIF from PH model.  No covariates specified for the baseline output
	so covariates held at reference (males) level
	CIF for event 1 (cancer cause) is requested;
title 'Proportional (sub)hazard model for cancer, with non-cancer the competing cause';
proc phreg data = colon ;
    model surv*status(0, 4) = sex/ eventcode = 1;
    hazardratio 'Subdistribution Hazards' sex / diff=pairwise;
    baseline  out = out0 cif = cif/rowid = sex;
run;

*	use PH assumption to compute baseline subhazard for females,
	using estimated subhazard ratio for females vs males;
data out0;
    set out0;
    cif_males = cif;
    cif_females = 1 - (1-cif_males)**exp(0.03394);
run;

*   Plot the associated cancer-specific cumulative incidence functions
	from 'first principles';
title 'First Principles (proportional subhazards plot of CIF by sex for cancer deaths)';
proc sgplot data = out0;
    step x = surv y = cif_males;
    step x = surv y = cif_females;
    yaxis label='Cause-specific cumulative incidence';
    xaxis label='Time since diagnosis (years)';
run;

*	now request separate CIF for males and females;
title ;
proc phreg data = colon ;
    model surv*status(0, 4) = sex/ eventcode = 1;
    hazardratio 'Subdistribution Hazards' sex / diff=pairwise;
    baseline covariates=cov out = out1 cif = cif/rowid = sex;
run;

title 'Baseline for covariates specified';
proc sgplot data = out1;
    step x = surv y = cif/group=sex;
    yaxis label='Cause-specific cumulative incidence';
    xaxis label='Time since diagnosis (years)';
run;


*   (f) Fit a competing risks model and plot the CIFs for deaths due to other causes than cancer ;
title ;
proc phreg data = colon ;
    model surv*status(0, 4) = sex/ eventcode = 2;
    hazardratio 'Subdistribution Hazards' sex / diff=pairwise;
    baseline covariates=cov out = out2 cif = cif/rowid = sex;
run;

title 'Other Causes';
proc sgplot data = out2;
    step x = surv y = cif/group=sex;
    yaxis label='Cause-specific cumulative incidence';
    xaxis label='Time since diagnosis (years)';
run;

data both_sex;
    set out1 out2 (in=a rename = (cif = cif_other));
    if a then cause = 'Other Causes';

    else cause = 'Cancer';
run;

title 'proportional subhazards CIF (Fine and Gray) by sex';
title2 'non-cancer causes are dashed lines';
proc sgplot data = both_sex;

    step x = surv y = cif_other/group=sex lineattrs = (pattern = dash);
    step x = surv y = cif/group=sex;
    yaxis label='Probability of death';
    xaxis label='Time since diagnosis (years)';
run;


*   (g) Fit a competing risks model using the flexible parametric approach;
*   (i) expand the data, one row for each cause;

*   Recode and set up data for competing risk analysis
    still need a unique ID for each row;
data colon_expanded;
    set colon;
    do cause = 1 to 2;
        uniqueID = 10*id+cause;
        cancer  = cause = 1;
        other   = cause = 2;
*   Event is the event indicator, coded like this:
    event is 1 if person died from cancer
    event is 1 if person died from other;
         event=(cause=status);
         output;
    end;
run;

title 'description of expanded colon data';
proc freq data = colon_expanded;
    table event;
run;

* Look at the created data;
proc print data = colon_expanded noobs;
    where id <= 4;
    by id;
    id id;
    var uniqueID status cause sex event ;
run;

*   (ii) Fit flexible parametric model for both causes simultaneously 
    constant effect of sex;
%stset(colon_expanded, event(1), surv, uniqueID);

*   Fit the stpm2 model assuming the effect of sex is the same for 
    both cancer and other causes
	separate baseline functions for cancer and non-cancer causes;
%stpm2( cancer other female, scale=hazard, tvc= cancer other, dftvc=4, 
    options = rcsbaseoff noint eform);
    
*   (iii) sex can have a different effect on cause;
data _events_;
    set _events_;
    fem_can = female*cancer;
    fem_other = female*other;
run;
%stpm2( cancer other fem_can fem_other, scale=hazard, tvc= cancer other, dftvc=4, 
    options = rcsbaseoff noint eform);

data test;
    b1 = 0.01879;
    b2 = -0.2195;
    var = 0.000582+0.002301 ;
    chi = (b1-b2)**2/var;
    df = 1;
    p = 1-probchi(chi, df);
run;
    
title 'Likelihood ratio test for addition of sex by cause of death interaction';
proc print data = test noobs;
    var chi df p;
run;

*   (h) predict CIFs for males and females;
*	first for males (sex covariates are held at zero by default, which is the Male code);
*	note output dataset produced:  cif_est  must be saved, as it will be overwritten by the next call;
%stpm2cif(cancermale othermale, cause1 = cancer:1, 
	cause2 = other:1);

%stpm2cif(cancerfemale otherfemale, cause1 = cancer:1 fem_can:1, 
    cause2 = other:1 fem_other:1);
    
*   combine with results from non-parametric analysis;
data cif_both; 
    set cif_sex (in=a keep = female surv cif failcode 
            rename = (cif = cancer_cif)
            where = (failcode = 1) ) 
        cif_sex (in=b keep = female surv cif failcode 
            rename = (cif = other_cif)
                where = (failcode=2))
        _events_ (in=c 
            keep = _newt cif_cancermale cif_othermale
            rename = (cif_cancermale = cif_cancer
                        cif_othermale = cif_other))
        _events_ (in=d 
            keep = _newt cif_cancerfemale cif_otherfemale
            rename = (cif_cancerfemale = cif_cancer
                        cif_otherfemale = cif_other));                  
    if c then female = 0;
    else if d then female = 1;
run;

title 'comparison of CIF estimtes: non-parametric vs flexible parametric';
proc sgpanel data = cif_both;
    panelby female/novarname;
    series x = _newt y = cif_cancer;
    series x = _newt y = cif_other/lineattrs=(color=red);
    step x = surv y = cancer_cif/
        lineattrs=(pattern=dash ) legendlabel='cancer CIF (n-p)';
    step x = surv y = other_cif/
        lineattrs=(pattern=dash color=red)  legendlabel='Other CIF (n-p)';
    colaxis label = 'Time Since Diagnosis (years)';
    rowaxis label= 'Probability of Death';
run;

*   (i) Stack the cumulative incidence functions and plot again;
data cif_comb;
    set cif_both 
        (keep = female _newt cif_cancer cif_other);
    if female = 0 then do;
        total1 = cif_cancer;
        total2 = total1+cif_other;
        output;
    end;    

    else do;
        total1 = cif_cancer;
        total2 = total1+cif_other;
        output;
    end;    
    keep female _newt total1 total2;
run;

title 'alternate display of flexible parametric estimates';
proc sgpanel data = cif_comb;
    panelby female/novarname;
    band x = _newt lower = total1 upper = total2/legendlabel='Other';
    band x = _newt lower = 0 upper = total1/legendlabel='Cancer';
    colaxis label = 'Time Since Diagnosis (years)';
    rowaxis values = (0 to 1 by .2) label= 'Probability of Death';
run; 

*   (j) Categorize age and create interactions with cause;
data _events_;
    set _events_;
    
    age0can = agegrp = 0 and cancer = 1;
    age1can = agegrp = 1 and cancer = 1;
    age2can = agegrp = 2 and cancer = 1;
    age3can = agegrp = 3 and cancer = 1;
    
    age0oth = agegrp = 0 and other = 1;
    age1oth = agegrp = 1 and other = 1;
    age2oth = agegrp = 2 and other = 1;
    age3oth = agegrp = 3 and other = 1;
run;

%stpm2( cancer other fem_can fem_other
        age1can age2can age3can age1oth age2oth age3oth,
        scale=hazard, tvc= cancer other, dftvc=3, 
    options = rcsbaseoff noint eform);

*   (k) predict CIFs;
%stpm2cif(cancermale_age0 othermale_age0, cause1=cancer:1, cause2=other:1); 

%stpm2cif(cancermale_age3 othermale_age3, cause1=cancer:1 age3can:1, 
    cause2=other:1 age3oth:1); 

data comb03;
	set _events_ ;

    age0 = cif_othermale_age0;
    age3 = cif_othermale_age3;
    death_cause= 'Other Causes';
    output;
    
    age0 = cif_cancermale_age0;
    age3 = cif_cancermale_age3;
    death_cause = 'Cancer';
    output;
    keep _newt age0 age3 death_cause;
run;

title 'including age by cause interactions';
proc sgpanel data = comb03;
    panelby death_cause/novarname;
    series x = _newt y = age0/legendlabel='<45';
    series x = _newt y = age3/legendlabel='75+';
    colaxis label = 'Time Since Diagnosis (years)';
    rowaxis values = (0 to .6 by .1) label= 'Probability of Death';
run;    

*   (l) allow for time-dependent effects for cancer;
title ;
%stpm2( cancer other fem_can fem_other
    age1can age2can age3can age1oth age2oth age3oth , scale=hazard,
    dftvc=cancer:4 other:4 3,
    tvc=cancer other fem_can age1can age2can age3can, 
    options = rcsbaseoff noint eform);
        
%stpm2cif( cancermale_age0_tvc othermale_age0_tvc, cause1=cancer:1,
   cause2=other:1);
    
%stpm2cif (cancermale_age3_tvc othermale_age3_tvc, cause1=cancer:1 age3can:1, 
   cause2=other:1 age3oth:1); 
 
data comb03_tvc;
	set _events_ ;

    if cif_othermale_age0_tvc ^= . then do;
        age0_tvc = cif_othermale_age0_tvc;
        age3_tvc = cif_othermale_age3_tvc;
        death_cause = 'Other Causes';
        output;
    end;
    

    if cif_cancermale_age0_tvc ^= . then do;
        age0_tvc = cif_cancermale_age0_tvc;
        age3_tvc = cif_cancermale_age3_tvc;
        death_cause = 'Cancer';
        output;
    end;

    keep _newt age0_tvc age3_tvc death_cause;
run;

data comb_tvc;
    set comb03 comb03_tvc;
run;
    
title 'comparison of effects of different age by cause interactions';
proc sgpanel data = comb_tvc;
    panelby death_cause/novarname;
    series x = _newt y = age0/legendlabel='<45';
    series x = _newt y = age3/legendlabel='75+';
    series x = _newt y = age0_tvc/legendlabel='<45' lineattrs=(pattern=dash);
    series x = _newt y = age3_tvc/legendlabel='75+' lineattrs=(pattern=dash);
    colaxis label = 'Time Since Diagnosis (years)';
    rowaxis values = (0 to .6 by .1) label= 'Probability of Death';
run;    


*   (m) separate knot placements
    m-(i) distribution of events;
title 'distribution of event times';
proc sgpanel data = _events_;
    where _death_ = 1;
    panelby status/novarname;
    histogram _t_;
	colaxis label = 'Time since diagnosis (years)';
run;

*   m-(ii) separate models;
*   no 'if' option in sas version of stpm2;
data _events_temp ; set _events_;run;
data _events_;set _events_temp (where = (cancer = 1));run;
title ;
%stpm2(fem_can age1can age2can age3can, scale=hazard, df=4, 
        dftvc=3, tvc=fem_can age1can age2can age3can , options = noprint);  
*   save knots points in macro strings from _model_ dataset;
proc sql noprint;
    select trim(desc) length = 64 into : knots_cancer from _model_ where comp = 'bh_knots';
    select trim(desc) length = 64 into : knots_cancer_tvc from _model_ where comp = 'knots_fem_can';
quit;

*	now create the _events_ dataset with just the non-cancer death records;

data _events_;set _events_temp (where = (other = 1));run;

%stpm2(fem_other age1oth age2oth age3oth, scale=hazard, df=4, options = noprint);  

proc sql noprint;
    select trim(desc) length = 64 into : knots_other from _model_ where comp = 'bh_knots';
quit;

*   m-(iii);
%put &=knots_cancer.;		*	to display the knot positions in the log file;
%put &=knots_other.;		
%put &=knots_cancer_tvc.;	

*	recover the full _events_ dataset;
data _events_;set _events_temp;run;

*	fit new model with cause-specific knot positions using the knotstvc option in %stpm2;
%stpm2 (cancer other fem_can fem_other 
    age1can age2can age3can age1oth age2oth age3oth , scale=hazard,
    options = rcsbaseoff noint eform,
    tvc = cancer other fem_can age1can age2can age3can,
    knotstvc = cancer: &knots_cancer.
    other: &knots_other. 
    fem_can: &knots_cancer_tvc.
    age1can: &knots_cancer_tvc.
    age2can: &knots_cancer_tvc.
    age3can: &knots_cancer_tvc.);

%stpm2cif (cancermale_age0_tvc2 othermale_age0_tvc2, cause1=cancer:1,
   cause2=other:1); 

%stpm2cif (cancermale_age3_tvc2 othermale_age3_tvc2, cause1=cancer:1 age3can:1,
   cause2=other:1 age3oth:1);  

*	combine CIF functions and rename for plotting later;
data comb03_tvc2;
	set _events_;

    if cif_othermale_age0_tvc2 ^= . then do;
        age0_tvc2 = cif_othermale_age0_tvc2;
        age3_tvc2 = cif_othermale_age3_tvc2;
        death_cause = 'Other Causes';
        output;
    end;
    
    if cif_cancermale_age0_tvc2 ^= . then do;
        age0_tvc2 = cif_cancermale_age0_tvc2;
        age3_tvc2 = cif_cancermale_age3_tvc2;
        death_cause = 'Cancer';
        output;
    end;

    keep _newt age0_tvc2 age3_tvc2 death_cause;
run;

*	combine with earlier version of CIF from the TVC model;
data comb_tvc2;
    set comb03_tvc comb03_tvc2;
run;

title 'comparison of tvc strategies by selected age group';
proc sgpanel data = comb_tvc2;
    panelby death_cause/novarname;
    series x = _newt y = age0_tvc2/legendlabel='<45' lineattrs=(pattern=dash);
    series x = _newt y = age3_tvc2/legendlabel='75+' lineattrs=(pattern=dash);
    series x = _newt y = age0_tvc/legendlabel='<45' ;
    series x = _newt y = age3_tvc/legendlabel='75+' ;
    colaxis label = 'Time Since Diagnosis (years)';
    rowaxis values = (0 to .6 by .1) label= 'Probability of Death';
run;    


*    (n) Now, estimate CIFs from a Cox model instead 
    - Assume that the effect of sex is the same on both outcomes;

*   use stratified cox model, which will allow CIF adjusted for covariates
	in much the same way as in the flexible parametric model;
*   expand data to one row for each patient and cause;
data colon_expanded;
    set colon;
    do cause = 1 to 2;
        cancer  = cause = 1;
        other   = cause = 2;
*   Event is the event indicator, coded like this:
    event is 1 if death from cause
    event is 0 if death from other;
         event=(cause=status);
         output;
    end;
    keep  event sex cause female surv;
run;

*   need a dataset to select which covariate patterns will appear in
    the cumulative hazard estimates;
data cov;
    sex = 1;
    output;
    sex = 2;
    output;
run;

*   stratified cox model (cause of death is stratum indicator) based on expanded data;
proc phreg data = colon_expanded;
    model surv*event(0) = sex/ties= breslow;
    strata cause;
    baseline covariates=cov out=strata_out cumhaz=ch ;
run;    

*   estimates are produced only at distinct failure times independently in each stratum;
data cif;
    merge strata_out    (rename = (ch = h_canc)     where = (cause = 1))
        strata_out      (rename = (ch = H_other)    where = (cause = 2));
        by sex surv;

    if first.sex then do;
        l_haz_c = h_canc;
        l_haz_oth = h_other;
        
        S0 = 1;
        cif_canc = 0;
        cif_oth = 0;
    end;
    
*   CIF calculation based on overall survival to beginning of interval,
    and cause-specific hazard within the interval.  Accounting for
    intervals where there has been no change in the hazard for that cause;
    else do;
    
*   hazard increments;
        canc_inc = ifn(h_canc =., 0, h_canc-l_haz_c);
        oth_inc  = ifn(h_other=., 0, h_other-l_haz_oth);
        
        cif_canc = cif_canc+S0*canc_inc;
        cif_oth = cif_oth+S0*oth_inc;

*   update survival estimate and save current hazard estimate for next interval;        
        S0 = S0*(1-(canc_inc+oth_inc));
        l_haz_c   = ifn(h_canc  = ., l_haz_c  , h_canc);
        l_haz_oth = ifn(h_other = ., l_haz_oth, h_other);
    end;
    
    retain l_haz_c l_haz_oth S0 cif_canc cif_oth;
    
    keep sex surv cif_canc cif_oth;
run;

title 'non-parametric CIF from stratified Cox model';
proc sgplot data = cif;
    series x = surv y = cif_canc/group = sex;
    series x = surv y = cif_oth/group = sex;
    xaxis label = 'TIme since diagnosis';
    yaxis label = 'Probability of death';
    format sex sex.;
run;

proc sgpanel data = cif;
    panelby sex/novarname;
    series x = surv y = cif_canc/legendlabel='Cancer';
    series x = surv y = cif_oth/legendlabel='Other Causes';
    rowaxis label = 'Probability of death';
    colaxis label = 'TIme since diagnosis';
    format sex sex.;
run;

*   (o) Relax the assumption that the effect of sex is the same on the two outcomes
*   model separate effcts of sex on cause of death;

data colon_expanded;
    set colon_expanded;
    
    fem_canc = female*(cause=1);
    fem_other = female*(cause=2);

*   variables will have values
    1       males, stratum of interest
    2       females, stratum of interest
    0       in other stratum;
    
    sex_canc = ifn(cause=1,sex,0);  
    sex_oth = ifn(cause=2,sex,0);   
run;

*   need all these covariate patterns for estimation of adjusted CIF;
data covm;      
    sex_canc = 1;
    sex_oth = 0;
    output;     
    sex_canc = 0;
    sex_oth = 1;
    output;     
    sex_canc = 2;
    sex_oth = 0;
    output;     
    sex_canc = 0;
    sex_oth = 2;
    output;     
run;

proc phreg data = colon_expanded;
    model surv*event(0) = sex_canc sex_oth;
    strata cause;
    baseline covariates=covm out=strata_outm 
        cumhaz=ch;
run;

*   calculations for males;
data cifm;
    merge strata_outm (rename = (ch = h_canc)   
            where = (cause = 1 and sex_canc = 1 and sex_oth = 0))
        strata_outm (rename = (ch = h_other)
            where = (cause = 2 and sex_canc = 0 and sex_oth = 1));
        by  surv;

    if _n_ = 1 then do;
        l_haz_c = h_canc;
        l_haz_oth = h_other;
        
        S0 = 1;
        cif_canc2 = 0;
        cif_oth2 = 0;
    end;
    
*   CIF calculation based on overall survival to beginning of interval,
    and cause-specific hazard within the interval.  Accounting for
    intervals where there has been no change in the hazard for that cause;
    else do;
    
*   hazard increments;
        canc_inc = ifn(h_canc =., 0, h_canc-l_haz_c);
        oth_inc  = ifn(h_other=., 0, h_other-l_haz_oth);
        
        cif_canc2 = cif_canc2+S0*canc_inc;
        cif_oth2 = cif_oth2+S0*oth_inc;

*   update survival estimate and save current hazard estimate for next interval;        
        S0 = S0*(1-(canc_inc+oth_inc));
        l_haz_c   = ifn(h_canc  = ., l_haz_c  , h_canc);
        l_haz_oth = ifn(h_other = ., l_haz_oth, h_other);
    end;
    
    retain l_haz_c l_haz_oth S0 cif_canc2 cif_oth2;
    
    keep sex_canc surv cif_canc2 cif_oth2 ;
run;
    
*   calculations for females;
data ciff;
    merge strata_outm (rename = (ch = h_canc)   
            where = (cause = 1 and sex_canc = 2 and sex_oth = 0))
        strata_outm (rename = (ch = h_other)
            where = (cause = 2 and sex_canc = 0 and sex_oth = 2));
        by  surv;

    if _n_ = 1 then do;
        l_haz_c = h_canc;
        l_haz_oth = h_other;
        
        S0 = 1;
        cif_canc2 = 0;
        cif_oth2 = 0;
    end;
    
*   CIF calculation based on overall survival to beginning of interval,
    and cause-specific hazard within the interval.  Accounting for
    intervals where there has been no change in the hazard for that cause;
    else do;
    
*   hazard increments;
        canc_inc = ifn(h_canc =., 0, h_canc-l_haz_c);
        oth_inc  = ifn(h_other=., 0, h_other-l_haz_oth);
        
        cif_canc2 = cif_canc2+S0*canc_inc;
        cif_oth2 = cif_oth2+S0*oth_inc;

*   update survival estimate and save current hazard estimate for next interval;        
        S0 = S0*(1-(canc_inc+oth_inc));
        l_haz_c   = ifn(h_canc  = ., l_haz_c  , h_canc);
        l_haz_oth = ifn(h_other = ., l_haz_oth, h_other);
    end;
    
    retain l_haz_c l_haz_oth S0 cif_canc2 cif_oth2;
    
    keep sex_canc surv cif_canc2 cif_oth2 ;
run;

*   combine adjusted and unadjusted CIF for plotting;
data cif_adj;
    set cifm (in=a)     /*  adjusted CIF for males      */
        ciff (in=b)     /*  adjusted CIF for females    */
        cif (in=c);     /*  unadjusted CIF, both sexes  */
    if a then sex = 1;
    else if b then sex = 2;
run;

title 'comparison of non-parametric CIF calculations';
title2 'separate effects of age on cause are dashed lines';
proc sgpanel data = cif_adj;
    panelby sex/novarname;
    series x = surv y = cif_canc/
        legendlabel='Cancer' lineattrs=(pattern=dash color=blue);
    series x = surv y = cif_canc2/
        legendlabel='Cancer'  lineattrs=(color=blue);
    
    series x = surv y = cif_oth/
        legendlabel='Other Causes'  lineattrs=(pattern=dash color=red);
    series x = surv y = cif_oth2/
        legendlabel='Other Causes'  lineattrs=(color=red);
        
    rowaxis label = 'Probability of death';
    colaxis label = 'TIme since diagnosis';
    format sex sex.;
run;

*   p) Test if the effect of sex differs between cancer and other causes;
proc phreg data = colon_expanded;
    model surv*event(0) = female|cause;
    strata cause;
run;

