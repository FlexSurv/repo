/*
	regression methods

	stpm2 macros 

    version:        3.0 May 2020

    macros required for the sas version of Paul Lamberts Stata program stpm2

    16 Oct 2015
    allow 'intercept-only' models.  That is, a model with no covariates.
    changes made to stpm2 and predict.

    1 - 18 April 2016
    PO models, theta model
    per option in predict for hazard and hazard difference estimates
    failure (= 1-S(t)) added as a measure for predict
    changed parameters for predict to accomodate new options (failure, cumhazard, sdiff1)

    6 June 2016
    added xbnobaseline, xb, dxb, meansurv to predict
    added normal scale for model
    
    13 Feb, 2017
    re-write to bring handling of user-supplied knots in line with stata
    rcsgen and stpm2
    
    18 March 2017 (released into the wild)
    added to cansurv repository.  renamed some modules, first draft of stpm2cif
    updated stom2cm.  lifelost added to predict
    
    20 March 2017   added grpd and by to stlifelost
                    (not sure that this actually works...)
    28 March 2017   added eform option to print 
                    exponentiated estimates from nlmixed
                    
    19 apr 2017     conthaz contmort hazard options for stpm2CIF    
    21 June 2017    use optsave/optload to restrict log output for both
                    stpm2 and predict

    7 July 2017     test for missing values in key variables and covariates in _events_ file
                    %stset now drops observations with time missing or <= 077
                    
    31 Oct 2017     vectorised processing for StLifeLost.  faster computation
                    of confidence intervals

    9 Nov 2017      fixed bug in hazard calc in stpm2CM
    
    16 Dec 2017     added check for non-numeric covariate

    28 Dec 2017     fixed bug in parsing specification of covariate-specific values for dftvc

    12 Jan 2018     implement lininit option to try if model does not converge
                    fits model with covariates and only the first of the requested
                    spline variables to obtain initial values for optimisation 
                    
    29 July 2018    added capability of merging on extra covariates to compute LifeLost
                    added check of uniqueness of popmort rows, given mergeby variables

    29 Aug 2018     added checks in StLifeLost
                    added ODS table to report on convergence status of model fit (_converge_)
                    so it can be reviewed in model building

	22 Jul 2019		now parses knotstvc properly.  If any tvc variables have custom knots, then ALL
					TVC variables must be so specified

	6 Aug 2019		attempts to detect linear dependencies in paramaters to be estimated.  If found,
					a report is given, and program is terminated

					estimation databases (_model_ _fit_ etc) are deleted at the beginning of a run
					if they exist from a previous run

	30 jan 2020		added reverse option to rcsgen
					added cure option to stpm2

	18 Feb 2020		added centile option for survival curves

	4 May 2020		fixed some bugs introduced by previous updates, specifically around prediction
					in models with TVC variables

	22 May 2020		re-write to bring the management of predictions more into line with Stata
					all predictions are now merged back on to the _events_ file

	21 Jan 2021		fixed bug in centile calculation when there are TVC variables

	25 Feb 2023		add error trap in IML calculations (error will still appear in log)
		
*/

    options NOQUOTELENMAX minoperator;
    
    %global save_knots;

/*
    stpm2

    fit survival data with restricted cubic splines for baseline cumulative log survival time

    parameters      usage
    ---------------------
    covar           list of covariates (numeric only)
                    this list can be empty (allow a null model)
    scale           scale for modeling (hazard odds theta normal)
    df              degrees for freedom for baseline splines 1 - 10 (3 - 11 for cure models)
    tvc             list of time varying covariates
    dftvc           degrees of freedom for splines describing time-varying covariates
    bknots          boundary knots (default depends on knscale)
    inittheta       = 1 (default) initial value for theta (used if scale is theta)
    constheta       value to be used if theta is to be considered fixed
    bhazard         name of mortality rate variable for excess hazard (relative survival) models
    knots           (internal) knot values to use instead of standard knot positions 
                    derived from df option
    knscale         scale that knot values are specified on
                    = centile 
                    = time (default)
                    = logtime
    knotstvc        specifies knot values to be used for time varying covariates (see knscale)
    bknotstvc       used to specify boundary knots for the tvc splines (otherwise,
                    boundary knots from bknots will be used)  * in development  *
    weights         =
                    name of variable holding weight for individual observations

    options         = list of options separated by blanks
					cure		to fit a cure model (requires scale = hazard
								and a bhazard variable, ie fit an excess hazard model)
                    eform       display exponentiated parameter estimates   
                    noprint     turn off printed reports
                    print       (default) print summary results from nlmixed and stset
                    noorthog    do not orthogonalise splines
                    orthog      (default) orthogonalise splines
                    noint       do not include a constant term in model
                    int         (default) include a constant term in model
                    rcsbaseoff  do not compute baseline splines
                    norcsbaseoff (default) compute baseline splines
                    nodebug     (default) turn off notes, verbose mode, macro print and logic, delete temporary files
                    debug1      turn on verbose mode
                    debug2      include debug1 and turn on notes
                    debug3      include debug2 and turn on macro print and logic, retain temporary files. 
                    lininit     use initial values from model with the 1st of the spline variables (the linear term)

*/
%macro stpm2(covar ,scale = hazard, df = , tvc =, dftvc =, bknots = ,
    bhazard = , knots = , knscale = time, knotstvc =, inittheta = 1,
    constheta = , bknotstvc=, 
    weights = , options = );

*   set option defaults;
    %let print = 1;
    %let orthog = 1;
    %let int = 1;
    %let rcsbaseoff = 0;
    %let debug = 0;
    %let verbose = 0;
    %let n_cov = 0;
    %let eform = 0;
    %let est = ;
    %let ods_nlmixed =  ;
    %let err = 0 ;
    %let lininit = 0;
    %let constraints = ;
	%let cure = 0;
	%let reverse = 0;
    
proc optsave out = work.myopts; run;    
options nonotes nomprint nomlogic ;

*   remove output datasets if they exist from a previous run;
    %if %sysfunc(exist(_model_)) %then %do;
        proc datasets library = work nolist nowarn force;
            delete _model_ _converge_ _fit_ _dim_ _cov_ _parms_;
        quit;
        run;
    %end;

*   evaluate options parameter;
    %if "&options." ne "" %then %do;
        %let options = %lowcase(&options.);
        %if %sysfunc(findw(&options., noprint)) ^= 0 %then %let print = 0;
        %if %sysfunc(findw(&options., eform)) ^= 0 %then %let eform = 1;
        %if %sysfunc(findw(&options., cure)) ^= 0 %then %do;
			%let cure = 1;
			%let reverse = 1;
			%let tvc_code = ;
			%let cure_model = 1;
		%end;
        %if %sysfunc(findw(&options., noorthog)) ^= 0 %then %let orthog = 0;
        %if %sysfunc(findw(&options., noint)) ^= 0 %then %let int = 0;
        %if %sysfunc(findw(&options., rcsbaseoff)) ^= 0 %then %let rcsbaseoff = 1;
        %if %sysfunc(findw(&options., debug1)) ^= 0 %then %let debug = 1;
        %if %sysfunc(findw(&options., debug2)) ^= 0 %then %let debug = 2;
        %if %sysfunc(findw(&options., debug3)) ^= 0 %then %let debug = 3;
        %if %sysfunc(findw(&options., lininit)) ^= 0 %then %let lininit = 1;
    %end;

    %if &debug. gt 0 %then %do;
        %if %sysfunc(findw(1 2 3,&debug.)) > 0 %then %let verbose = 1;
        %if %sysfunc(findw(2 3,&debug.)) > 0 %then %do;
            options notes;
        %end;
        %if &debug. eq 3 %then %do;
            options spool mprint mlogic;
        %end;
    %end;

*	conditions for cure models;
	%if &cure. %then %do;
		%if &bhazard. eq %then %do;
			%err_mess(bhazard parameter must be set for cure models); 
			%goto fin;
    	%end;
		%if &scale. ne hazard %then %do;
			%err_mess(scale must be 'hazard;' for cure models); 
			%goto fin;
    	%end;
		%if &df. lt 3 or &df. gt 11 %then %do;
			%err_mess(df must be between 3 and 11 for cure models); 
			%goto fin;
    	%end;
	%end;

*   check case file for specified variables;
    %if not %sysfunc(exist(_events_)) %then %do;
            %err_mess(_events_ file is missing.  Have you run stset?); 
            %goto fin;
        %end;
        
*   file exists.  check for key variables;
    %let dsid = %sysfunc(open(_events_));   
        
    %if not %sysfunc(varnum(&dsid.,_t_)) 
        %then %err_mess(time variable _t_ is not in _events_ file.  Have you run stset?); 

    %if not %sysfunc(varnum(&dsid.,_study_id_)) 
        %then %err_mess(unique ID variable _study_id_ is not in _events_ file.  Have you run stset?); 

    %if not %sysfunc(varnum(&dsid.,_death_)) 
        %then %err_mess(event status variable _death_ is not in _events_ file.  Have you run stset?); 

    %if not %sysfunc(varnum(&dsid.,cons)) 
        %then %err_mess(intercept variable cons is not in _events_ file.  Have you run stset?); 

*   test covariate list (if present) each variable must be in the _events_ file, and of numeric type ;
    %if &covar. ne %then %do;
        %do nc = 1 %to %sysfunc(countw(&covar.," "));
            %let varn = %scan(&covar.,&nc.);    
            %let vnum = %sysfunc(varnum(&dsid.,&varn.)) ;
            %if &vnum eq 0 %then %err_mess(covariate &varn. is not in _events_ file); 
            %else %if %sysfunc(vartype(&dsid., &vnum.)) ne N %then %err_mess(covariate &varn. is not numeric); 
        %end;
    %end;
    %let rc = %sysfunc(close(&dsid.));  
%if &err. %then %goto fin;

*   test for missing values in specified covariates;
    %if &covar. ne %then %do;
        %let miss_var =;
        %let cov_list = %sysfunc(compbl(&covar.));
        %let cov_list = %sysfunc(translate(&cov_list.,","," "));
        proc sql noprint;
            select nmiss(&cov_list.) into : miss_var
            from _events_;
        quit;
        %if &miss_var. ne 0 %then 
            %err_mess(_events_ file has missing values for some of the specified covariates);
    %end;

*   check for missing values in key variables (including base hazard rate, if specified);
    %let miss_var =;
    %if &bhazard. ne %then %let key_var = _t_, _study_id_, _death_, cons, &bhazard.;
    %else %let key_var = _t_, _study_id_, _death_, cons;
    proc sql noprint ;
        select cmiss(&key_var.) into : miss_var
        from _events_ ;
    quit;
    %if &miss_var. ne 0 %then 
        %err_mess(detected missing values for key variables);

%if &err. %then %goto fin;

    %let scale = %lowcase(&scale.);
    %if %sysfunc(findw(hazard odds theta normal, &scale.)) = 0 %then %do;
        %put &=scale.;
        %err_mess(Only the hazard, odds, normal or theta scales are allowed);
        %goto fin ;
    %end;

    %if "&df." ^= "" %then %do;
        %if %eval(&df.) < 1 or %eval(&df.) > 10 %then %do;
            %put &=df.;
            %err_mess(degrees of freedom for baseline must be between 1 and 10);
            %goto fin;
        %end;
    %end;

    %if "&df." = "" and "&knscale." ne "" %then %do;
        %if %sysfunc(findw(CEN TIM LOG,%substr(%upcase(&knscale.),1,3))) = 0 %then %do;
        %err_mess(%str(Knot scale specified was &knscale
        scale for baseline and tvc knots must be one of TIMe, CENtile or LOGtime));
        %goto fin;
        %end;
    %end;
    
    %if "&df." = "" and "&knscale." = "" and &rcsbaseoff. eq 0 %then %do;
        %err_mess(one of 'df' or 'knscale' must be specified if baseline splines are required);
        %goto fin;
    %end;


    %if &rcsbaseoff. eq 1 and ("&df." ne "" or "&knots." ne "" or "&bknots." ne "" ) %then %do;
        %put &=rcsbaseoff.  &=df.  KNOTS = &knots.  &=bknots.;
        %err_mess(%str(When baseline splines are not required, degrees of freedom,
        baseline knots and boundary knots must be null));
        %goto fin ;
    %end;

    %if "&df." ne "" and "&knots." ne "" %then %do;
        %put &=df.   &=knots.;
        %err_mess(provide only ONE of DF and KNOTS);
        %goto fin;
    %end;

    %if "&dftvc." ne "" and "&knotstvc." ne "" %then %do;
        %err_mess(provide only ONE of DFTVC and KNOTSTVC);
        %goto fin;
    %end;

    %if "&tvc." ne "" and "&dftvc." eq "" and "&knotstvc." eq "" %then %do;
        %put &=tvc   &=dftvc.  &=knotstvc.;
        %err_mess(If time-varying covariates are specified, you must supply one of DFTVC or knotstvc);
        %goto fin;
    %end;

    %if "&tvc." ne "" %then %do;
        %let n_tvc = %sysfunc(countw(&tvc.," "));
        %let all_found = 0;
        %do _itvc = 1 %to &n_tvc.;
            %let found = %sysfunc(findw(&covar.,%scan(&tvc.,&_itvc.)," ",E));
            %let all_found = %eval(&all_found. + (&found. > 0));
        %end;
        %if &all_found. ne &n_tvc. %then %do;
            %err_mess(All time-varying covariates must be declared in covariate list);
            %goto fin;
        %end;
    %end;

    %if "&knotstvc." ne "" %then %do;
        %let n_tvc = %sysfunc(countw(&tvc.," "));
        %let all_found = 0;
        %do _itvc = 1 %to &n_tvc.;
            %let found = %sysfunc(findw(&knotstvc.,%scan(&tvc.,&_itvc.)," :",E));
            %let all_found = %eval(&all_found. + (&found. > 0));
        %end;
        %if &all_found. ne &n_tvc. %then %do;
            %err_mess(if specifying knot positions for TVC variables, ALL TVCs must be specified);
            %goto fin;
        %end;
    %end;

    %if %str(&dftvc.) ne %then %do;
        %let i = 1;
        %do %while (%length(%scan(&dftvc.,&i.))>0);
            %if %index(%scan(&dftvc.,&i.," "),:) > 0 %then %do;
                %if %sysfunc(findw(&tvc., %scan(%scan(&dftvc.,&i.),1,:)," ",E)) = 0 %then %do;
                    %err_mess(variables in DFTVC must be declared in TVC list);
                    %goto fin;
                %end;
                %if %scan(%scan(&dftvc.,&i.),2,:) < 1 or  %scan(%scan(&dftvc.,&i.),2,:) > 9 %then %do;
                    %err_mess(degrees of freedom in DFTVC must be between 1 and 9);
                    %goto fin;
                %end;
            %end;
            %if %index(%scan(&dftvc.,&i.," "),:) = 0 %then %do;
                %if %scan(&dftvc.,&i.) < 1 or %scan(&dftvc.,&i.)  > 9 %then %do;
                    %err_mess(degrees of freedom in DFTVC must be between 1 and 10);
                    %goto fin;
                %end;
            %end;
        %let i = %eval(&i.+1);
        %end;
    %end;
    
*   process Boundary and user-supplied knots options;
%let u_knots =;
    %let i_knots =;         *   this will hold internal user knots transformed to log time;
    %if "&df." = "" and &rcsbaseoff. = 0 %then %do;
        %if "&bknots." = "" %then %do;
            proc sql noprint;
                select 
                    min(_ln_t) ,
                    max(_ln_t) 
                into    
                    :lowerknot, 
                    :upperknot 
                from _events_
                where _death_ = 1;
            quit;
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = t %then %do;
            %let lowerknot = %sysfunc(log(%scan(&bknots.,1)));
            %let upperknot = %sysfunc(log(%scan(&bknots.,2)));
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = l %then %do;
            %let lowerknot = %scan(&bknots.,1);
            %let upperknot = %scan(&bknots.,2);
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = c %then %do;
            %let pct1 = %scan(&bknots.,1);
            %let pct2 = %scan(&bknots.,2);
            proc univariate data = _events_
                (where = (_death_ = 1)) noprint;
                var _ln_t;
                output out = _knots_ PCTLPRE=P pctlpts = &pct1. &pct2.;
                %if "&weights." ^= "" %then weight &weights.%str(;);
            run;
            proc sql noprint;
                select 
                    P&pct1., 
                    P&pct2. 
                into 
                    :lowerknot, 
                    :upperknot 
                from _knots_ ;
            quit;
            proc datasets nolist ;
                delete _knots_;
            quit;
            run;    
        %end;

*   now process user-supplied (internal) knots; 
        %if %substr(%lowcase(&knscale.),1,1) = t %then 
            %let i_knots = %ln(&knots.) ;
        %else %if %substr(%lowcase(&knscale.),1,1) = l %then
            %let i_knots =  &knots. ;   
        %else %if %substr(%lowcase(&knscale.),1,1) = c %then %do;
            proc univariate data = _events_
                (where = (_death_ = 1)) noprint;
                var _ln_t;
                output out = _knots_ PCTLPRE=P pctlpts = &knots.;
                %if "&weights." ^= "" %then weight &weights.%str(;);
            run;
            proc transpose
                data = _knots_
                out = _t_knots_ ;
            run;
            proc sql noprint;
                select col1 format best12. into :i_knots
                separated by ' '
                from _t_knots_;
            quit;
            proc datasets nolist ;
                delete _knots_ _t_knots_;
            quit;
            run;    
        %end;
        %let u_knots = &lowerknot. &i_knots. &upperknot.;   
    %end;   

*   if tvc option specified, then check for tvc-specific boundary knots.  If none, use 
    those specified for the baseline splines;
    %if %str(&tvc.) ne %then %do;
        %if "&bknotstvc." eq "" %then %let bknotstvc = &bknots.;
        %if "&bknotstvc." = "" %then %do;
            proc sql noprint;
                select 
                    min(_ln_t) ,
                    max(_ln_t) 
                into    
                    :lowertvcknot, 
                    :uppertvcknot 
                from _events_
                where _death_ = 1;
            quit;
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = t %then %do;
            %let lowertvcknot = %sysfunc(log(%scan(&bknotstvc.,1)));
            %let uppertvcknot = %sysfunc(log(%scan(&bknotstvc.,2)));
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = l %then %do;
            %let lowertvcknot = %scan(&bknotstvc.,1);
            %let uppertvcknot = %scan(&bknotstvc.,2);
        %end;
        %else %if %substr(%lowcase(&knscale.),1,1) = c %then %do;
            %let pct1 = %scan(&bknotstvc.,1);
            %let pct2 = %scan(&bknotstvc.,2);
            proc univariate data = _events_
                (where = (_death_ = 1)) noprint;
                var _ln_t;
                output out = _knots_ PCTLPRE=P pctlpts = &pct1. &pct2.;
                %if "&weights." ^= "" %then weight &weights.%str(;);
            run;
            proc sql noprint;
                select 
                    P&pct1., 
                    P&pct2. 
                into 
                    :lowertvcknot, 
                    :uppertvcknot 
                from _knots_ ;
            quit;
            proc datasets nolist ;
                delete _knots_;
            quit;
            run;    
        %end;
    %end;
*   check for delayed-entry data (if _t0_ is present and at least one value is non-zero)
    if present, then compute indicator and log transform entry time;
    %let del_entry = 0;
    %let dsid = %sysfunc(open(work._events_));
    %if %sysfunc(varnum(&dsid.,_t0_)) ne 0 %then %do;
        %let rc = %sysfunc(close(&dsid.));
        proc sql noprint;
            select sum(_t0_) into :s 
            from _events_ where _t0_ >0;
        quit;   
        %if &s. > 0 %then %do;
            %let del_entry = 1;
            data _events_;
                set _events_;
                _de_ = (_t0_ > 0);
                if _de_ then _ln_t0 = log(_t0_);
                else _ln_t0 = 0;
            run;
        %end;
    %end;
    %else %let rc = %sysfunc(close(&dsid.));

*   get count of valid data;
    proc sql noprint;
        select count(*) into :Nt
        from _events_ where _t_ ^= .  and _t_ > 0;
    quit;

*   if fixed value for theta is supplied, check that scale is theta;
    %if "&constheta." ne "" %then %do;
        %let ln_theta = %sysfunc(log(&constheta.));
        %if &scale. eq theta %then %let constraints = %str(bounds &ln_theta. <= ln_theta <= &ln_theta.;);
        %else %do;
            %err_mess(constheta option is only valid for theta scale);
            %goto fin;
        %end;
    %end;

*   relative survival;
    %if %str(&bhazard.) ne %then %let relsurv = 1;
    %else %let relsurv = 0;

*   linear predictor and its derivative;
    %if &int. = 1 %then %do;
        %let xb =  cons*_cons;
        %let dxb = 0;
        %let parms = cons &covar.;
        %let xb0 = cons*_cons;
    %end;
    %else %do;
        %let xb = 0;
        %let dxb = 0;
        %let xb0 = 0;
        %let parms = &covar.;
    %end;

    %let del_var =;
    %let deriv = ;

*   process covariate list, if it is not null;
    %let n_cov = 0;
    %if &covar. ne  %then %do;
        %let n_cov = %sysfunc(countw(&covar.));
        %let i = 1;
        %do %while (&i <= &n_cov.);
            %let xb = &xb. + %scan(&covar.,&i.)*_%scan(&covar.,&i.);
            %if &del_entry. %then %let xb0 = &xb0. + %scan(&covar.,&i.)*_%scan(&covar.,&i.);
            %let i = %eval(&i.+1);
        %end;
    %end;

*   add basline spline terms if requested;
    %if &rcsbaseoff. = 0 %then %do;
        %if "&u_knots." = "" %then %let d = &df.;   
        %if %str(&df.) eq %then %let d = %eval(%sysfunc(countw(&u_knots., ' '))-1);
		%if &cure. %then %do;
			%if &constraints. eq %then %let constraints = bounds  0 <= rcs&d. <= 0;
			%else %let constraints = &constraints. , 0 <= rcs&d. <= 0;
		%end;
        %do i = 1 %to &d.;
            %let xb = &xb. + rcs&i.*_rcs&i.;
            %let dxb = &dxb. + rcs&i.*_drcs&i.;
            %let parms = &parms. rcs&i.;
            %let deriv = &deriv. drcs&i.;
            %if &del_entry. and &rcsbaseoff. = 0 %then %do;
                %let xb0 = &xb0. + rcs&i.*_s0_rcs&i.;
                %let del_var = &del_var. s0_rcs&i.;
            %end;
        %end;

*   generate baseline splines;
        %rcsgen(_ln_t, if1 = _death_ = 1, dgen = drcs, orthog = &orthog., df = &df., reverse = &reverse.,
            knots = &u_knots., fw = &weights.);

*   save the knots on the log time scale for later use by predict
    and its kin;
        %let ln_bhknots = &save_knots;
        %ms(ln_bhknots, %str(knots for baseline splines, log scale, with boundary knots));              

*   save internal knots on time scale.  trim boundary knots;
        %let nn = %sysfunc(countw(&save_knots.,' '));
        %let bh_knots = ;
        %if &nn > 2 %then %do;
            %do in = 2 %to %eval(&nn.-1);
                %let bh_knots = &bh_knots %scan(&save_knots.,&in.,' ');
            %end;
            %let bh_knots = %exp(&bh_knots.);
        %end;
        %ms(bh_knots, %str(internal knots for baseline splines, time scale));

*   if requested orthogonal, then save the T matrix;        
        %if &orthog. %then %do;
            proc datasets library = work nolist nowarn force;
                delete _T_bh_;
                change Tmat = _T_bh_;
            quit;
            run;
        %end;
    %end;

*   if delayed entry, then compute the baseline splines using the non-zero entry times;
    %if &del_entry. = 1 and &rcsbaseoff. = 0 %then %do;
        %rcsgen(_ln_t0,
            gen=s0_rcs, dgen = s0_drcs, tmatrix = _T_bh_, orthog = &orthog., reverse = &reverse.,
            knots = &save_knots., if2 = _ln_t0 ^= ., fw = &weights.);
    %end;


*   add list of time-varying parameters;
*   time-varing covariates, generate splines and create the interaction variables
    and first derivatives;

*   df for time-varying covariate splines;
    %let tvc_var =;
    %if %str(&dftvc.) ne %then %do;
        %let tvc_gen = 0;
        %let i = 1;
        %do %while (%length(%scan(&dftvc.,&i.)) > 0);
            %if %index(%scan(&dftvc.,&i.),:) = 0 %then %let tvc_gen = %scan(&dftvc.,&i.);
            %else %let tvc_var = &tvc_var %scan(&dftvc.,&i.) ;
            %let i = %eval(&i.+1);
        %end;
        %do i = 1 %to &n_tvc.;
            %if %sysfunc(findw(&dftvc.,%scan(&tvc.,&i.)," :",E)) = 0 %then
                %let tvc_var = &tvc_var. %scan(&tvc.,&i.):&tvc_gen.;
        %end;
    %end;

    %if "&tvc." ne "" %then %do;
        %do j = 1 %to &n_tvc.;
            %let _loc = %eval(2*&j.-1);
            %let tname =  %scan(&tvc_var.,&_loc., ' :');
/*          %let tname =  %scan(&tvc_var.,&j.);*/
            %let k =;
            %let t_knots =;
            %if %str(&dftvc.) ne %then %let d = %scan(%scan(&tvc_var.,&j.),2,:);
            %else %if "&knotstvc." ne "" %then %do;

  				%let tname = %scan(&tvc.,%eval(&j.));
				%let tvc_k_p = %index(&knotstvc.,&tname.);  
				%let tvc_k_s = %substr(&knotstvc., &tvc_k_p.);   
				%let k = %substr(&tvc_k_s.,%eval(%index(&tvc_k_s.,:)+1)); 
				%if %sysfunc(scan(&tvc_k_s.,-2,":")) = &tname. %then                     
					%let d = %eval(%sysfunc(countw(&k., ' '))+1); 

				%else %do;		/*	find next tvc knot specification, since this is not the last	*/
					%let next_tvc_p = %index(&k.,:); 		
					%let rev = %sysfunc(reverse(%substr(&k., 1, %index(&k.,:))));
					%let term = %sysfunc(anyspace(&rev.));		/*	end of the current tvc knot specification	*/
					%let k = %substr(&k.,1,%eval(%index(&k.,:)-&term.));   
					%let d = %eval(%sysfunc(countw(&k., ' '))+1);    
				%end;          
                  
                %let t_knots=;
                %if %substr(%lowcase(&knscale.),1,1) = t %then 
                    %let t_knots = &lowertvcknot. %ln(&k.) &uppertvcknot.;
                %else %if %substr(%lowcase(&knscale.),1,1) = l %then
                    %let t_knots = &lowertvcknot. &k. &uppertvcknot.;           
                %if %substr(%lowcase(&knscale.),1,1) = c %then %do;
                    proc univariate data = _events_
                        (where = (_death_ = 1)) noprint;
                        var _ln_t;
                        output out = _knots_ PCTLPRE=P pctlpts = &k.;
                        %if "&weights." ^= "" %then weight &weights.%str(;);
                    run;
                    proc transpose
                        data = _knots_
                        out = _t_knots_ ;
                    run;
                    proc sql noprint;
                        select col1 format best12. into :t_knots
                        separated by ' '
                        from _t_knots_;
                    quit;
                    proc datasets nolist ;
                        delete _knots_ _t_knots_;
                    quit;
                    run;    
                    %let t_knots = &lowertvcknot. &t_knots. &uppertvcknot.; 
                %end;   
            %end;
            %do i = 1 %to &d.;
                %let parms = &parms. rcs_&tname.&i.;
                %let deriv = &deriv drcs_&tname.&i.;
                %let xb = &xb + rcs_&tname.&i.*_rcs_&tname.&i.;
                %let dxb = &dxb + rcs_&tname.&i.*_drcs_&tname.&i.;
                %if &del_entry. %then %do;
                    %let del_var = &del_var. s0_rcs_&tname.&i.;
                    %let xb0 = &xb0. + rcs_&tname.&i.*_s0_rcs_&tname.&i.;
                %end;
            %end;
*	constraints for cure model, code to force initial values for nlmixed;
			%if &cure. %then %do;
				%let constraints = &constraints. , 0<= rcs_&tname.&d. <= 0;
				%let tvc_code = &tvc_code. %str(rcs_&tname.&d. = 0;);
			%end;
            %if "&t_knots." ne "" %then %let d1 = ;
            %else %let d1 = &d.;
            %rcsgen(_ln_t, if1 = _death_ = 1, reverse = &reverse.,
                gen= tvc, dgen= dtvc, orthog = &orthog.,
                df = &d1., knots = &t_knots.,  fw = &weights.);

*   save knots on log time scale for later use by predict and its kin;
        %let lnknots_&tname = &save_knots.;
        %ms(lnknots_&tname, %str(knots for this TVC, log scale, including boundary knots));

*   save internal knots on time scale;
        %let nn = %sysfunc(countw(&save_knots.,' '));
        %let knots_&tname = ;
        %if &nn > 2 %then %do;
            %let k = ;
            %do in = 2 %to %eval(&nn.-1);
                %let k = &k %scan(&save_knots.,&in.,' ');
            %end;
            %let knots_&tname = %exp(&k.);
        %end;
        %ms(knots_&tname, %str(internal knots for TVC splines, time scale));

            %if &del_entry. %then %do;
                %rcsgen(_ln_t0, reverse = &reverse.,
                    gen=s0_tvc, dgen=s0_dtvc, tmatrix = tmat, orthog = &orthog.,
                    knots = &save_knots., if2 = _ln_t0 ^= ., 
                    fw = &weights.);
            %end;

*   save the t matrix;
            %if &orthog. = 1 %then %do;
            proc datasets library = work nolist nowarn force;
                delete _T_&tname._;
                change tmat = _T_&tname._;
                quit;
            run;
            %end;

*   build the regressor variables;
            data _events_;
                set _events_;
                array new(*) rcs_&tname.1 - rcs_&tname.&d.;
                array dew(*) drcs_&tname.1 - drcs_&tname.&d.;
                array t(*) tvc1 - tvc&d.;
                array d(*) dtvc1 - dtvc&d.;

                do _i_ = 1 to &d.;
                    new(_i_) = t(_i_)*&tname.;
                    dew(_i_) = d(_i_)*&tname.;
                end;

                drop _i_ tvc1 - tvc&d. dtvc1 - dtvc&d.;
            run;

            %if &del_entry. = 1 %then %do;
                data _events_;
                    set _events_;
                    array new(*) s0_rcs_&tname.1 - s0_rcs_&tname.&d.;
                    array t(*) s0_tvc1 - s0_tvc&d.;

                    do _i_ = 1 to &d.;
                        new(_i_) = t(_i_)*&tname.;
                    end;

                    drop _i_ s0_tvc1 - s0_tvc&d.;
                run;
            %end;
        %end;
    %end;

*   build estimates specification if we want to see the exponentiated estimates
    from nlmixed;
    
    %if &eform %then %do;
        %let ods_nlmixed = %str(ods output AdditionalEstimates = _e_parms_;);
        %let n_parms = %sysfunc(countw(&parms));
        %do n_est = 1 %to &n_parms.;
            %let p = %scan(&parms., &n_est, ' ');
            %let est = &est. estimate "&p." exp(&p.)%str(;);
        %end;
    %end;

    %if &verbose. %then %do;
        %put  &=scale.;
        %if &rcsbaseoff. = 0 %then %put DF for baseline spline = &df.;
        %if "&knots" ne "" %then %do;
            %put Baseline knots = &bh_knots.;
            %put Knot Scale = &knscale.;
        %end;
        %if "&constheta." ne "" %then %put log of fixed value for theta = &ln_theta.;
        %if "&inittheta." ne "" %then %put initial value for theta = &inittheta.;
        %put  Covariate list = &covar.;
        %put  Parameters = &parms.;
        %put &=est.;
        %put  Linear predictor = &xb.;
        %put  Derivative = &dxb.;
        %put  Derivative Variables = &deriv.;
		%put &=constraints.;
        %put  Constant term required = &int.;
        %if &relsurv. %then %do;
            %put  Excess Hazard model;
            %put  Population mortality = &bhazard.;
        %end;
        %if %str(&tvc.) ne %then %do;
            %put  No. Time Dependent covariates = &n_tvc.;
            %put  Time Dependent covariates : df = &tvc_var.;
            %if "&knotstvc." ne "" %then %do;
                %put Supplied knots for TVC = &knotstvc.;
                %put Knot Scale = &knscale.;
            %end;
        %end;
        %put  Orthogonal splines = &orthog.;
        %if &del_entry. %then %do;
            %put Delayed Entry data;
            %put Linear Predictor for delayed entry = &xb0.;
            %put Delayed Entry variables = &del_var.;
        %end;
    %end;
    


%odsoff;

*	get sensible bounds for centile estimation in IML;
	proc means data = _events_ noprint ;
		where _death_ = 1;
		output out = _bounds_ p1(_t_) = p1 p99(_t_) = p99;
	run;
	data _bounds_; set _bounds_;
		call symput("p1", p1);
		call symput("p99", p99);
	run;

	%let bounds = &p1. &p99.;
	%ms(bounds, upper and lower bounds for centile estimation);

*   generate initial values;
proc phreg data = _events_ (keep = _study_id_ _t_ _death_  &parms. &bhazard. &weights.) ;
    model _t_*_death_(0) = &covar.;

    %if &weights. ^= %str() %then weight &weights.%str(;);
    output out = _cox1_ xbeta = xb;
run;



*   mean of cumulative sum of linear predictor;
proc means data = _cox1_ noprint;
    %if &weights. ^= %str()  %then weight &weights.%str(;);
    output out = _mxb_ (drop = _type_ _freq_)
        mean (xb) = mean_lp;
run;

*   center linear predictor;
data _cox2_;
    merge _cox1_
        _mxb_;
    if _n_ = 1 then mean = mean_lp;
    retain mean;
    xb = xb - mean;

    drop mean mean_lp;
run;

*   fit centered linear predictor and get cummulative baseline hazard funcion;
proc phreg data = _cox2_ noprint;
    model _t_*_death_(0) = xb;

    %if &weights. ^= %str()  %then weight &weights.%str(;);
    baseline out = _cum_ cumhaz = s;
    output out = _cox2out_ xbeta = xb2;
run;

*   censored events are dropped as a result of this join, as they are not informative for
    the linear regression of cumulative hazard on spline variables and covariates;
proc sql;
    create table _cox3_ as (select c1.*,
        c2.s,
        c3.xb2
        from _cox2_ c1, _cum_ c2, _cox2out_ c3
        where c1._t_ = c2._t_
            and c1._study_id_ = c3._study_id_
            and c1._death_ = 1);
quit;

data _cox4_;
    set _cox3_;
    if &relsurv. then s = s-.1*&bhazard.*_t_;
    s = exp(-s);
    Sadj = exp(xb2);
    Sadj = S**Sadj;
    if Sadj < 1 then do;
        if "&scale." = 'hazard' then            z = log(-log(Sadj));
        else if "&scale." = 'odds' then         z = log((1-Sadj)/Sadj);
        else if "&scale." = 'theta' then        z = log((Sadj**(-&inittheta.)-1)/(&inittheta.));
        else if "&scale." = 'normal' then       z = probit((&nt.*(1-Sadj)-3/8)/(&nt.+1/4));
    end;
    cons = 1;
run;

proc reg data = _cox4_ 
    outest = _s_ (drop =  _model_ _type_ _depvar_ _rmse_ z);
    model z = &parms. /noint;
    %if &weights. ^= %str()  %then weight &weights.%str(;);
    ods output ParameterEstimates = _regest;
run;

*	if this is a cure model, force the initial value for the last baseline parameter to be zero;
*	similarly for the last of any tvc splines;
%if &cure. and not &rcsbaseoff. %then %do;
	data _s_; 
		set _s_;
		rcs&df. = 0;

		if "&tvc." ^= "" then do;
			&tvc_code.;
		end;		
	run;
%end;

*	test for aliased covariates (linear dependency);
%let no_df =;
proc sql noprint;
	select variable into :no_df separated by ' '
		from _regest 
	where df = 0;
quit;

%if &no_df. ne %then %do;
	%err_mess(covariate(s) &no_df.  may be aliased.  look for linear dependencies); 
            %goto fin;
%end;

*   add the initial values for variates for first derivatives of baseline splines
    and any tvc variables;
%let baseon = ;
%let th_add = ;
%if &rcsbaseoff. = 0 %then 
    %let baseon= %str(array d1(*) _drcs1 - _drcs&d.;
        array dd2(*) _rcs1 - _rcs&d.;

        do _i_ = 1 to &d.;
            d1(_i_) = dd2(_i_);
        end;
        _type_ = 'parms';);

%if "&scale." = "theta" %then 
    %let th_add = %str(if "&constheta." = "" then ln_theta = log(input("&inittheta.",8.));
        else ln_theta = log(input("&constheta.",8.)););

    data _s_;
        set _s_;
    &baseon.;
    &th_add.;
        
    run;

*   add an underscore to each declared parameter variable name,
    so the variable name can be the representation of the parameter to be estimated in proc nlmixed;
    %add_h(_events_, &parms. );
    %add_h(_events_, &deriv.);
    %if &del_entry = 1 %then %do;
        %add_h(_events_, &del_var.);
    %end;
    
/*  processing of computed initial values   */  

*   log likelihood;
    %if &scale. eq hazard %then %do;
        %let st = exp(-exp(&xb.));
        %let ht = (&dxb.)*exp(&xb.);
        %if &relsurv. %then
            %let ht = &bhazard. + 1/_t_ *(&dxb.)*exp(&xb.);
        %else  %let ht = (&dxb.)*exp(&xb.);
    %end;
    %if &scale. eq odds %then %do;
        %let st = ((1 + exp(&xb.))**(-1));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_ *(&dxb.)*exp(&xb.)/((1 + exp(&xb.))));
        %else %let ht = ((&dxb.)*exp(&xb.)/(1 + exp(&xb.)));
    %end;

    %let th_parm = ;
    %if &scale. eq theta %then %do;
        %let st = (exp(ln_theta)*exp(&xb.) + 1)**(-1/exp(ln_theta));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_*((&dxb.)*exp(&xb.))/(exp(ln_theta)*exp(&xb.) + 1));
        %else %let ht = ((&dxb.)*exp(&xb.))/(exp(ln_theta)*exp(&xb.) + 1);
        %let th_parm =  ln_theta;
    %end;

    %if &scale. eq normal %then %do;
        %let st = cdf('normal',-(&xb.));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_*(&dxb.)*pdf('normal',&xb.)/cdf('normal',-(&xb.)));
        %else %let ht = ((&dxb.)*pdf('normal',(&xb.)))/cdf('normal',-(&xb.));
    %end;

    %let ll = _death_*log(&ht.) + log(&st.);

    %if &del_entry = 1 %then %let ll = &ll. + _de_*exp(&xb0.);

    %if "&weights" ne "" %then %let ll = &weights.*(&ll.);



/*  lininit option in force.  fit model with just the first of the rcs variables    */
*   have to collect:
    constant (if present)
    covariates 
    first baseline spline (if present)
    first delayed entry baseline spline (if present)
    first tvc spline (if any);
    
    %if &lininit. %then %do;

*       if constant;
        %if &int. = 1 %then %do;
            %let xb_l =  cons*_cons;
            %let dxb_l = 0;
            %let parms_l = cons &covar.;
            %let xb0_l = cons*_cons;
        %end;
        %else %do;
            %let xb_l = 0;
            %let dxb_l = 0;
            %let xb0_l = 0;
            %let parms_l = &covar.;
        %end;
        
*       add basline spline term if requested;
        %if &rcsbaseoff. = 0 %then %do;
            %let xb_l = &xb_l. + rcs1*_rcs1;
            %let dxb_l = &dxb_l. + rcs1*_drcs1;
            %let parms_l = &parms_l. rcs1;
            %if &del_entry. and &rcsbaseoff. = 0 %then %do;
                %let xb0_l = &xb0_l. + rcs1*_s0_rcs1;
            %end;
        %end;
        
*   process covariate list, if it is not null;
        %let n_cov = 0;
        %if &covar. ne  %then %do;
            %let n_cov = %sysfunc(countw(&covar.));
            %let i = 1;
            %do %while (&i <= &n_cov.);
                %let xb_l = &xb_l. + %scan(&covar.,&i.)*_%scan(&covar.,&i.);
                %if &del_entry. %then %let xb0_l = &xb0_l. + %scan(&covar.,&i.)*_%scan(&covar.,&i.);
                %let i = %eval(&i.+1);
            %end;
        %end;
        
*   time-varying covariate splines;
        %if "&tvc." ne "" %then %do;
            %do j = 1 %to &n_tvc.;
                %let _loc = %eval(2*&j.-1);
                %let tname =  %scan(&tvc_var.,&_loc., ' :');
                %let parms_l = &parms_l. rcs_&tname.1;
                %let xb_l = &xb_l. + rcs_&tname.1*_rcs_&tname.1;
                %let dxb_l = &dxb_l. + rcs_&tname.1*_drcs_&tname.1;
                %if &del_entry. %then %do;
                    %let xb0_l = &xb0_l. + rcs_&tname.1*_s0_rcs_&tname.1;
                %end;
            %end;
        %end;
        
    %if &scale. eq hazard %then %do;
        %let st = exp(-exp(&xb_l.));
        %let ht = (&dxb_l.)*exp(&xb_l.);
        %if &relsurv. %then
            %let ht = &bhazard. + 1/_t_ *(&dxb_l.)*exp(&xb_l.);
        %else  %let ht = (&dxb_l.)*exp(&xb_l.);
    %end;
    %if &scale. eq odds %then %do;
        %let st = ((1 + exp(&xb_l.))**(-1));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_ *(&dxb_l.)*exp(&xb_l.)/((1 + exp(&xb_l.))));
        %else %let ht = ((&dxb_l.)*exp(&xb_l.)/(1 + exp(&xb_l.)));
    %end;

    %let th_parm = ;
    %if &scale. eq theta %then %do;
        %let st = (exp(ln_theta)*exp(&xb_l.) + 1)**(-1/exp(ln_theta));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_*((&dxb_l.)*exp(&xb_l.))/(exp(ln_theta)*exp(&xb_l.) + 1));
        %else %let ht = ((&dxb_l.)*exp(&xb_l.))/(exp(ln_theta)*exp(&xb_l.) + 1);
        %let th_parm =  ln_theta;
    %end;

    %if &scale. eq normal %then %do;
        %let st = cdf('normal',-(&xb_l.));
        %if &relsurv. %then
            %let ht = (&bhazard. + 1/_t_*(&dxb_l.)*pdf('normal',&xb_l.)/cdf('normal',-(&xb_l.)));
        %else %let ht = ((&dxb_l.)*pdf('normal',(&xb_l.)))/cdf('normal',-(&xb_l.));
    %end;

    %let ll_lininit = _death_*log(&ht.) + log(&st.);

    %if &del_entry = 1 %then %let ll_lininit = &ll_lininit. + _de_*exp(&xb0_l.);

    %if "&weights" ne "" %then %let ll_lininit = &weights.*(&ll_lininit.);


proc nlmixed data = _events_ 
        (where = (_t_ ^=.))
        tech = newrap cov;
    &constraints.;;
    parms  / data = _s_ (keep = &parms_l. &th_parm.);
    model _t_ ~ general(&ll_lininit.);
    ods output ParameterEstimates = _parms_lininit;
run;

*   initial values from above fit.  All other initial values are defaulted to zero;
    
    proc transpose data = _parms_lininit 
        (rename = (parameter = _name_)) 
        out = _transposed_ name = parameter;
        var estimate;
    run;

*   set all initial parameter values to zero, then copy in the new values;
    data _s_;
        array p &parms.;
        do over p;
            p = 0;
        end;

            output;
        keep &parms.;
    run;
    
    data _s_;
        merge _s_ _transposed_;
    run;        

    %end;

%if &verbose. %then %put Log Likelihood = &ll.;

*   fit model, save parameters;
proc nlmixed data = _events_
        (where = (_t_ ^=.))
        tech = newrap cov;
    &constraints.;
    parms  / data = _s_ (keep = &parms. &th_parm.);
    model _t_ ~ general(&ll.);
    &est.;
    ods output ParameterEstimates = _parms_;
    &ods_nlmixed.;
    ods output FitStatistics = _fit_;
    ods output dimensions = _dim_;
    ods output CovMatParmEst = _cov_;
    ods output ConvergenceStatus = _converge_;
run;

*   print results of fit:  AIC, etc  observations used  parameter estimates;
%if &print. %then %do;
%odson;
proc print data = _fit_ noobs label;
    var descr value;
    label descr = 'Criterion'
        value = 'Value';
    title 'results of fit';
run;

proc print data = _dim_ noobs label;
    var descr value;
    label descr = 'Criterion'
        value = 'Value';
    title 'Description of data used';
run;


*   some condtional formatting and choice of what to print, depending on
    whether we have raw (log) estimates or exponentiated estimates;
%if &eform %then %do;
    %let p_list = _e_parms_;
    %let p_name = label;
    %let l_line = 'e(b)';
%end;
%else %do;
    %let p_list = _parms_;
    %let p_name = parameter;
    %let l_line = 'Estimate';
%end;
proc print data = &p_list. noobs label;
    var &p_name. estimate standarderror tvalue probt lower upper;
    label &p_name. = 'Parameter'
        estimate  = &l_line.
        standarderror = 'SE'
        tvalue = 't-value'
        probt = 'Prob t'
        lower = 'lower 95% ci'
        upper = 'upper 95% ci';
    title 'estimated parameters';
run;
%odsoff;
%end;

*   save current model definitions;
    %ms(scale,      scale in use for estimation);
    %ms(parms,      list of parameter names in log likelihood);
    %ms(int,        do we need an intercept);
    %ms(deriv,      list of variable names that are derivatives);
    %ms(covar,      list of covariates in model (not including constant));
    %ms(n_cov,      number of covariates in model);
    %ms(del_entry,  is there delayed entry data);
    %ms(orthog,     flag to request orthogonalisation of computed spline variables);
    %ms(relsurv,    is this excess hazard model);
    %ms(rcsbaseoff, turn off baseline splines);
    %if &constheta. ne %then
        %ms(constheta,  theta held at this value);
    %if &relsurv. = 1 %then
        %ms(bhazard,        name of population-based prob of survival for relative survival);
    %if %str(&tvc.) ne %then %do;
        %ms(tvc,        list of time-varying covariats);
        %ms(n_tvc,      number of time-varying covariates);
    %end;
	%if &cure. %then
		%ms(cure_model,	this is a cure model);

*   return the covariate names to their original form.  ;
    %rem_h(_events_, &parms.);
    %rem_h(_events_, &deriv.);
    %if &del_entry. = 1 %then %do;
        %rem_h(_events_, &del_var.);
    %end;

*   clean up;
    %if &debug. eq 0 %then %do;
        proc datasets library = work force nolist ;
        delete _start_ _cox1_ _cox2out_ _cox2_ _cox3_ 
        _cox4_ _cum_ _mxb_ _s_ _knots_  _regest
        _transposed_  _parms_lininit_ ;
        quit;
        run;
    %end;

*   error conditions branch to here;

%fin:
*   clean up;
    %if &debug. eq 0 %then %do;
        proc datasets library = work force nolist ;
        delete _start_ _cox1_ _cox2out_ _cox2_ _cox3_ 
        _cox4_ _cum_ _mxb_ _s_ _knots_  _regest
        _transposed_  _parms_lininit_ _bounds_;
        quit;
        run;
    %end;


proc optload data = work.myopts; run;
%odson;
%mend stpm2;


/*
    stpm2cif
    
    compute cumulative incidence functions after appropriate fit
    
    29 Jan 2017
    
*/
%macro stpm2cif( newvar, cause1 = , cause2 = , cause3 = , cause4 = ,
    cause5 = , cause6 = , cause7 =, cause8 =, cause9 =, cause10 =,
    maxt =, timename =_newt, obs = 1000, options = noci);
    
*   get model descripton from last fit and set local macros;

proc optsave out = work.myopts; run;
options nonotes nomprint nomlogic;

%let set = _events_;
%let ltime = _ln_t;
%let n_cov = 0;
%let xb = ;
%let dxb = ;
%let iml_vars = ;
%let ci = 0;
%let debug = 0;
%let verbose = 0;
%let contmort = 0;
%let conthaz = 0;
%let hazard = 0;

data _null_;
    set _model_;
    if comp ^= '' then call symput(comp, desc);
run;
    
*   unpack options paramter;
    %if "&options." ne "" %then %do;
        %let options = %lowcase(&options.);
        %if %sysfunc(findw(&options., ci)) >0 %then %let ci = 1;
        %if %sysfunc(findw(&options., contmort)) >0 %then %let contmort = 1;
        %if %sysfunc(findw(&options., hazard)) >0 %then %let hazard = 1;
        %if %sysfunc(findw(&options., conthaz)) >0 %then %let conthaz = 1;
        %if %sysfunc(findw(&options., debug1)) ^= 0 %then %let debug = 1;
        %if %sysfunc(findw(&options., debug2)) ^= 0 %then %let debug = 2;
        %if %sysfunc(findw(&options., debug3)) ^= 0 %then %let debug = 3;
    %end;

    %if &debug. gt 0 %then %do;
        %if %sysfunc(findw(1 2 3,&debug.)) > 0 %then %let verbose = 1;
        %if %sysfunc(findw(2 3,&debug.)) > 0 %then %do;
            options notes;
        %end;
        %if &debug. eq 3 %then %do;
            options spool mprint mlogic;
        %end;
    %end;   
    

        
*   count number of new variables to be created;
    %let _n_var = %sysfunc(countw(&newvar.," "));
    %if &_n_var. < 2 or &_n_var. > 10 %then %do;
        %err_mess(number of new variables must be between 2 and 10);
        %goto cif_fin;
    %end;   
    
*   at least the first 2 causes must be specified;
    %if "&cause1." = "" or "&cause2." = "" %then %do;
        %err_mess(cause1 and cause2 MUST be specified);
        %goto cif_fin;
    %end;

*   count total number of causes specified;
    %let _n_cause = 2;
    %do _nc= 3 %to 10 ;
        %if &&cause&_nc. ^=  %then %let _n_cause = %eval(&_n_cause. + 1);
    %end;
    
*   number of new variables should match the number of causes;
    %if &_n_var. ^= &_n_cause. %then %do;
        %err_mess(number of new variables must match number of causes specified);
        %goto cif_fin;
    %end;
    
*   set time limits for time points to be generated;
    %if "&maxt." = "" %then %do;
        proc sql noprint;
            select 
            max(_t_) into :maxt from _events_
			where _death_ = 1;
        quit;
    %end;
    %let step = %sysevalf(&maxt./&obs.);
    
*   create dataset for estimation;  
    %let set = _events_t_;
    %let ltime = l_&timename.;
    data _events_t_;
        set _events_;
        if _n_ <= &obs. then do;
            &timename. = _n_*&step.;
            &ltime. = log(&timename.);
            output;
        end;
    run;

*   get prediction functions for each survival and subhazard function;
    %let s_all = 1;
    %let measure = cif;
*   first, the survival all prediction function;
    %let timevar = &timename.;
    %do i1 = 1 %to &_n_cause.;
        %let at = &&cause&i1. zero;
        %if &i1. = 1 %then %p_iml(&at.,1);
        %else %p_iml(&at.,0);
        %let s_all = (&s_all.)*(exp(-exp(&xb.)));
    %end;
    
*   now, the individual integrand functions;
    %do i2 = 1 %to &_n_cause.;
        %let at = &&cause&i2. zero;
        %p_iml(&at.,0);
        %let pred&i2. = (&s_all.)*(exp(-log(x[1]) + log(&dxb.) + (&xb.)));
/*      %let p_haz&i2. = (exp(-log(x[1]) + log(&dxb.) + (&xb.))); */
        %let p_haz&i2. = (1/x[1])*(&dxb.)*exp(&xb.);
    %end;

*   generate spline variables and any tvc variables;
    %if &rcsbaseoff. = 0 %then %do;
        %rcsgen(&ltime., dgen=drcs, orthog = &orthog., set = _events_t_,    
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

    %if %str(&tvc) ne %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(&ltime., gen = rcs_&t., dgen = drcs_&t., knots = %superq(lnknots_&t.),  
                set = _events_t_, orthog = &orthog.,  tmatrix = _T_&t._, if2 = &ltime ^= . );
        %end;
    %end;
    
    %if &verbose %then %do;
        %put IML vars = &iml_vars.;
        %put;
        %put parameters = &parms.;
        %put;
        %put &=pred1;   
        %put &=pred2;   
        %put &=p_haz1;  
        %put &=p_haz2;  
    %end;

    %do i_cause = 1 %to &_n_cause.;
        %let pf = &&pred&i_cause.;
        %let ph = &&p_haz&i_cause.;
        %let varn = %scan(&newvar., &i_cause.); 
    
        proc iml;
*   variance - covariance matrix from last model fit;
        use _cov_;
        read all var {&parms.} into cov;

*   parameter estimates from model fit;
        use _parms_;
        read all var {estimate} into est;

*   data points (and required variables) at which to estimate function of interest;
        use &set.;
        read all var {&iml_vars.} into obs;

*   maxtrix to hold predicted values, with elements time, cif, lcl, ucl, SE 
    and hazard prediction
    requesting confidence intervals increases the processing time;

    varnames = {"&timename." "cif_&varn."};
    n_col = 2;
    if "&ci." = "1" then do;
        varnames = varnames ||{"cif_&varn._lci" "cif_&varn._uci" "se_cif&varn."};   
        n_col = n_col+3;
    end;

    if "&hazard." = "1" | "&conthaz." = "1" then do;    
        n_col = n_col+1;
        varnames = varnames || "h_&varn.";
    end;

    pred = j(nrow(obs),n_col,0);

*   temporary vector of computed interval-specific predictions;
        Pr = j(nrow(obs),1,0);

*   matrix to hold observation-specific derivatives;
*       grad = j(nrow(obs), nrow(est), 0);

*   lower triangular matrix with interval width
    half interval width on diagonal;
        L = j(nrow(obs),nrow(obs),0);
            do i = 1 to nrow(obs);
                do j = 1 to i;
                    L[i,j] = ifn(i=j,&step./2,&step.);
                end;
            end;
        L[1,1] = &step.;

*   functions of parameters to estimate
        x is vector of observed values
        e is vector of parameter estimates;

*	added code to trap errors;
        start pred(e) global (x);
			errFlag = 1;           /* set flag. Assume we will get an error :-) */   
			on_error = "if errFlag then do;  p = .;  resume; end;";
			call push(on_error);   /* PUSH code that will be executed if an error occurs */
            p = &pf.;
			errFlag = 0;
            return (p);
        finish pred;

        start phaz(e) global (x);
            ph = &ph.;
            return (ph);
        finish phaz;

*   call finite difference solver to get predictions, 1st derivatives of functions
    wrt parameter estimates for each point of interest;
        do i = 1 to nrow(obs);
            x = obs[i,];    

    if "&ci." = "0" then do;                            
            pred[i,1] = obs[i,1];                   *   time variable;
            Pr[i]     = pred(est);                  *   prediction;
    end;

*   cif estimate if confidence intervals requested;
    if "&ci" = "1" then do;
        pred[i,1] = obs[i,1];                   *   time variable;
        call nlpfdd(p, grd, h, "pred", est);    *   returns estimate (p) and 1st derivatives (grd);
        Pr[i]     = p;                          *   prediction;
*       grad[i,] = grd;                         *   observation-specific derivatives;
    end;
            
*   hazard estimate;
        if "&hazard." = "1" | "&conthaz." = "1"
            then pred[i,n_col] = phaz(est);
    end;    

*   matrix calculations to generate estimate and SE of estimate, based on delta method
    (if requested);
    if "&ci."= "1" then do;
        pred[,2] = L*Pr;
        pred[,5] = sqrt(vecdiag(L*grd*cov*t(grd)*t(L)));
        pred[,3] = pred[,2] - 1.96*pred[,5];
        pred[,4] = pred[,2] + 1.96*pred[,5];
    end;
    
    if "&ci."= "0" then pred[,2] = L*Pr;

*   create sas dataset for further manipulation;
            
        create _pred_ from pred[colname = varnames];
        append from pred;
    quit;
        
*   set up data step commands to compute summary measures (conthaz, contmort)
    and to retain hazard estimates, if required by the hazard option;   
        %if &conthaz. %then %do;
            %let h1 = tot_haz = h_&varn.;
            %let h2 = tot_haz = tot_haz + h_&varn.;
        %end;
        %else %do;
            %let h1 =;
            %let h2 =;
        %end;

        %if &contmort. %then %do;
            %let c1 = tot_cif = cif_&varn.%str(;);
            %let c2 = tot_cif = tot_cif + cif_&varn.%str(;);
        %end;   
        %else %do;
            %let c1 =;
            %let c2 =;
        %end;
        
        %let s1 =;
        %let s2 =;
    
        %if &i_cause. = 1 %then %do;
            data _events_;
				merge _events_ _pred_;
                &h1.; &c1.; &s1.
            run;
        %end;
        %else %do;
            data _events_;
				merge _events_ _pred_;
                &h2.; &c2.; &s2.
            run;
        %end;
    %end;

*   processing for hazard, conthazard and contmort options;

    %if &conthaz. or &contmort. %then %do;

        %let h3 = ;
        %let d3 =;
        %if &conthaz. %then %do;
            %if not &hazard. %then %let d3 = drop;      
            %do i_cause = 1 %to &_n_cause.; 
                %let h3 = &h3. conthaz_%scan(&newvar., &i_cause.) = h_%scan(&newvar., &i_cause.)/tot_haz %str(;);
                %if not &hazard. %then %let d3 = &d3. h_%scan(&newvar., &i_cause.);
            %end;
            %let h3 = &h3. drop tot_haz %str(;);
            %if not &hazard. %then %let d3 =&d3. %str(;);
        %end;
        %let c3 = ;     
        %if &contmort. %then %do;
            %do i_cause = 1 %to &_n_cause.; 
                %let c3 = &c3. contmort_%scan(&newvar., &i_cause.) = cif_%scan(&newvar., &i_cause.)/tot_cif %str(;);
            %end;
            %let c3 = &c3. drop tot_cif %str(;);
        %end;
            
        data _events_;
			merge _events_ _pred_;
                &c3.;
                &h3.;
                &d3.;
        run;
    %end;   
        
*   clean up;
    %if &debug. eq 0 %then %do;
        proc datasets library = work force nolist ;
        delete _pred_ ;
        quit;
        run;
    %end;

*   error exit;
    %cif_fin:
proc optload data = work.myopts; run;   
%mend stpm2cif;

/*
    stpm2CM
    compute crude probabilities of mortality (due to cancer/ other causes) after
    a relative survival regression fit using %stpm2

	added conditional survival 
	CIF after an initial survival time given in the paramter tcond (default = 0)

	July 2020

    parmaters:
    var         the prefix to be used for the returned variables 
                variables returned in cif_dat :
                the time variable created = var_t if tgen parameter not specified
                var_c
                var_o
                var_all
                var_s_all
                var_lamda
                var_rate
                var_St_star
                plus the following if the ci option has been chosen:
                var_c_lci
                var_c_uci
                var_c_se
                var_o_lci
                var_o_uci
                var_o_se
    at          covariate pattern (if any) for calculations
    using       population life table (default = popmort)
    popprob     variable in using dataset that holds the death probability (default = prob)
    mergeby     variables on the using dataset to merge to case data, in the following
                order: age, sex, year  (default = _age sex _year)
    sex         estimation is made for subjects with this value (default = 1)
    diagyear    year of diagnosis (a value, default = 2000)
    diagage     age at diagnosis (a value, default = 50)
    maxage      maximum age in the popmort file (default = 99)
    maxyear     maximum year available in popmort file (default = 2015)
	tcond		conditional survival time (a value, default = 0)
    nobs        number of time points to use in estimation (default = 1000)
    maxt        largest time point (years after diagnosis) to use (default = 10)
    tgen        time variable to be generated
    mergegen    other variables (and values) to be used in the merge to the 
                popmort file (eg. geographic region, SES, etc) 
    options     ci / noci  (default = noci )  compute confidence intervals using delta method
                nodebug     (default) turn off notes, verbose mode, macro print and logic
                debug1      turn on verbose mode
                debug2      include debug1 and turn on notes
                debug3      include debug2 and turn on macro print and logic
    
*/

%macro stpm2CM(var , at =,  popprob = prob, nobs = 1000, maxt = 10,  maxage = 99, tgen = ,
    using = popmort, mergeby = _age sex _year, sex = 1, diagyear = 2000,  diagage = 50,
    mergegen = , maxyear = 2015, tcond = 0, options = noci nodebug );

proc optsave out = work.myopts; run;
options nonotes nomprint nomlogic;

*   get model descripton from last fit and set local macros;
%let tvc = ;
%let iml_vars = ;
%let xb = ;
%let dxb = ;
%let ci = 0;
%let debug = 0;
%let err = 0;
%let measure = hazard;
%let verbose = 0;


data _null_;
    set _model_;
    if comp ^= '' then call symput(comp, desc);
run;

*   evaluate options parameter;
    %if "&options." ne "" %then %do;
        %let options = %lowcase(&options.);
        %if %sysfunc(findw(&options., ci)) >0 %then %let ci = 1;
        %if %sysfunc(findw(&options., debug1)) ^= 0 %then %let debug = 1;
        %if %sysfunc(findw(&options., debug2)) ^= 0 %then %let debug = 2;
        %if %sysfunc(findw(&options., debug3)) ^= 0 %then %let debug = 3;
    %end;

    %if &debug. gt 0 %then %do;
        %if %sysfunc(findw(1 2 3,&debug.)) > 0 %then %let verbose = 1;
        %if %sysfunc(findw(2 3,&debug.)) > 0 %then %do;
            options notes;
        %end;
        %if &debug. eq 3 %then %do;
            options spool mprint mlogic;
        %end;
    %end;

*   check for relative survival fit;
    %if &relsurv ne 1 %then %err_mess(stpm2CM can only be used after a fit with the bhazard option set (excess hazard model));

*   some checking of population life table;
    %if "&using." eq "" %then %err_mess(population life table parameter is required); 
    %if not %sysfunc(exist(&using.)) %then %err_mess(life table file &using. not found); 

*   must be at least 3 variables in the mergeby parameter;
    %let wc = %sysfunc(countw(&mergeby.,' '));
    %if &wc. lt 3 %then %err_mess(at least 3 variables must be specified (age/sex/year) in mergeby parameter ); 

*   are all variables present in the life table?;
    %let dsid=%sysfunc(open(&using.));
    %if %str(&mergeby.) ne %then %do;
        %let wc = %sysfunc(countw(&mergeby.,' '));  
        %do c = 1 %to &wc.;
            %if not %sysfunc(varnum(&dsid., %scan(&mergeby.,&c.))) %then 
                %err_mess(variable %scan(&mergeby.,&c.) is not in &popmort.);
            %else %let merge&c. = %scan(&mergeby.,&c.);
        %end;
    %end;
    
    %if not %sysfunc(varnum(&dsid., &popprob.)) %then %err_mess(variable &popprob. is not in &using. file); 
    %let rc = %sysfunc(close(&dsid.));
    
*	check for non-negative conditioning time;
	%if &tcond. lt 0 %then %err_mess(condition time must be non-negative (you chose tcond = &tcond.));		

*   check for valid maxt and nobs;
    %if &maxt. < 0 %then %err_mess(stpm2CM requires a positive maxt option (you chose maxt = &maxt.));
    %if &nobs. < 0 %then %err_mess(stpm2CM requires a positive Nobs option (you chose nobs = &nobs.));

    %if &err %then %goto cm_fin;
    
*   build data step code for fixed_sex, as it must be able to cope with both character and numeric coding;
    %let m2_code = fixed_sex = input("&sex.",1.);
    %if %sysfunc(anyalpha(&sex.)) gt 0 %then %let m2_code = fixed_sex = "&sex.";

*   set name for time variable to be generated;
    %if "&tgen." = "" %then %let timevar = &var._t; 
    %else %let timevar = &tgen.;

*   check for list of extra merge variables (and values to be used);
    %let mg_extra = ;
    %let code_extra = ;
    %if "&mergegen." ne "" %then %do;
        %let n_mergegen = %sysfunc(countw(&mergegen.));
        %do _in = 1 %to &n_mergegen;
            %let w = %scan(&mergegen.,&_in," ");
            %let ex_var = %scan(&w.,1,':');
            %let val = %scan(&w.,2,':');
            %let mg_extra = &mg_extra. and a.%scan(&w.,1,':') = b.%scan(&w.,1,':');
            %let code_extra = &code_extra. %str(&ex_var. = &val.;);
        %end;
    %end;
    
%if &verbose. %then %do;
%put *** code resulting from mergegen paramater;    
    %put &=mg_extra.;
    %put &=code_extra.;
%end;
    
    
*   build dataset for alternate time variable analysis.  assumed constant width is taken from last interval;
    %let ltime = l_&timevar.;
    %let width = %sysevalf(&maxt./&nobs.);
    data _events_t_;
        set _events_;
        if _n_ <= &nobs.;
        &timevar. = _n_*&width.;
        start = (_n_-1)*&width.;        *   start of interval, as in strs;
        l_&timevar. = log(&timevar.);
        fixed_age = input("&diagage.",4.);
        fixed_year = input("&diagyear.",4.);
        &m2_code.;
/*      if anyalpha("&sex.") then fixed_sex = "&sex."; */
/*      else fixed_sex = input("&sex.",1.);  */
        &code_extra.;   
    run;

*   get expected survival functions
        St_star         expected survival to time t
        rate            expected hazard at time t (aka lambda_star);
   	 proc sql noprint;
   	 create table _exp_ as select 
        a.&timevar. ,
        b.&popprob.**&width.        as p_star,
        -1*log(b.&popprob.)         as rate
        from _events_t_ a left join &using. b
            on min(floor(a.fixed_age + start), &maxage.) = b.&merge1.
            and a.fixed_sex = b.&merge2.
            and min(floor(a.fixed_year + start), &maxyear.) = b.&merge3.
            &mg_extra.
    	order by a.&timevar. ;

    select count(*) into : missed from _exp_ where p_star = .;
    quit;   

    %if &missed. > 0 %then %err_mess(&missed. records from _events_t_ file did not link to &using. database. PredictCM processing stops);
    %if &err. %then %goto cm_fin; 

    data _exp_;
        set _exp_;

        if _n_ = 1 then do;
			St_star = 1;
			st_star_cond = 0;
			rate_cond = 0;
		end;

        retain  St_star st_star_cond rate_cond;

        St_star = St_star*p_star;
        rate = -log(p_star)/&width.;

*	if this is conditional processing, then compute expected survival and rate at
	the condition time;
		if &tcond. > 0 then do;
			if &timevar. <= &tcond. then do
				st_star_cond = st_star;
				rate_cond = rate;
			end;
		end;

        drop p_star;
    run;

*   get components of prediction functions and variable list;

    %p_iml(&at.,1);

*   for CIF ca;
        %let pred_ca = (1/x[1]) * (&dxb.) * exp(&xb. -exp(&xb.))*y[1];              

*   for CIF other;
        %let pred_oth = exp(-exp(&xb.))*y[1]*y[2];

*   for overall survival;
        %let pred_s_all = exp(-exp(&xb.))*y[1];

*   for hazard (lambda);
        %let pred_lambda =  (1/x[1]) * (&dxb.) * exp(&xb.);                     
    
*   generate spline variables and any tvc variables;
    %if &rcsbaseoff. = 0 %then %do;
        %rcsgen(&ltime., dgen=drcs, orthog = &orthog., set = _events_t_,    
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

    %if %str(&tvc) ne %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(lnknots_&t.)," "))-1);
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(&ltime., gen = rcs_&t., dgen = drcs_&t., knots = %superq(lnknots_&t.),  
                set = _events_t_, orthog = &orthog.,  tmatrix = _T_&t._ );
        %end;
    %end;

    %if &verbose. %then %do;
        %put *** local macros in use ***;
        %put &=pred_ca;
        %put &=pred_oth;
        %put &=ci;
        %put &=width;
    %end;

%if &tcond. gt 0 %then %do;
*	numerator of conditional measure;	
	%let xb_n = &xb;				
	%let dxb_n = &dxb.;		
	
*	get the denominator (conditioning time)
	returns description of xb and dxb with references to the s[] matrix loaded
	from temp_one with spline values for the condition time tcond;
	%p_iml(&at.,0, cond = 1);	

*   for CIF ca;
        %let pred_ca = (1/x[1]) * (&dxb_n.) * exp(&xb_n. -exp(&xb_n.))*y[1]/(yc[1]*exp(-exp(&xb.)));              

*   for CIF other (y[1] is cumulative expected survival to time t, y[2] is the rate at time t
	yc[] holds the same values for t = tcond;
        %let pred_oth = exp(exp(&xb.) -exp(&xb_n.))*y[1]*y[2]/yc[1];

*   for overall survival;
        %let pred_s_all = exp(exp(&xb.)-exp(&xb_n.))*y[1]/yc[1];

*   for hazard (lambda);
        %let pred_lambda =  (1/x[1]) * (&dxb_n.) * exp(&xb_n.);                     

*	build dataset to provide condition;
	data _temp_one_;
		set _events_ (obs = 1);
		&timevar. = &tcond.;
		l_&timevar. = log(&tcond.);
		call symput("l_cond", l_&timevar.);
	run;

*   generate spline variables and any tvc variables;
    %if &rcsbaseoff. = 0 %then %do;
        %rcsgen(l_&timevar.,  dgen=drcs, orthog = &orthog., set = _temp_one_, scalar = &l_cond.,
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

    %if %str(&tvc) ne %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(l_&timevar., gen = rcs_&t., dgen = drcs_&t., knots = %superq(lnknots_&t.),  
                scalar = &l_cond.,set = _temp_one_, orthog = &orthog.,  tmatrix = _T_&t._ );
        %end;
    %end;
%end;


*   iml code cloned from mainline predict, and tooled for this particular situation
    two distinct prediction functions are being evaluated;
    proc iml;

*   variance - covariance matrix from last model fit;
    use _cov_;
    read all var {&parms.} into cov;

*   parameter estimates from model fit;
    use _parms_;
    read all var {estimate} into est;
    
*   data points (and required variables) at which to estimate function of interest;
    use _events_t_;
    read all var {&iml_vars. } into obs;

	read all var {_study_id_} into ID;

*   crude probabilities.  read in expected values;
    use _exp_; read all var{St_star rate} into exp;

*	if conditional estimates then bring in the conditioning data, otherwise
	just create a dummy dataset with one row.  It will only be referenced
	in the prediction function if this is actually a conditional estimation
	similarly, create dummy yc[] matrix and load it with expected survival at
	time = tcond and rate at same time;
	s = obs[1,1:ncol(obs)];
	yc = j(1,2,0);
	if &tcond. > 0 then do;
		use _temp_one_;
		read  var {&iml_vars.} into s;

		use _exp_ where (&timevar. <= &tcond.);
		read all var {st_star_cond  &timevar.} into ytemp;

*	the last value (closest to but not more than the condition time) is the one we want;
		yc = ytemp[nrow(ytemp),];
	end;


*   maxtrix to hold predicted values, with elements time, CIF for cancer, lcl, ucl, se, CIF for other causes lcl, ucl, se
        or just the time and two CIF values, if confidence intervals were not requested
		add a column for study ID, so we can merge back with the _events_ file;
    if &ci. then pred = j(nrow(obs),14,0);
    else pred = j(nrow(obs),8,0);

*   temporary vector of computed interval-specific predictions;
    Pr = j(nrow(obs),1,0);

*   matrix to hold observation-specific derivatives;
    grad = j(nrow(obs), nrow(est), 0);

*   lower triangular matrix with interval width
    half interval width on diagonal
    first column also 1/2 width ([1,1] element is 1/4 width);
    L = j(nrow(obs),nrow(obs),0);
            do i = 1 to nrow(obs);
                do j = 1 to i;
                    L[i,j] = ifn(i=j,&width./2,&width.);
                end;
            end;
    L[,1] = L[,1]/2;

    *   function of parameters to estimate 
        x is vector of observed values
        y is a vector holding constants related to expected rates
        e is vector of parameter estimates;

    *   the preduction measure for cancer cause of death;
	*	calls pred_ca ;
    start pred_ca(e) global (x, y, s, yc);
		if x[1] < &tcond. then p = 0;
		else p = &pred_ca.;
        return (p);
    finish pred_ca;

    *   the preduction measure for other causes of death;
    start pred_oth(e) global (x, y, s, yc);        
		if x[1] < &tcond. then p = 0;
        else p = &pred_oth.;
        return (p);
    finish pred_oth;
    
    *   prediction function for overall survival;
    start pred_s_all(e) global (x,y, s, yc);
		if x[1] < &tcond. then p = 0;
        else p = &pred_s_all.;
        return p;
    finish pred_s_all;

    *   prediction function for hazard;
    start pred_lambda(e) global (x,y);
        p = &pred_lambda.;
        return p;
    finish pred_lambda;

*   cancer causes;
    do i= 1 to nrow(obs);
        x = obs[i,];
        y = exp[i,];

        if &ci. then do;
            call nlpfdd(p, grd, h, "pred_ca", est); 
            Pr[i]     = p;  
            grad[i,] = grd; 
        end;
        else Pr[i] = pred_ca(est);
    end;


    pred[,1] = obs[,1];         *   time variable;
    pred[,2] = L*Pr;            *   CIF for cancer;
    pred[,4] = pred[,2];        *   will hold CIF for total;

*   confidence intervals for cancer CIF, if requested;
    if &ci. then do;
        pred[,11] = sqrt(vecdiag(L*grad*cov*t(grad)*t(L)));
        pred[,9] = pred[,2] + 1.96*pred[,11];
        pred[,10] = pred[,2] - 1.96*pred[,11];
    end;

*   other causes;
    do i= 1 to nrow(obs);
        x = obs[i,];
        y = exp[i,];

        if &ci. then do;
            call nlpfdd(p, grd, h, "pred_oth", est);
            Pr[i]     = p;  
            grad[i,] = grd; 
        end;
        else Pr[i] = pred_oth(est);
    end;

    if &ci. then do;
        pred[,3] = L*Pr;                *   CIF for other causes (if confidence intervals);
        pred[,14] = sqrt(vecdiag(L*grad*cov*t(grad)*t(L)));
        pred[,12] = pred[,3] + 1.96*pred[,14];
        pred[,13] = pred[,3] - 1.96*pred[,14];
    end;
    else    pred[,3] = L*Pr;            *   CIF for other causes (no confidence intervals);

    pred[,4] = pred[,4] + pred[,3];

*   overall survival;
    do i = 1 to nrow(obs);
        x = obs[i,];
        y = exp[i,];
        Pr[i] = pred_s_all(est);
    end;
    pred[,5] = Pr;

*   hazard (lambda);
    do i = 1 to nrow(obs);
        x = obs[i,];
        y = exp[i,];
        Pr[i] = pred_lambda(est);
    end;

    pred[,6] = Pr;              *       hazard (lambda);
    pred[,7] = exp[,2];         *       rate;
    pred[,8] = exp[,1];         *       St_star;

*	if this is conditional analysis, then some of the values will be 0, when they should be missing;
	if &tcond. > 0 then do;
/*		times = pred[,1];*/
		idx = loc(pred[,1] < &tcond.);
		pred[idx,2:4] = .;
		if &ci. then pred[idx,9:14] = .;
	end;
			

*   create sas dataset for further manipulation;
    varnames =  {"&timevar." "&var._d" "&var._o" "&var._all" "&var._s_all" "&var._lambda" "&var._rate" "&var._St_star"};
    if &ci. then varnames = varnames|| {"&var._d_lci" "&var._d_uci" "&var._d_se" "&var._o_lci" "&var._o_uci" "&var._o_se" };

	_pred_ = tablecreate (varnames, pred, {'_study_id_'},ID);
	call TableWriteToDataSet(_pred_, "work", "_pred_"); 

quit;


data _events_;
	merge _events_ _pred_;
	by _study_id_;
run;

%if &debug. eq 0 %then %do;
	proc datasets library = work force nolist ;
    	delete _pred_ _events_t_;
	quit;
	run;
%end;

*   end of processing (also error return);
%cm_fin:

proc optload data = work.myopts; run;
%mend stpm2CM;



/*  predict macro

    compute survival functions post fit, at specified time points, with specified
    covariate values

    parameter
    name        usage
    -----------------
    var         variable to create in output dataset
    measure     hazard  (requires ^at^ parameter if this is not a null model)
                survival (ditto)
                failure (= 1-S(t))
				cure
				uncured
                cumhazard
                xbnobaseline linear predictor (without baseline splines)
                    also known as prognostic index
                xb  linear predictor (includes baseline splines)
                dxb derivative of linear predictor
                deviance    request deviance residuals
                martingale  request martingale residuals
                lifelost
    meansurvwt  variable holding optional weights to be applied to each subject
                for meansurv measure (only valid if meansurv measure is in effect)      
    at          covariate pattern for ^hazard^ and ^survival^ measures
    hdiff1      covariate pattern for first hazard in hazard difference
    hdiff2      default = zero (if hdiff1 is supplied)
    per         scale factor (default = 1) to multiply hazard and hazard difference estimates
    hrnum       Hazard ratio will be predicted
    hrdenom     default zero (if hrnum supplied)
    sdiff1      difference in survival curves will be computed
    sdiff2      defualt zero (if sdiff1 is supplied)
    timevar     if supplied, a time variable to be used instead of the standard _t_
                this variable will be log transformed.
	tcond		(default = 0) non-negative integer. 
    options 	ci / noci (default)
                nodebug (default)   
                debug (1 2 3) different levels of verbosity 

parameters specific to the estimation of lifelost:
    diagage     age at diagnosis (required for lifelost measure)
                default = agediag   
    diagyear    year of diagnosis (required for lifelost measure)
                default = yeardiag
	maxage		default = 99 life table survival estimates for follow-up at ages 
				greater than this value will be replaced with estimates for this age
	maxyear		(default = 2050) life table survival estimates for follow-up at 
				years greater than this value will be replaced with estimates for this year
	by			if grpd option in effect, mean survival expectations will be grouped by 
				unique combinations of the variables passed, rather than based 
				on individual covariate values
    mergeby     list of variables on life table (default is _age, sex, _year)
                that are used to assign probabilities of death
	nodes		(default = 40) number of nodes to use in gaussian quadrature calculations
	tinf		(default = 50) expected and observed (modelled) survival will be calculated
				up tp a minimum of diagage+tinf and maxage
    using       life table data file
	mergeby		specification of variable names on the using file to match with the 
				values specified in diagage, sex and diagyear in the case file. 
				must appear in the order age/sex/year.  Other covariates to be used in
				merging with life table information are appended to these first 3.

    options specific to lifelost measure:
	grpd / nogrpd (the default) use grouped mean expected survival in lifelost estimation
                requires the by variable specification if grpd is chosen

    at (also hdiff1, hdiff2, hrdenom, hrnum, sdiff1, sdiff2  if they are present) will have the following form:
    at = var1:val1 var2:val2  ... (more variable:value pairs) option
        variable var1 is to be held at value val1, etc.
        all other variables not in the above list are held at the ^option^ value
        if option is a number, or retain their values held in the dataset, if option
        equals the missing flag ^.^  If option zero is set , then all variables
        in the parameter list not explicitly set to a value will be held at the value 0
        the constant (if present in the model), is never changed

    if at present, then one of failure, survival, hazard, cumhazard must be present (and vice versa)
    if hdiff1, hrnum or sdiff1 present, then none of at, failure, survival cumhazard or hazard can be present
    per is only valid if hdiff1 or hazard is present
    
    at is not allowed with lifelost option

    ifP = optional where clause to be applied to prediction dataset before estimation
    21 June 2017  fixed bug in meansurv when ifP option in use

	1 Feb 2020
	added cure and uncured paramters
	cure is a measure
	uncured is a qualifier of a measure (which can be survival or hazard)

	4 JUne 2020
	add conditional paramter for survival measures

	25 Feb 2023
	add error traping in proc IML 
*/


%macro predict(var, measure, at=, hdiff1=, hdiff2=zero, hrdenom = zero, hrnum =,
        sdiff1 = , sdiff2 = zero, per = 1, ifP = 1, timevar =, stub = , 
		uncured = 0, cent = 0,
        diagage = agediag, diagyear = yeardiag, survprob = prob,
        maxage = 99, maxyear = 2050,  by =, using =, mergeby=,
        nodes = 40, TINF = 50, TCOND = 0, meansurvwt =,
        options = nogrpd);

%global pred_p iml_vars;        *   we may need the prediction function in a calling program;

    proc optsave out = work.myopts;run;
    options nonotes nomprint nomlogic ;

*   evaluate options chosen;
*   default options;
    %let ci = 0;
    %let debug = 0;
    %let verbose = 0;
	%let err = 0;
    %if "&options." ne "" %then %do;
        %let options = %lowcase(&options.);
        %if %sysfunc(findw(&options., ci)) ^= 0 %then %let ci = 1;      
        %if %sysfunc(findw(&options., debug1)) ^= 0 %then %let debug = 1;
        %else %if %sysfunc(findw(&options., debug2)) ^= 0 %then %let debug = 2;
        %else %if %sysfunc(findw(&options., debug3)) ^= 0 %then %let debug = 3;
    %if &debug. gt 0 %then %do;
        %if %sysfunc(findw(1 2 3,&debug.)) > 0 %then %let verbose = 1;
        %if %sysfunc(findw(2 3,&debug.)) > 0 %then %do;
            options notes;
        %end;
        %if &debug. eq 3 %then %do;
            options spool mprint mlogic;
        %end;
    %end;
    %end;

%let measure = %lowcase(&measure.);
%if "&hdiff1." ne "" %then %let measure = hdiff;
%if "&sdiff1." ne "" %then %let measure = sdiff;
%if "&hrnum." ne "" %then %let measure = hratio;
%let cure = 0;
%if &measure. eq cure %then %let cure = 1;

%if %sysfunc(findw(hazard survival cure dxb failure cumhazard meansurv, &measure.)) = 0 and "&at." ne "" %then 
    %err_mess(the AT parameter is only valid if MEASURE is one of hazard, survival, cumhazard, meansurv or failure);

%if &tcond. ne 0 and %sysfunc(findw(survival meansurv lifelost, &measure.)) = 0 %then
    %err_mess(the tcond parameter is only valid if MEASURE is lifelost, meansurv or survival);

%if %sysfunc(findw(hdiff lifelost sdiff hratio, &measure.)) gt 0 and "&at." ne "" %then 
    %err_mess(AT parameter is not valid with hdiff, sdiff or hratio measures);

%if %sysfunc(findw(hdiff sdiff hratio, &measure.)) gt 0 and "&at." ne "" %then 
    %err_mess(AT parameter is not valid with hdiff, sdiff or hratio measures);
 
%if &meansurvwt. ne  and &measure. ne meansurv %then 
    %err_mess(meansurvwt parameter is only valid if MEASURE is meansurv);

*	if this is lifelost, but using file not found, all bets are off;
%if &measure. eq lifelost and not %sysfunc(exist(&using.)) %then 
	%err_mess(life table file &using. not found); 

%if &err. %then %goto pred_fin ;

*   get model descripton from last fit and set local macros;

%let tvc = ;
%let iml_vars = ;
%let set = _events_;
%let ltime = _ln_t;
%let xb = ;
%let dxb = ;
%let wt = ;
%let wt_str =;
%let covpat =;
%let n_cov = 0;
%let iml_str1 =;
%let iml_str2 =;
%let covar = ;
%let grpd = ;
%let cure_model = 0;
%let tvc = ;
%let n_tvc = 0;


data _null_;
    set _model_;
    if comp ^= '' then call symput(comp, desc);
run;

*	if we are using cure models, then splines are computed in reverse order;
	%let reverse = 0;
	%if &cure_model. %then %let reverse = 1;

*	only use uncured if a cure model has been fit;
	%if &uncured. and not &cure_model. %then %do;
 		%err_mess(uncured option is only available for cure models);
    	%goto pred_fin ;
	%end;	

*	restrictions on use of centile option;
	%if &cent. gt 0 %then %do;
		%if %lowcase(&measure.) ne survival %then %do;
			%err_mess(centile option is only available when requesting survival prediction);
    		%goto pred_fin ;
		%end;	
		%if &cent. ge 100 %then %do;
			%err_mess(centile option must be greater than 0 and less than 100);
    		%goto pred_fin ;
		%end;	
/*		%if &ci. %then %do;*/
/*			%err_mess(confidence limits are not currently available for centile estimation);*/
/*    		%goto pred_fin ;*/
/*		%end;	*/
	%end;

*   if lifelost option, then call specific routine, exit when done;
    %if %substr(%upcase(&measure.),1,2) = LI %then %do;
        %if %sysfunc(findw(&options., grpd)) ^= 0  %then %let grpd = 1;
        %stlifelost (mergeby=&mergeby.,  diagage=&diagage., ci = &ci.,
            diagyear=&diagyear., maxage=&maxage., maxyear=&maxyear.,
            survprob=&survprob., using=&using., grpd=&grpd., by = &by.,
            stub=&stub. , tinf=&tinf., nodes=&nodes., tcond=&tcond.);
        %goto pred_fin;
    %end;

*   residuals (both use cumulative hazard measure);
    %let t_measure =;
    %if %sysfunc(findw(martingale deviance, &measure.)) > 0 %then %do;
        %let t_measure = &measure.;
        %let measure = cumhazard;
    %end;

%if (&cure. or &uncured. or %sysfunc(findw(hazard survival failure cumhazard, &measure.)) gt 0 )
    and "&at." eq ""
    and "&covar." ne "" %then %let at = .;
    
*   create dataset limited to observations with alternate time variable if a variable has been specified;
%let drop_str = ;
%if %str(&timevar.) ne %then %do;
	%let drop_str = drop = &timevar.;
    %let set = _events_t_;
    %let ltime = l_&timevar.;
    data _events_t_;
        set _events_ (where = (&timevar. ^= . and &timevar. > 0));
        l_&timevar. = log(&timevar.);
    run;
%end;
%else %let timevar = _t_;

*   if this is meansurv, then build dataset for analysis;
%if &measure. eq meansurv %then %do;

*   get count of covariate patterns;
    proc means data = _events_ (where = (&ifp.)) noprint nway ;
        class &covar. &meansurvwt. ;
        output out = msurvpop (where =(_stat_ = 'N')
            keep = &covar. _freq_ _stat_  &meansurvwt. 
            rename = (_freq_ = n));
    run;

*   need to know the number of distinct covariate patterns
    and total number of observations;
    proc sql noprint;
        select count(*) into :covpat
        from msurvpop ;

        select sum(n) into :Nt		/*	total count 	*/
        from msurvpop;
    quit;

*   need the covariate list (excluding spline function names);
    %let cv_list =;
    %do iv = 1 %to %sysfunc(countw(&covar.));
        %let cv_list = &cv_list. a.%scan(&covar.,&iv.),;
    %end;

*   build dataset for analysis (each covariate pattern is replicated
    for each unique time point in the analytic dataset);
*   this is where meansurvwt is applied;
    %let t_wt = 1;
    %if &meansurvwt. ne %then %let t_wt = a.&meansurvwt.;
    proc sql noprint;
        create table ms_t as select distinct 
		&cv_list.
		%if not &rcsbaseoff. %then 1 as cons,;
        a.n*&t_wt. as ms_wt,
        log(b.&timevar.) as &ltime.,
        b.&timevar.
        from msurvpop a , 
			(select distinct &timevar. from _events_ where &timevar. not in(0 .)) b

        order by b.&timevar.;
    quit;
    %let set = ms_t;
    %let wt = ms_wt;
    %let wt_str = , covfreq;

*   IML code for meansurv analysis;
    %let iml_str1 = %str(
        use iml_set;
        read all var {&wt.} into w;
        covfreq = t(w[1:&covpat.]);
        pred = j(nrow(obs)/&covpat.,4,.););

    %let iml_str2 = %str(
        do i = 1 to nrow(obs) by &covpat.;
            s0 = i;
            s1 = s0+&covpat.-1;
            p0 = 1+(s0-1)/&covpat.;
            x = obs[s0:s1,] ;
			if obs[s0,1] >= &tcond. then do; 
            	call nlpfdd(p, grd, h, "pred_no_err", est);
            	SE_p =  sqrt(grd*cov*t(grd));
            	pred[p0,2] = exp(p); 
            	pred[p0,3] = exp(p-1.96*se_p); 
            	pred[p0,4] = exp(p+1.96*se_p); 
			end;
			else do;
            	pred[p0,2] = .; 
            	pred[p0,3] = .; 
            	pred[p0,4] = .; 
			end;

            pred[p0,1] = obs[s0,1];         *   time variable;
       end;);
%end;

*   if measure is xbnobaseline, xb, cure, uncured dxb or meansurv, and at is not present,
    then the at subcommand must be constructed with all covariates;
%if (&uncured. or %sysfunc(findw(xbnobaseline xb dxb cure meansurv, &measure.)) > 0)
    and "&at." eq "" %then %do;
    %do i = 1 %to %sysfunc(countw(&covar.));
        %let at = &at. %scan(&covar.,&i.):. ;
    %end;
%end;

*	if we are estimating for a cure, then parameter list is only the covariates
	including constant, if present.  also need the same definition to be
	retained if the the interest is uncured;
%if &cure. or &uncured. %then %do;
	%let th_var = ;
	%if &cure. %then %do;
		%let parms = ;
		%if &int. %then %let parms = cons;
		%let parms = &parms. &covar.; 
		%p_iml(&at., 1, cntl = cure);
		%let pred_p = &xb.;
	%end;

	%if &uncured. %then %do;
		%p_iml(&at., 1, cntl = cure);
		%let pi = exp(-exp(&xb.));
		%p_iml(&at., 0, cntl = uncured);
/*
		uncured survival function from stata
		ln(-(ln(`pi'^(`exprcs') - `pi') - ln(1 - `pi')))

*/
		%if &measure. eq survival %then 
			%let pred_p = log(-(log((&pi.)**exp(&xb.) - &pi.) - log(1 - &pi.)));	
/*
	 	uncured hazard function from stata
	ln(-ln(`pi')*((`drcslist')/`t')*`exprcs'*`pi'^(`exprcs'))- ln(`pi'^(`exprcs') - `pi')
*/

		%if &measure. eq hazard %then 
			%let pred_p = log(-log(&pi.)*((&dxb.)/x[2])*exp(&xb.)*(&pi.)**(exp(&xb.))) 
				- log((&pi.)**exp(&xb.) - (&pi.));
	%end;
%end;			

%if not &cure. and not &uncured. %then %do;
    
*   single AT parameter is present, or no AT parameter is present, and
    neither is any of hrnum, hrdenom, hdiff1, hdiff2, sdiff1, sdiff2;
    %if "&at." ne ""
        or ("&at." eq "" and %sysfunc(findw(hratio hdiff sdiff, &measure.)) = 0)
        %then %do; 
        %p_iml(&at.,1);
    %end; 

*   xb, dxb, xbnobaseline;
    %if %sysfunc(findw(xbnobaseline xb cure, &measure.)) > 0 %then %let pred_p = &xb.;
    %else %if &measure. = dxb %then %let pred_p = &dxb.;

/*     	ln(-ln(`pi')*((`drcslist')/`t')*`exprcs'*`pi'^(`exprcs'))- ln(`pi'^(`exprcs') - `pi') 		*/

*   hazard difference;
    %if &measure. eq hdiff %then %do;
        %p_iml(&hdiff1.,1);
        %let xb0 = &xb.;
        %let dxb0 = &dxb.;
        %p_iml(&hdiff2.,0);
        %let xb1 = &xb.;
        %let dxb1 = &dxb.;
/*      %let measure = hdiff; */
    %end;

*   survival difference;
    %if &measure. eq sdiff %then %do;
        %p_iml(&sdiff1.,1);
        %let xb0 = &xb.;
        %let dxb0 = &dxb.;
        %p_iml(&sdiff2.,0);
        %let xb1 = &xb.;
        %let dxb1 = &dxb.;
/*      %let measure = sdiff; */
    %end;

*   Hazard Ratio;
    %if &measure. eq hratio %then %do;
        %p_iml(&hrnum.,1);
        %let xb0 = &xb.;
        %let dxb0 = &dxb.;
        %p_iml(&hrdenom.,0);
        %let xb1 = &xb.;
        %let dxb1 = &dxb.;
/*      %let measure = hratio; */
    %end;

*   some special strings for theta models;
    %let th_var = ;
    %if &scale. eq theta %then %do;
        %let theta_pos = %eval(1+%sysfunc(countw(&parms.)));
        %let th_var = theta;
    %end;

*   build prediction functions;

    %if %sysfunc(findw(cumhazard survival failure meansurv, &measure.)) gt 0 %then %do;
        %if &scale. ne theta %then %let pred_p = &xb.;
        %else %let pred_p = log(log(e[&theta_pos.]*exp(&xb.)+1)/e[&theta_pos.]);
    %end;

    %if &scale. eq hazard %then %do;
        %if &measure. eq hazard %then
            %let pred_p = -log(x[1]) + log(&dxb.) + (&xb.);
        %else %if &measure. = hdiff %then
            %let pred_p = (1/x[1])*(&dxb0.)*exp(&xb0.) - (1/x[1])*(&dxb1.)*exp(&xb1.);
        %else %if &measure. = hratio %then
            %let pred_p = log(&dxb0.) - log(&dxb1.) + (&xb0.) - (&xb1);
        %else %if &measure. = sdiff %then
            %let pred_p = exp(-exp(&xb0.)) - exp(-exp(&xb1.));
    %end;
    %else %if &scale. eq odds %then %do;
        %if &measure. eq hazard %then
            %let pred_p = -log(x[1]) + log(&dxb.) + (&xb.) -log(1+exp(&xb.));
        %else %if &measure. = hdiff %then
            %let pred_p = (1/x[1]*(&dxb0.)*exp(&xb0.)/((1 + exp(&xb0.))) -
                            1/x[1] *(&dxb1.)*exp(&xb1.)/((1 + exp(&xb1.))));
        %else %if &measure. = hratio %then
            %let pred_p = log(&dxb0.) - log(&dxb1.) + ((&xb0.) - (&xb1.)) - log(1+exp(&xb0.)) + log(1+exp(&xb1.));
        %else %if &measure. = sdiff %then
            %let pred_p = 1/(exp(&xb0.)+1) - 1/(exp(&xb1.)+1);
    %end;
    %else %if &scale. eq normal %then %do;
        %if &measure. eq hazard %then
            %let pred_p = -log(x[1]) + log(&dxb.)
            + log(pdf('normal',(&xb.)))
            - log(cdf('normal',-(&xb.)));
        %else %if &measure. = hdiff %then
            %let pred_p = (1/x[1]*(&dxb0.)*pdf('normal',(&xb0.)/(cdf('normal',-(&xb0.)))
                - 1/x[1] *(&dxb1.)*pdf('normal',(&xb1.))/(cdf('normal',-(&xb1.)));
        %else %if &measure. = hratio %then
            %let pred_p = log(&dxb0.) - log(&dxb1.)
            + (log(pdf('normal',-(&xb0.) - (&xb1.)))
            - log(pdf('normal',(&xb0.)))
            + log(pdf('normal',(&xb1.)));
        %else %if &measure. = sdiff %then
            %let pred_p = cdf('normal',-(&xb0.)) - cdf('normal',-(&xb1.));
    %end;
    %else %if &scale. eq theta %then %do;
        %if &measure. eq hazard %then
            %let pred_p = -log(x[1]) + log(&dxb.) + (&xb.) - log(e[&theta_pos.]*exp(&xb.) + 1);
        %else %if &measure. = hdiff  %then
            %let pred_p = (1/x[1]*((&dxb0.)*exp(&xb0.))/((e[&theta_pos.])*exp(&xb0.) + 1)
                        - 1/x[1]*((&dxb1.)*exp(&xb1.))/((e[&theta_pos.])*exp(&xb1.) + 1));
        %else %if &measure. = hratio %then
            %let pred_p = log(&dxb0.) - log(&dxb1.) + ((&xb0.) - (&xb1.))
                        -log(e[&theta_pos.]*exp(&xb0.) + 1) + log(e[&theta_pos.]*exp(&xb1.) + 1);
        %else %if &measure. = sdiff %then
            %let pred_p = (e[&theta_pos.]*exp(&xb0.) + 1)**(-1/e[&theta_pos.])
                        -(e[&theta_pos.]*exp(&xb1.) + 1)**(-1/e[&theta_pos.]);
    %end;

    %if &measure. = meansurv %then %do;
        %if &scale. = hazard %then %let pred_p = log(covfreq*(exp(-exp(&pred_p.)))/&Nt.);
        %else %if &scale. = odds %then %let pred_p = log(covfreq*(1/(1+exp(&pred_p.)))/&Nt.);
        %else %if &scale = theta %then %let pred_p = log(covfreq*((&theta.*exp(&pred_p.)**(1/&theta.))/&Nt.);
    %end;
%end;

*	conditional predictor for survival
	prediction function transformed to cumulative hazard scale);
%if &tcond. gt 0 and (&measure. eq survival or &measure. eq meansurv) %then %do;
	%let xb_n = &xb;					*	numerator of conditional measure;
	%p_iml(&at.,0, cond = 1);			*	gives us the denominator xb;
	
	%if &measure. = survival %then %let pred_p = exp(&xb.) - exp(&xb_n.);

	%if &measure. = meansurv %then %do ;
        %if &scale. = hazard %then %let pred_p = log(covfreq*(exp(exp(&xb)-exp(&xb_n.)))/&Nt.);
        %else %if &scale. = odds %then %let pred_p = log(covfreq*((1+exp(&xb.))/(1+exp(&xb_n.)))/&Nt.);
        %else %if &scale = theta %then %let pred_p = log(covfreq*((&theta.*exp(&xb_n.-&xb.)**(1/&theta.))/&Nt.);
    %end;

*	build dataset to provide condition;
	data _temp_one_;
		set _events_ (obs = 1);
		&timevar. = &tcond.;
		l_&timevar. = log(&tcond.);
		call symput("l_cond", l_&timevar.);
	run;

*   generate spline variables and any tvc variables;
    %if &rcsbaseoff. = 0 %then %do;
        %rcsgen(l_&timevar.,  dgen=drcs, orthog = &orthog., set = _temp_one_, scalar = &l_cond.,
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

    %if %str(&tvc) ne %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(l_&timevar., gen = rcs_&t., dgen = drcs_&t., knots = %superq(lnknots_&t.),  
                scalar = &l_cond.,set = _temp_one_, orthog = &orthog.,  tmatrix = _T_&t._ );
        %end;
    %end;
%end;


%if &verbose. %then %do;
    %put measure to be estimated is &measure.;
    %put Prediction function = &pred_p.;
    %put IML_var list = &iml_vars.;
    %put Parameter list = &parms.;
    %put first derivative of linear predictor = &dxb.;
	%put first meansurv IML code = &iml_str1.;
	%put second meansurv IML code = &iml_str2.;
%end;

*   generate spline variables and any tvc variables;
    %if &rcsbaseoff. = 0 and not &cure. %then %do;
        %rcsgen(&ltime., dgen=drcs, orthog = &orthog., set = &set., reverse = &reverse., 
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

    %if %str(&tvc) ne and not &cure. %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(&ltime., gen = rcs_&t., dgen = drcs_&t., knots = %superq(lnknots_&t.),  
				reverse = &reverse.,
                set = &set., orthog = &orthog.,  tmatrix = _T_&t._, if2 = &ltime ^= . );
        %end;
    %end;

*   proc IML to compute function of interest, using delta method for SE of estimate
    requires the following macro strings:
    IML_vars    names of observed variables in prediction equation (subset of parms)
    pred_p      prediction equation (expressed in terms of matrix elements
    parms       names of variables in model fit;
%if &measure. ne meansurv %then %do;
	data iml_set;
    	set &set. (where = (&ifp.));
	run;
%end;
%else %if &measure. eq meansurv %then %do;
	data iml_set;
    	set &set. ;
	run;
%end;

%if &cent. = 0 %then %do;
proc iml;

*   variance - covariance matrix from last model fit;
    use _cov_;
    read all var {&parms. &th_var.} into cov;

*   parameter estimates from model fit;
    use _parms_;
    read all var {estimate} into est;

*	if this is a cure prediction, then we only need the rows for the
	estimates and rows and columns for the variance-covariance matrix 
	for the covariates in the model (ie, not the baseline or tvc variables);
	if &cure. then do; 
		cov = cov[1:ncol(cov),1:ncol(cov)];
		est = est[1:ncol(cov),1];
	end;

*   data points (and required variables) at which to estimate function of interest;
    use iml_set ;
    if "&measure." ^= "meansurv" then do;
        read all var {&iml_vars.} into obs;
        read all var {_study_id_} into ID;
	end;
    else read all var {&iml_vars. } into obs;   

*	if conditional survival then bring in the conditioning data, otherwise
	just create a dummy dataset with one row.  It will only be referenced
	in the prediction function if this is actually a conditional estimation;
	s = obs[1,1:ncol(obs)];
	if &tcond. ^= 0 then do;
		use _temp_one_;
		read  var {&iml_vars. } into s;
	end;

*   create diagonal matrix from weights if this is meansurv analysis;
    if "&measure." = "meansurv" then do;
        &iml_str1.;
    end;
    if "&measure." ^= "meansurv" then do;
*       maxtrix to hold predicted values, with elements time, predicted, lcl, ucl;
        pred = j(nrow(obs),4,.);
    end;

*   function of parameters to estimate
        x is vector of observed values
            plus theta constant and weight matrix, if required
        e is vector of parameter estimates;
    start pred(e) global (x , s &wt_str.);
		errFlag = 1;
   		on_error = "if errFlag then do; p = .;  resume; end;";
   		call push(on_error);   /* PUSH code that will be executed if an error occurs */
		p = &pred_p.;
		errFlag = 0;
        return (p);
    finish pred;

*	version of prediction function that can be called by finite difference function nlpfdd();
    start pred_no_err(e) global (x , s &wt_str.);
		p = &pred_p.;
        return (p);
    finish pred_no_err;

*   call finite difference solver to get prediction, 1st derivatives of function
    wrt parameter estimates for each point of interest;
*   if this is not meansurv analysis, then each row is an observation;

*	get last knot, in case we are estimating for uncured;
		knots = {&ln_bhknots.};  
		lastk = knots[ncol(knots)];
		lastk = exp(lastk);

    if "&measure." ^= "meansurv" then do;
        do i= 1 to nrow(obs);
			if obs[i,1] > &tcond. then do;		*	estimates prior to tcond are missing;
            	x = obs[i,] ;
            	if &ci. then do;
					if &uncured. & x[1] >= lastk then do;
						pred[i,2] = .;
                		pred[i,3] = .;
                		pred[i,4] = .;
					end;
*	test that prediction is estimable at this point;
					else if pred(est) ^= . then do;
	                	call nlpfdd(p, grd, h, "pred_no_err", est);    *   returns estimate (p) and 1st derivatives (grd);
	                	SE_p =  sqrt(grd*cov*t(grd));   *   delta method to estimate SE;
	                	pred[i,2] = p;                  *   prediction;
	               		pred[i,3] = p-1.96*se_p;        *   Lower 95% conf. limit;
	               		pred[i,4] = p+1.96*se_p;        *   Upper 95% conf. limit;
					end;
					else do;							*	prediction function returned error;
	               		pred[i,2] = .; 
	               		pred[i,3] = .;
	               		pred[i,4] = .;
					end;
                	pred[i,1] = obs[i,1];           *   time variable;
            	end;
            	else do;
  					if &uncured. & x[1] >= lastk then pred[i,2] = .;
					else pred[i,2] = pred(est);          *   prediction;
              		pred[i,1] = obs[i,1];           *   time variable;
            	end;  
			end; 
        end;
    end;

*   if this is meansurv, then build the matrix of observations;
    if "&measure." = "meansurv" then do;
        &iml_str2.;
        create _pred_ from pred[colname = {"&timevar." 'Pred' 'Lower' 'Upper' }];
        append from pred;
    end;

*   create sas datasets for further manipulation
	need to use table for export because ID variable is allowed to be character;   
    else do;
       	_pred_ = tablecreate ({"&timevar." 'Pred' 'Lower' 'Upper' }, pred, {'_study_id_'},ID);
		call TableWriteToDataSet(_pred_, "work", "_pred_"); 
    end;
quit;

%if &ci. %then %do;
	%let d_ci =;
	%let k_ci = &var._lci &var._uci;
%end;
%else %do;
	%let d_ci = &var._lci &var._uci;
	%let k_ci = ;
%end;

*	new for 3.0.  append new variables to the _events_ file  22 May 2020
	meansurv merges on the time variable.  all other predictions merge on _study_id_;
data _pred_;
   set _pred_;

*   survival (or failure) scale (limits swap because of the - sign);
    if &uncured. or "&measure." in("survival" "failure" "cumhazard" "meansurv" "cure") then do;
        if "&scale." = 'hazard' and &tcond. = 0 then do;
            &var. = exp(-exp(pred));
			if &uncured. and pred = . then pred = 0;
            if &ci. then do;
                &var._uci = exp(-exp(Lower));
                &var._lci = exp(-exp(Upper));
            end;
        end;
        if "&scale." = 'hazard' and &tcond. ^= 0 then do;
            &var. = exp(pred);
            if &ci. then do;
                &var._uci = exp(upper);
                &var._lci = exp(lower);
            end;
        end;
        else if "&scale." = 'odds' then do;
            &var. = 1/(1+exp(pred));
            if &ci. then do;
                &var._uci = 1/(1+exp(Lower));
                &var._lci = 1/(1+exp(Upper));
            end;
        end;
        else if "&scale." = 'theta' then do;
            &var. = exp(-exp(pred));
            if &ci. then do;
                &var._uci = exp(-exp(Lower));
                &var._lci = exp(-exp(Upper));
            end;
        end;
        else if "&scale." = 'normal' then do;
            &var. = cdf('normal',-pred);
            if &ci. then do;
                &var._uci = cdf('normal',-Lower);
                &var._lci = cdf('normal',-Upper);
            end;
        end;
        if "&measure." = "failure" then do;
            &var = 1-&var.;
            if &ci. then do;
                t = &var._uci;
                &var._uci = 1-&var._lci;
                &var._lci = 1-t;
                drop t;
            end;
        end;
        else if "&measure." = "cumhazard" then do;
            &var. = -log(&var.);
            if &ci. then do;
                &var._uci = -log(&var._uci);
                &var._lci = -log(&var._lci);
            end;
        end;
        else if "&measure." = "meansurv" then do;
            &var. = pred;
            if &ci. then do;
                &var._uci = Upper;
                &var._lci = Lower;
            end;
        end;

    end;

*   hazard scale;
    if "&measure." = "hazard"  then do;
        &var. = exp(pred)*&per.;
        if &ci. then do;
            &var._uci = exp(upper)*&per.;
            &var._lci = exp(lower)*&per.;
        end;
    end;

*   differences on the hazard scale;
    if "&measure." = "hdiff"    then do;
        &var. = pred*&per.;
        if &ci. then do;
            &var._lci = lower*&per.;
            &var._uci = upper*&per.;
        end;
    end;

*   hazard ratio scale;
    if "&measure." = "hratio"    then do;
        &var. = exp(pred);
        if &ci. then do;
            &var._lci = exp(lower);
            &var._uci = exp(upper);
        end;
    end;

*   survival difference scale;
    if "&measure." = "sdiff"    then do;
        &var. = pred;
        if &ci. then do;
            &var._lci = lower;
            &var._uci = upper;
        end;
    end;

*   linear predictor, derivative of linear predictor or prognostic index;
    if "&measure." in("xbnobaseline", "xb", "dxb") then do;
        &var. = pred;
        if &ci. then do;
            &var._lci = lower;
            &var._uci = upper;
        end;
    end;

    drop  pred lower upper &d_ci.;;
run;

*	meansurv estimates are merged on time variable.  all other measures are merged 
	on _study_id_;
	%if &measure. ne meansurv %then %do;
		data _events_;
			merge _events_ _pred_;
			by _study_id_;

*	martingale and deviance residuals are computed from the death indicator and the
	cumulative hazard;
		 	if "&t_measure." = "martingale" then
            	&var. = _death_ - &var.;
        	else if "&t_measure." = "deviance" then
				&var. = sign(&var.)*sqrt( -2*(&var. + _death_*(log(_death_ -&var.))));
		run;
	%end;

*	merge meansurv estimates on time variable then re-sort by study ID;
	%if &measure. eq meansurv %then %do;			
		proc sort data = _events_ ; by &timevar.;run;
		data  _events_;
			merge _events_ (in=a)
				_pred_ ;
			by &timevar.;

			if a;
		run;

		proc sort data = _events_; by _study_id_;run;
	%end;  
        
%end;

%if &cent. > 0 %then %do;

	%let tvc_mat_names = ;
	%let tvc_mat_list = ;
	%if &tvc. eq %then %centile(&cent.);
	%else %do;			*	build names of T matrices if orthogonal is on;
						*	also need a comma-separated  list of T matrices;
		%if &orthog. %then %do;
			%do _i_ = 1 %to &n_tvc.;
				%let w = %scan(&tvc.,&_i_.);
				%let tvc_mat_names = &tvc_mat_names. _t_&w._;
				%let tvc_mat_list = &tvc_mat_list.,  _t_&w._;
			%end; 
		%end;
		%centile(&cent.);
	%end;
*	append computed centiles to the _events_ dataset;
	data _events_;
		merge _events_ _pred_;
		by _study_id_;
	run;

%end;		

*   clean up;
    %if &debug. lt 1 %then %do;
        proc datasets library = work nolist ;
        delete _pred_ _events_t_ iml_set _temp_one_ ms_t msurvpop  ;
        quit;
        run;
    %end;

%pred_fin:

proc optload data = work.myopts; run;
%mend;          *    predict;




*   compute spline values for specified scalar time point, and store
    them in global macro strings;
%macro ll_set (ag, _t, _deg =, tmat =, kn = );
*   ag      = scalar value to generate spline values for
    _t      = variable stub to hold spline values generated
    _deg    = number of degrees of freedom
    tmat    = orthogonalisation matrix
    kn      = log knot values;
    
    %global &_t.1 &_t.2 &_t.3 &_t.4 &_t.5 &_t.6 &_t.7 &_t.8 &_t.9 ;
    %local _ind;
    
    %rcsgen(time, gen = &_t. ,  knots = &kn., tmatrix = &tmat., set = _single, scalar = &ag.);
    proc sql noprint;
        %do _ind = 1 %to &_deg-1;
            select &_t.&_ind. into :&_t.&_ind. from _single ;
        %end;
    quit;
%mend;


/*
    program for estimating loss in future expectation of life
    modelled on Stlifelost from Therese Andersson
    
    updated to allow for matrix algebra in IML
    20 Oct 2017
    
    
    limitations:
    models with 60 or fewer parameters.  Can easily be extended to more.  however, see the next point

    there is a hard limit of 64k on the length of macro strings in sas.  The macro strings that
    are built in this routine to compute analytic derivatives of the prediction function approach
    this limit for a model with 29 parameters, using 40 nodes for the numerical integration.  Actual
    length is effected more by TVC variables, and especially for conditional survival models.
    
    if the user gets the message in the log about maximum macro length exceeded, 
    try a model with fewer nodes or modify the system options below 
    (given with their current settings as available from proc options)

    MEXECSIZE=65536    Specifies the maximum macro size that can be executed in memory
    MSYMTABMAX=4194304 Specifies the maximum amount of memory available to the macro variable symbol
                       table or tables
    MVARSIZE=65534     Specifies the maximum size for a macro variable that is stored in memory
*/
%macro stlifelost(mergeby=,  diagage=, 
         diagyear=, maxage=, maxyear=, 
        survprob = prob, using=, grpd=, by=,
        stub=, tinf=, nodes=, tcond=, ci=);
        
%local _i _t save_xb _lt _nt _d _w _px _k  nc pred;

*   variables names for mean observed and mean expected survival, if stub option used;

%if "&stub." eq "" %then %let stub = surv;
%let meanexp = &stub.exp;
%let meanobs = &stub.obs;
%let err = 0;
    
*   install checks as per Stlifelost;
%if &scale. ne hazard or &bhazard. eq %then %do;
    %err_mess(you must fit a relative survival model with the hazard scale before using Lifelost option);
    %goto ll_fin;
%end;   

%if "&grpd." eq "" and "&by." ne "" %then %do;
    %err_mess(You cannot specify the by option without the grpd option);
    %goto ll_fin;
%end;   

%if &nodes. < 0 %then %do;
    %err_mess(You need non-negative number of nodes);
    %goto ll_fin;
%end;

%if &tcond. < 0 %then %do;
    %err_mess(Tcond must be non-negative);
    %goto ll_fin;
%end;

%if &tinf. < 0 %then %do;
    %err_mess(Tinf must be non-negative);
    %goto ll_fin;
%end;


*   must be at least 3 variables in the mergeby parameter;
    %let wc = %sysfunc(countw(&mergeby.,' '));
    %if &wc. lt 3 %then %err_mess(at least 3 variables must be specified (age/sex/year) in mergeby parameter ); 

*   are all variables present in the life table?;
    %let dsid=%sysfunc(open(&using.));
    %if %str(&mergeby.) ne %then %do;
        %let wc = %sysfunc(countw(&mergeby.,' '));  
        %do c = 1 %to &wc.;
            %if not %sysfunc(varnum(&dsid., %scan(&mergeby.,&c.))) %then 
                %err_mess(variable %scan(&mergeby.,&c.) is not in &popmort.);
            %else %let merge&c. = %scan(&mergeby.,&c.);
        %end;
    %end;
    
    %if not %sysfunc(varnum(&dsid., &survprob.)) %then %err_mess(variable &survprob. is not in &using. file); 
    %let rc = %sysfunc(close(&dsid.));  
    
*   test for duplicate rows in popmort file;    
    proc sort data = &using. nodupkeys out = _null_ dupout= _d_mort;
        by &mergeby.;
    run;
    
    %let DSID = %sysfunc(open(_d_mort));
    %let dup_count = %sysfunc(attrn(&DSID, NLOBS)); 
    %let rc = %sysfunc(close(&dsid.));
    %if &dup_count > 0 %then 
        %err_mess(mergeby variables do not uniquely determine rows in &using. file);


*   are all extra merge variables present in the case file?;
    %if &wc. > 3 %then %do;
        %let dsid=%sysfunc(open(_events_));
        %do c = 4 %to &wc.;
            %if not %sysfunc(varnum(&dsid., %scan(&mergeby.,&c.))) %then 
                %err_mess(variable %scan(&mergeby.,&c.) is not in _events_);
        %end;
        %let rc = %sysfunc(close(&dsid.));
    %end;
    

*   check case file for specified variables;
*   may be necessary for some uses of Stlifelost;
    %if not %sysfunc(exist(_events_)) %then %do;
            %err_mess(_events_ file is missing); 
            %goto ll_fin;
        %end;
        
*   file exists.  check for key variables;
    %let dsid = %sysfunc(open(_events_));   
        
    %if not %sysfunc(varnum(&dsid.,_study_id_)) 
        %then %err_mess(unique ID variable _study_id_ is not in _events_ file); 

    %if not %sysfunc(varnum(&dsid.,cons)) 
        %then %err_mess(intercept variable cons is not in _events_ file); 

*   test covariate list (if present) each variable must be in the _events_ file, and of numeric type ;
    %if &covar. ne %then %do;
        %do nc = 1 %to %sysfunc(countw(&covar.," "));
            %let varn = %scan(&covar.,&nc.);    
            %let vnum = %sysfunc(varnum(&dsid.,&varn.)) ;
            %if &vnum eq 0 %then %err_mess(covariate &varn. is not in _events_ file); 
            %else %if %sysfunc(vartype(&dsid., &vnum.)) ne N %then %err_mess(covariate &varn. is not numeric); 
        %end;
    %end;
    %let rc = %sysfunc(close(&dsid.));  

*   test for missing values in specified covariates;
    %if &covar. ne %then %do;
        %let miss_var =;
        %let cov_list = %sysfunc(compbl(&covar.));
        %let cov_list = %sysfunc(translate(&cov_list.,","," "));
        proc sql noprint;
            select nmiss(&cov_list.) into : miss_var
            from _events_;
        quit;
        %if &miss_var. ne 0 %then 
            %err_mess(_events_ file has missing values for some of the specified covariates);
    %end;

    %if &err %then %goto ll_fin;

    
*   loop over each year for each observation until tinf and 
    append the interval-specific expected survival at every year.
    calculate cumulative expected survival at each year;

%let m = %eval(1+&tinf.);
%let b = %eval(&tcond.+1);

*   parse mergeby paramter;
%let mby1 = %scan(&mergeby.,1);
%let mby2 = %scan(&mergeby.,2);
%let mby3 = %scan(&mergeby.,3);

%let merge_add =;
%let extra_var = ;
%if &wc. > 3 %then %do;
    %do i = 4 %to &wc.;
        %let merge_add = &merge_add. 
            and a.%scan(&mergeby., &i.) = b.%scan(&mergeby., &i.);
        %let extra_var = &extra_var. %scan(&mergeby., &i.);
    %end;
%end;

*   first, create a row for every potential year of survival;
*   even if we are not interested in the first few (conditional survival);
data _eventsPlus_;
    set _events_ (keep = _study_id_ &diagage. &diagyear. sex &extra_var.);
    do i = 1 to &tinf.;
        attage = floor(min(&diagage.+ i-1, &maxage.));
        attyr = floor(min(&diagyear. + i-1, &maxyear.));
        output;
    end;
    keep _study_id_ i attage attyr sex &extra_var.;
run; 
        
proc sql noprint;
    create table _events_p_ as
    select a._study_id_,
        a.i,
        b.&survprob.
    from  _eventsPlus_ a left join &using. b
        on a.attage = b.&mby1.
        and a.sex = b.&mby2.
        and a.attyr = b.&mby3
        &merge_add.
    order by a._study_id_, a.i;
    
    select count(*) into : nmiss
    from _events_p_ where &survprob. = .;
quit;

*   check for any unlinked events;
%if &nmiss. > 0 %then %do;
    %err_mess(&nmiss. records failed to match with the population file &using.);
    %goto ll_fin;
%end; 

data _events_s_;
    set _events_p_ (drop = i);
    by _study_id_;

    array S(*) s_star_1 - s_star_&m.;
    if first._study_id_ then do;
        do i = 1 to &m.;
            S(i) = 0;
        end;
        S(&b.) = 1;
        pos = &b.+1;
        c = 1;
    end;
    
    retain s_star_1 - s_star_&m. pos c;

*   only interested in the probabilities after the conditional limit;
    if c >= &b. then do;
        S(pos)=S(pos-1)*&survprob.; 
        pos = pos+1;
    end;
    c = c+1;
    
    if last._study_id_ then output;
    keep _study_id_ s_star_1 - s_star_&m.;
run;

*   grouped mean expected survival;
%let s_set = _events_s_;
%if &grpd. = 1 %then %do;

*   go back and get the by variables;
    data _events_s_;
        merge _events_s_
        _events_ (keep = _study_id_ &by.);
    by _study_id_;
    run;
    
    proc means data = _events_s_ noprint nway;
        class &by.;
    output out = mean_exp_surv mean(s_star_1 - s_star_&m.)=;
    run;

*   replace into observations;
    proc sort data = _events_s_;
        by &by.;
    run;
    
    data exp_s;
        merge  _events_s_ (in=a 
            drop = s_star_1 - s_star_&m.) 
            mean_exp_surv ;
        by &by.;
        if a;
    run;
    
    proc sort data = exp_s;
        by _study_id_;
    run;    
    
    %let s_set = exp_s;
%end;
        
*   find the nodes and the weights for  numerical integration
    using guassian quadrature and then calculate the corresponding time points;
    
proc iml;
    n =  &nodes.;

    muzero = 2;

    A= j(n,n,0);
    do j=1 to n-1;
        b =j/sqrt(4*j**2 -1);
        A[j,j+1] = b;
        A[j+1,j] = b;
    end;
        
    call eigen(nodes, evec, A); 
    
*   calculate time points of interest;
    T = J(n,2,.);
    weights = (evec[1,]##2#muzero)`;
    
    do i = 1 to n;
        T[i,2] = i;
        T[i,1] = (&tinf.-&tcond.)*.5*nodes[i]+(&tinf.+&tcond.)*0.5;
    end;
        
*   sort the times, use the second column to re-order the weights accordingly;
    call sort(T,1);
    w = weights[T[,2],1];
            
    t=T`;   
    w = w`; 
*   need the node times and weights in macro strings;
    s = rowcat(char(t[1,]) +" " );       
    call symputx("times", s); 

    r = rowcat(char(w) +" " );       
    call symputx("weights", r); 
    
*   output time dataset;    
    l_times = log(T[1,]`);
    create _times_ from l_times[colname="l_time"];
    append from l_times;
quit;

*   need a dataset with the times expressed as spline variables;
*   add a final row for the conditional survival time, if present;
%if &tcond. ne 0 %then %do;
    data _times_;
        set _times_ end=last;
        output;
        if last then do;
            l_time = log(&tcond.);
            output;
        end;
    run;
%end;

*   first, for the baseline splines;
    %if &rcsbaseoff. = 0  %then %do;
        %rcsgen(l_time, gen = rcs ,orthog = &orthog., set = _times_,
            knots = &ln_bhknots., tmatrix = _T_bh_);
    %end;

*   and now, for each TVC variable;
    %if %str(&tvc.) ne %then %do;
        %let knots = knots;
        %do itvc = 1 %to &n_tvc.;
            %let t = %scan(&tvc., &itvc.);
            %let df_tvc = %eval(%sysfunc(countw(%superq(knots_&t.)," "))-1);
            %rcsgen(l_time, gen = rcs_&t., knots = %superq(lnknots_&t.),  
                set = _times_, orthog = &orthog.,  tmatrix = _T_&t._);
        %end;
    %end;
            
*   Calculate cumulative expected survival at every time point of interest,
    and multiply with the weights;  
    
data Weighted_s;
    set &s_set.;
    
    array s(*)  S_star_1-S_star_&m.;
    array sw(*) S_w_1 - S_w_&nodes.;

    &meanexp. = 0;

    do is = 1 to &nodes.;
        wgt         = input(scan("&weights.",is,' '),best12.);
        time        = input(scan("&times.",is,' '),best12.);
        floort      = floor(time);
        ceilt       = ceil(time);
        dist        = time-floort;
        Sw(is)      = (s(1+floort)-(s(1+floort)-s(1+ceilt))*dist)*wgt;
        &meanexp.   = &meanexp. + Sw(is);
    end;
    &meanexp. = &meanexp.*(&tinf.-&tcond.)*.5;
    
    keep _study_id_ &meanExp. S_w_1 - S_w_&nodes.; 
run;

*   ready to start building the predictor and gradient strings; 

*   number of main effects to include (including constant, if present);
*   and covariables to load into IML (including constant term);
    %if &covar. eq %then %do;
        %let n_var = 1;
        %let l_covar = cons;
    %end;
    %else %do;
        %let n_var = %eval(%sysfunc(countw(&covar.))+&int.);
        %let l_covar = cons &covar.;
    %end;

*   number of baseline splines;
    %let bh_df =;
    %if &rcsbaseoff. = 0 %then %do;
        %let bh_df = %eval(%sysfunc(countw(&ln_bhknots.," "))-1);
    %end;

*   number of parameters;
    %let n_param = %sysfunc(countw(&parms.));
    
*   number of distinct time points, including conditional survival time, if specified;
    %let nc = &nodes.;
    %if &tcond. > 0 %then %let nc = %eval(&nodes.+1);

*   need a spline pointer.  for this, we need the list of spine names in the _times_ dataset;
    proc contents data = _TIMES_
        out = _vars_ (keep = varnum name)
        noprint;
    run;
    
    proc sql noprint;
        select name, varnum
        into :spl_names separated by ' ', :num separated by ' '
        from _vars_
        where name ^= 'l_time'
        order by varnum;
    quit;  

*   build place holders for up to 60 gradient functions;
*   increase this list if more are required;
%do i = 1 %to 60;
    %let g&i. = ;
%end;

*   if this we have conditional survival, then we need a contirubtion from tcond for each predictor and
    gradient string;
%let xc =;

*   loop over the nodes and build the prediction function;
%do _nt = 1 %to &nodes.;
    
*   contribution from the main effects of each covariate;
    %let xb = obs[,1:&n_var.]*e[1:&n_var.]; 
    %if &tcond. ne 0 %then %let xc = &xb.;  
    
*   spline pointers;
    %let s1 = 1;
    %let s2 = 0;
    
*   estimates pointers;
    %let e1 = %eval(&n_var.+1);
    %let e2 = 0;

*   now, the baseline hazard contribution;  
    %if &rcsbaseoff. = 0 %then %do;
        %let s2 = %eval(&s1.+&bh_df.-1);
        %let e2 = %eval(&e1.+&bh_df.-1);
        %let xb = &xb. + repeat(sp[&_nt.,&s1.:&s2.],nrow(obs))*e[&e1.:&e2.];
        %if &tcond. ne 0 %then %let xc = &xc. + repeat(sp[&nc.,&s1.:&s2.],nrow(obs))*e[&e1.:&e2.];
    %end;
    
*   loop over all the TVC variables;
    %if %str(&tvc.) ne %then %do;
        %do itvc = 1 %to &n_tvc.;
            %let w = %scan(&tvc.,&itvc.);
            %let p = %eval(%sysfunc(findw(&covar., &w.,": ",e))+1);
            %let df_tvc = %eval(%sysfunc(countw(%superq(lnknots_&w.), " "))-1); 
            %let s1 = %eval(&s2.+1);
            %let s2 = %eval(&s1.+&df_tvc.-1);
            %let e1 = %eval(&e2.+1);
            %let e2 = %eval(&e1.+&df_tvc.-1);
            %let xb = &xb. + repeat(sp[&_nt.,&s1.:&s2.],nrow(obs))*e[&e1.:&e2.]#obs[,&p.];
            %if &tcond. ne 0 %then %let xc = &xc. + repeat(sp[&nc.,&s1.:&s2.],nrow(obs))*e[&e1.:&e2.]#obs[,&p.];        
        %end;
    %end;

*   multiply the node-specific prediction function by this nodes weight value
    and accumulate the sum prediction functions over all nodes into one (long) function;

*   this is also the best place to compute the first analytic derivatives of the function with
    respect to each parameter;

    %if &tcond. eq 0 %then %do;
        %if &_nt. eq 1 %then %let pred = exp(-exp(&xb.)+log(s[,&_nt.]));
        %else %let pred = &pred. + exp(-exp(&xb.)+log(s[,&_nt.]));
    %end;
    %else %do;
        %if &_nt. eq 1 %then %let pred = exp(exp(&xc.)-exp(&xb.)+log(s[,&_nt.]));
        %else %let pred = &pred. + exp(exp(&xc.)-exp(&xb.)+log(s[,&_nt.]));
    %end;

*   now for the partial derivatives wrt each parameter; 
    %if &ci. %then %do;
    %do p = 1 %to &n_param;
        %if &_nt = 1 %then %let g&p. = ;
        %let parm = %scan(&parms.,&p.);

        /*  main effect of all covariates   */
        %if &p. le &n_var. %then %do;
            %if &tcond. eq 0 %then
                %let g&p. = &&g&p.-obs[,&p.]#exp(-exp(&xb.))#exp(&xb.)#s[,&_nt];
            %else 
                %let g&p. = &&g&p.+exp(exp(&xc.)-exp(&xb.)+log(s[,&_nt.]))#(obs[,&p.]#(exp(&xc.) - exp(&xb.)));
        %end;

        /*  parameter for baseline spline   */
        %if &parm. in rcs1 rcs2 rcs3 rcs4 rcs5 rcs6 rcs7 rcs8 rcs9 %then %do;
            %let s = %sysfunc(findw(&spl_names., &parm.," ", E));           
            %if &tcond. eq 0 %then
                %let g&p. =  &&g&p.-repeat(sp[&_nt.,&s],nrow(obs))#exp(-exp(&xb.))#exp(&xb.)#s[,&_nt];
            %else 
                %let g&p. = &&g&p.+exp(exp(&xc.)-exp(&xb.)+log(s[,&_nt.]))#(repeat(sp[&nc.,&s],nrow(obs))#exp(&xc.)-repeat(sp[&_nt.,&s],nrow(obs))#exp(&xb.));
        %end;
        
        /*  parameter is tvc spline variable            */
        /*  also need the position of base tvc variable in parms list   */
        %if &p. gt &n_var. and not (&parm. in rcs1 rcs2 rcs3 rcs4 rcs5 rcs6 rcs7 rcs8 rcs9) %then %do;
            %let s = %sysfunc(findw(&spl_names., &parm.," ", E));
            %let base = %substr(&parm., %index(&parm.,_)+1, %length(&parm.)-5);
            %let p1 = %sysfunc(findw(&parms., &base.," ", E));  
            %if &tcond. eq 0 %then          
                %let g&p. =  &&g&p.-repeat(sp[&_nt.,&s],nrow(obs))#obs[,&p1.]#exp(-exp(&xb.))#exp(&xb.)#s[,&_nt];
            %else
                %let g&p. =  &&g&p.+exp(exp(&xc.)-exp(&xb.)+log(s[,&_nt.]))#(repeat(sp[&nc.,&s],nrow(obs))#obs[,&p1.]#exp(&xc.)-repeat(sp[&_nt.,&s],nrow(obs))#obs[,&p1.]#exp(&xb.));
        %end;
        /*  there.  not so hard after all, was it?  */
    %end;
    %end;
%end;

*   adjust for the total time frame and conditional expectation of observed survival
    and create final form of gradient vector for insertion into IML code;
%if &ci. %then %do;
    %do p = 1 %to &n_param;
        %let g&p. = g[,&p.]=(&&g&p.)*(&tinf.-&tcond.)*.5 %str(;);
    %end;
%end;

%let pred = (&pred.)*(&tinf.-&tcond.)*.5 %str(;);

*   IML to carry out the predictions, given the prediction function, the covariate estimates
    (and the var/covar matrix, if needed) and each subjects covariate values;
proc iml;

*   variance - covariance matrix from last model fit;
    use _cov_ (drop = parameter row);
    read all into cov;

*   parameter estimates from model fit;
    use _parms_;
    read all var {estimate} into est;

*   data points (and required variables) at which to estimate function of interest;
    use _events_ ;
    read all var {&l_covar.} into obs;
    read all var { _study_id_} into id;

*   get weighted expected surival at each node (for each subject);
    use Weighted_s;
    read all var("s_w_1" : "s_w_&nodes.") into s;

*   get spline values for node times;
    use _times_ (drop = l_time);
    read all  into sp;

*   maxtrix to hold predicted values, with elements predicted, lcl, ucl;
    pred = j(nrow(obs),3,0);
    
*   function of parameters to estimate
        x   vector of observed values for this subject
        sw  vector of weighted expected surival estimates
        e   vector of parameter estimates;
    start pred(e) global (obs, s, sp);
        p = &pred.;
        return (p);
    finish pred;
    
*   gradient routine will handle a maximum of 60 parameter functions    
    unused gradient functions are null macro strings    ;               
    start grad(e) global (obs, s, sp);
        g = j(nrow(obs),nrow(e),0);     
        &g1. &g2. &g3. &g4. &g5. &g6. &g7. &g8. &g9. &g10.
        &g11. &g12. &g13. &g14. &g15. &g16. &g17. &g18. &g19. &g20.
        &g21. &g22. &g23. &g24. &g25. &g26. &g27. &g28. &g29. &g30.
        &g31. &g32. &g33. &g34. &g35. &g36. &g37. &g38. &g39. &g40.
        &g41. &g42. &g43. &g44. &g45. &g46. &g47. &g48. &g49. &g50.
        &g51. &g52. &g53. &g54. &g55. &g56. &g57. &g58. &g59. &g60.
        return (g);
    finish grad;

*   compute estimates;
    pred[,1] = pred(est);
/*    pred[,4] = id;  */

    
*   confidence intervals, if requested;
    %if &ci. %then %do;
        grd         = grad(est);                        *   compute gradient matrix;
        se_p        = sqrt(vecdiag(grd*cov*t(grd)));    *   delta method to estimate SE;    
        pred[,2]    = pred[,1]-1.96*se_p;               *   Lower 95% conf. limit;
        pred[,3]    = pred[,1]+1.96*se_p;               *   Upper 95% conf. limit;
    %end;

*   output dataset; 
	_pred_ = tablecreate ({'Pred' 'Lower' 'Upper'}, pred, {'_study_id_'},ID);
	call TableWriteToDataSet(_pred_, "work", "_pred_"); 
/*    create _pred_ from pred[colname = { 'Pred' 'Lower' 'Upper' '_study_id_'}];*/
/*    append from pred;*/
quit;

*   merge predicted (observed and expected) survival onto standard dataset
    and compute difference (future years of life lost);
data _events_;
    merge _events_
    weighted_s (keep = _study_id_ &meanexp.)
    _pred_;
    by _study_id_;
        
    &meanobs. = pred;
    &var. = &meanexp. - &meanobs.;

    %if &ci. %then %do;
        &meanobs._lci = lower;
        &meanobs._uci = upper;
        &var._lci = &meanexp. - &meanobs._uci;
        &var._uci = &meanexp. - &meanobs._lci;
    %end;
        
    %if &ci. %then drop lower upper;;
        drop pred;
run;

*   clean up and get out;   
proc datasets  noprint;   
    delete  _pred_ weighted_s _eventsPlus_ _events_p_ _events_s_ 
        _times_ _vars_ exp_s _d_mort;  
quit; run;  
    
*   error exit;
%ll_fin:

%mend;

/*
    p_iml

	version 2.1
	4 june 2020

    common code for prediction functions

    parameter       usage
    ---------------------
    at              describes conditions on covariates
    clear           = 0 retain current values of stored string (iml_vars)
                    = 1 initialise value of stored string
	cntl controls whether to append the covariate prediction strings
					= both 	most measures
					= cure	for cure (and the covariate strings for uncured), just the covariates
					= uncured	for uncured, ony the spline terms (including tvc splines)
	cond			= 0 ordinary use (and default)
					= 1 when called to return an xb where the spline values are referenced from the conditioning matrix

    note use of elementwise matrix product (#) for prediction strings needed for mean survival
    analysis.

	changes occasioned by the necessities of estimation in cure models, where a linear predictor
	composed of covariates only is required.  So, in this version, time variable, constant (if present)	
	and covariates appear first in xb, then are followed by baseline splines (interleaved with their derivatives)
	and the TVC spline terms (and derivatives)

	allow for conditional survival with a scalar time > 0  (use cond = 1)

*/
%macro p_iml(at, clear, cntl = both, cond = 0);

*	variable names for conditional processing;
	%let x = x;						*	matrix name for observed values if not conditional call;
	%if &cond. %then %let x = s;	*	only needed for observed (conditional value of) spline terms;

    %if &clear. %then %let iml_vars = &timevar.;

    %let xb = 0;
    %let dxb = 0;

%if &cntl. ne uncured %then %do;
*   add term for constant if IML_vars has been reset.  Do not do this if the measure
    to be computed is xbnobaseline;
    %if &int. = 1 and &measure. ne xbnobaseline %then %do;
        %if &clear. %then %do;
            %let iml_vars = &iml_vars. cons;
        %end;
        %let xb =  e[1]*x[2];
        %if &measure. eq meansurv %then %let xb =  e[1]#x[,2];
    %end;

*	add the covariate terms;
    %do i = 1 %to &n_cov.;                    *   add covariate patterns (if present);

        %let w = %scan(&covar.,&i.);
        %let pos = %sysfunc(findw(&at., &w.,": ",e));

        %if &pos. > 0 %then %do;
            %let val = *%scan(&at.,%eval(&pos.+1),': ');
            %if "&val." eq "*." %then %do;
                %if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                    %let iml_vars = &iml_vars. &w.;
                %end;
                %let px = %sysfunc(findw(&iml_vars., &w.," ",E));
                %if &measure. eq meansurv %then %let val = #x[,&px.];
                %else %let val = *x[&px.];

            %end;

            %let pe = %sysfunc(findw(&parms., &w.," ",E));
            %if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                %let iml_vars = &iml_vars. &w.;
            %end;
            %let px = %sysfunc(findw(&iml_vars., &w.," ",E));

            %let xb = &xb. + e[&pe.]&val.;
        %end;

        %else %if %sysfunc(findw(&at.,zero," ",E)) = 0 %then %do;   
            %if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                %let iml_vars = &iml_vars. &w.;
            %end;
            %let pe = %sysfunc(findw(&parms., &w.," ",E));
            %let px = %sysfunc(findw(&iml_vars., &w.," ",E));
            %if &measure. eq meansurv %then %let xb = &xb.+ e[&pe.]#x[,&px.];
            %else %let xb = &xb. + e[&pe.]*x[&px.];
        %end;

        %else %if %sysfunc(findw(&at.,zero," ",E)) ^= 0 %then %do;   
            %if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                %let iml_vars = &iml_vars. &w.;
            %end;
            %let pe = %sysfunc(findw(&parms., &w.," ",E));
            %let xb = &xb. + e[&pe.]*0;
        %end;
    %end;
%end;

*   add baseline splines if needed;
%if &cntl. ne cure %then %do;
    %if &rcsbaseoff. = 0 and &measure. ne xbnobaseline %then %do;
        %do i = 1 %to %eval(%sysfunc(countw(&ln_bhknots., " "))-1);

            %let w = rcs&i.;
			%if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then
                    %let iml_vars = &iml_vars. &w.;
            %let pe = %sysfunc(findw(&parms, &w.," ", E));
            %let px = %sysfunc(findw(&iml_vars, &w," ", E));

            %if &measure. eq meansurv %then %let xb = &xb. +e[&pe.]#&x.[,&px.];
            %else %let xb = &xb. +e[&pe.]*&x.[&px.];

            %let w = drcs&i.;
			%if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then
                    %let iml_vars = &iml_vars. &w.;
            %let px = %sysfunc(findw(&iml_vars., &w.," ", E));

            %let w = rcs&i.;
            %let pe = %sysfunc(findw(&parms., &w.," ", E));

            %if &measure. eq meansurv %then %let dxb = &dxb.+ e[&pe.]#x[,&px.];
            %else %let dxb =  &dxb. + e[&pe.]*x[&px.];
        %end;
    %end;

*   add tvc patterns (if present);
    %if %str(&tvc.) ne %then %do;
    	%do i = 1 %to &n_cov.;                    
        	%let w = %scan(&covar.,&i.);
			%if %sysfunc(findw(&tvc., &w.,": ",e)) > 0 %then %do;;
        		%let pos = %sysfunc(findw(&at., &w.,": ",e));
				
       			%if &pos. > 0 %then %do;
            		%let val = *%scan(&at.,%eval(&pos.+1),': ');
            		%if "&val." eq "*." %then %do;
                		%if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                    		%let iml_vars = &iml_vars. &w.;
                		%end;
                		%let px = %sysfunc(findw(&iml_vars., &w.," ",E));
                		%if &measure. eq meansurv %then %let val = #x[,&px.];
                		%else %let val = *x[&px.];
            		%end;

            		%let pe = %sysfunc(findw(&parms., &w.," ",E));
            		%if %sysfunc(findw(&iml_vars., &w.," ",E)) = 0 %then %do;
                		%let iml_vars = &iml_vars. &w.;
      		      	%end;
            		%let px = %sysfunc(findw(&iml_vars., &w.," ",E));
       		 	%end;

            	%else %if "&at." eq "." %then %do;
               	 	%if %sysfunc(findw(&iml_vars., &w.," ", E)) = 0 %then %let iml_vars = &iml_vars. &w.;
               		%let px = %sysfunc(findw(&iml_vars., &w.," ", E));
            		%let val = *x[&px.];
				%end;

				%else %if %sysfunc(findw(&at.,zero," ",E)) ^= 0 %then %let val = *0;   

				%let df_tvc = %eval(%sysfunc(countw(%superq(lnknots_&w.), " "))-1);
				%do n = 1 %to &df_tvc.;
					%let t = rcs_&w.&n.;                    *   this var gets created later;
					%let wr = rcs_&w.&n.;                   *   this is the parameter associated;

					%let pe = %sysfunc(findw(&parms., &wr., " ", e));
					%if %sysfunc(findw(&iml_vars., &t.," ", e)) = 0 %then %do;
						%let iml_vars = &iml_vars. &t.;
						%let pt = %sysfunc(findw(&iml_vars., &t.," ", e));

						%let iml_vars = &iml_vars. drcs_&w.&n.;
					%end;

					%if %sysfunc(findw(&tvc., &w.," ", e)) = 0 %then %do;       *   just add the rcs_var*_tvc..*var term;

						%if &measure. eq meansurv  %then %let xb = &xb.+ e[&pe.]#&x.[,&pt.]#x[,&px.];
						%else %let xb = &xb + e[&pe.]*&x.[&pt.]*x[&px.];      *   just recalculate;

						%let pt = %sysfunc(findw(&iml_vars., &t.," ", e));
						%if &measure. eq meansurv  and &cntl. ne uncured %then %let dxb = &dxb. + e[&pe.]#x[,&pt.]#x[,&px.];
						%else %let dxb = &dxb. + e[&pe.]*x[&pt.]*x[&px.];
					%end;
					%else %do;                  *       add the term rcs_<var>n*_tvc<k>n*val;
						%let pt = %sysfunc(findw(&iml_vars., &t.," ", e));
						%if &measure. eq meansurv  %then %let xb = &xb + e[&pe.]*&x.[,&pt.]&val.;
						%else %let xb = &xb + e[&pe.]*&x.[&pt.]&val.;

						%let t = drcs_&w.&n.;       *   derivative has to be added, too;
						%let pt = %sysfunc(findw(&iml_vars., &t.," ", e));
						%if &measure. eq meansurv  %then %let dxb = &dxb. + e[&pe.]*x[,&pt.]&val.;
						%else  %let dxb = &dxb. + e[&pe.]*x[&pt.]&val.;
					%end;
				%end;
            %end;
        %end;
    %end;
%end;
%mend;          *    p_iml;


/*
    rcsgen

    generate restricted cubic splines.  based on Stata command rcsgen from Paul Lambert

    parameter           usage
    -------------------------
    var                  input variable.  generated splines will describe the cumulative
                         distribution of this variable

    allowed options:

    gen = stub          stubname for generated spline variables (default = rcs)
    dgen = stub         stubname for generated derivatives of spline (default = drcs)
    knots               location of knots on the scale of var (includes boundary knots)
    percentiles         knot points on the percentile scale of the supplied var
                        including boundary percentiles
    orthog              0 (default) do not orthogonalize generated spline variables
                        1  compute orthogonal spline variables  
    tmatrix             use supplied matrix for orthogonalization
    if2                 condition to apply when generating knots
    fw                  name of variable containing weights when generating knots
    scalar              a single value to calculate the spline basis for
    set                 data set name for splines (defaut is _events_)

    df                  degrees of freedom (for default knot positions)
    bknots              boundary knots (default to 0 100)
                        Ony used when specifying df
*/

/*  to account for cure models:  (31 July, 2015 - not implemented)

(in stpm2 module)
allow ^cure^ as parameter in stpm2
scale should be hazard if cure is in use
force ^reverse^ option in rcsgen if cure is in use
with cure, only 3 <= df <= 9
    standard knot placement is df-1, plus add a knot at 95th percentile
extra constraint on last baseline parameter (=0) and for last tvc parameter
(in rcsgen)
*/

/*
	rcsgen macro

	version with reverse option, used by cure models
	24 jan 2020

*/


%macro rcsgen( var, gen = rcs , dgen =  , knots = , orthog = 1, tmatrix = ,
    if1 = 1, percentiles =, if2 =  &var. ^= ., fw = , scalar = , 
    set = _events_, df =, bknots =, reverse = 0);
    
%local m;

%if "&df." ^= "" and ("&knots." ^= "" or "&percentiles." ^= "") %then %do;
    %put Error in specifying required knots.  Do not supply df;
    %put with either of knots or percentiles;
    %goto rcs_fin;
%end;   

%if "&tmatrix." ^= "" and ("&df." ^= "" or "&percentiles." ^= "") %then %do;
    %put Error in specifying parameters.  Do not supply tmatrix;
    %put with either of df or percentiles;
    %goto rcs_fin;
%end;   

*   it is clumsy to get IML to do conditional execution.  Instead, set up conditional
    IML code here as code to be included when proc iml is run
    if a tmatrix has been supplied, then it is used for orthogonalisation, if none has
    been supplied, then weuse an implementation of Paul Lamberts Stata code
    to do the orthogonalisation, and the
    resulting t matrix will be created to be re-used later, if required;

    %if "&tmatrix." eq "" %then %let iml_code = %str(
        meanz = mean(z);
    v = z-meanz || j(nrow(z),1,1);
    q = j(nrow(v),ncol(z),.);
    T = j(ncol(v),ncol(v),0);
    T[ncol(v),] = meanz||1;
    
    do i = 1 to ncol(z);
        r = norm(v[,i])/sqrt(nrow(v));
        q[,i] = v[,i]/r;
        T[i,i] = r;
        do j = i + 1 to ncol(z);
            r = q[,i]` * v[,j]/nrow(v);
            v[,j] = v[,j]-r*q[,i];
            T[i,j] = r;
        end;    
    end;

    T1 = T[1:ncol(z), 1:ncol(z)];
    if "&dgen." ^= "" then dgen = zp*inv(T1);

	if &reverse. then do;
		use _t_knots2_;
		read all into k;

		nk = ncol(k);
		kmin = k[1];
		kmax= k[nk];
		interior = nk-2;
		rcs = j(1,nk,1);

        do _j_= 1 to interior;
			_h_ = nk - _j_;
            _phi_ = (k[_h_]-kmin)/(k[nk]-kmin);
            rcs[_j_] = (k[_h_] > kmax)*(k[_h_] - kmax)**3
               	 - (kmax > kmax)*_phi_*(kmax - kmax)**3
                - (kmin > kmax)*(1-_phi_)*(kmin - kmax)**3;
        end;

		rcs[interior+1] = kmax;
		z1 = rcs*inv(T);
		Q = Q-z1[1,1:ncol(z)];

	end;

    create Tmat from T [colname = varnames];
    append from T;
    
    create _ortho_ from Q [colname = varnames] ;
    append from Q;
);

    %else %let iml_code = %str(
                use &tmatrix.;
    read all into T;
    T_inv= inv(T);
    T1 = inv(T[1:ncol(z),1:ncol(z)]);

    z1 = z||J(nrow(z),1,1);

    Q = z1*T_inv;
    if "&dgen." ^= "" then dgen = zp*T1;

	if &reverse. then do;
		k = {&knots.}; 
		nk = ncol(k);
		kmin = k[1];
		kmax= k[nk];
		interior = nk-2;
		rcs = j(1,nk,1);

        do _j_= 1 to interior;
			_h_ = nk - _j_;
            _phi_ = (k[_h_]-kmin)/(k[nk]-kmin);
            rcs[_j_] = (k[_h_] > kmax)*(k[_h_] - kmax)**3
               	 - (kmax > kmax)*_phi_*(kmax - kmax)**3
                - (kmin > kmax)*(1-_phi_)*(kmin - kmax)**3;
        end;

		rcs[interior+1] = kmax;
		z2 = rcs*inv(T);
		Q = Q-z2;
	end;

    q1 =  Q[,1:ncol(z)];
    create _ortho_ from q1 [colname=varnames];
    append from q1;
        );

/*  standard knot points, based on degrees of freedom and boundary knots        */
/*  there are difficulties with the test on knots if there are decimals in the knot points  */
/*  have to be careful with the test here       */
    %let cen_knot =;
    %if &df. ne %then %do;
        %if %str(&bknots.) eq %then %do;
            %let lower = 0;
            %let upper = 100;
        %end;
        %else %do;
            %let lower = %scan(&bknots.,1);
            %let upper = %scan(&bknots.,2);
        %end;

        %let cen_knot = &lower.;
        %let ik = 1;

/*	----------------------------------------------------------------------------------	*/
/*	use standard knot assignments if this is a tvc variable								*/
/*	not an ideal solution.  Should really have this done in the calling program (stpm2)	*/
/*	----------------------------------------------------------------------------------	*/

		%if &reverse. and "&gen." ne "tvc" %then %let lim = %eval(&df.-1);
		%else %let lim = &df.;
        %do %while(&ik. lt &lim.);
            %let cen_knot = &cen_knot. %sysfunc(round(&ik.*100/&lim.));
            %let ik = %eval(&ik. + 1);
        %end;
		%if &reverse. and "&gen." ne "tvc" %then %let cen_knot = &cen_knot. 95;
        %let cen_knot = &cen_knot. &upper.;
    %end;
    
    %else %if "&percentiles." ne "" %then %let cen_knot = &percentiles.;
    
    %if "&cen_knot." ne "" %then %do;
        proc univariate data = &set.
            (where = (&if1.))
            noprint;
            var &var.;
            output out = _knots_
            PCTLPRE=P
            pctlpts = &cen_knot.;
            %if &fw. ^= %str() %then weight &fw.%str(;);
        run;

        proc transpose
            data = _knots_
            out = _t_knots_ ;
        run;

        proc sql noprint;
            select col1 format best12. into :knots
            separated by ' '
            from _t_knots_;
        quit;

    %end;

*   The number of new spline variables will be one less than the number of knots supplied;
    %let m =  %eval(%sysfunc(countw(&knots.," ")));
    %let n_spl = %eval(&m.-1);
    %let save_knots = &knots.;

*	save knots in temporary dataset for use in IML if needed;
	data _t_knots2_;
    	array k(*) _p1 - _p&m. (&knots.);       *   knot points;
		output;
	run;
		

*   scalar option;
    %if "&scalar." eq "" %then %let v = &var;
    %else %let v = &scalar.;

/*	if spline derivatives required		*/
%if &dgen. ne %then %do;
data &set. ;
    set &set.;

*   arrays for computation of spline functions;

    array z(*)  &gen.1 - &gen.&n_spl. ;                 *   spline variables;
    array zp(*) &dgen.1 - &dgen.&n_spl.;                *   first derivatives;
    array k(*) _p1 - _p&m. (&knots.);       *   knot points;

	kmin = _p1;
	kmax = _p&m.;

	if not &reverse. then do;
    	if &var ^= . then do;
        	z(1)    = &v.;
        	zp(1)   = 1;

        do _j_=2 to &n_spl.;
            _phi_ = (kmax-k(_j_))/(kmax-kmin);
            z(_j_) = (&v. > k(_j_))*(&v. - k(_j_))**3
                - (&v. > kmin)*_phi_*(&v. - kmin)**3
                - (&v. > kmax)*(1-_phi_)*(&v. - kmax)**3;

*       first derivatives of above spline functions;
           	 	zp(_j_) =  (&v. > k(_j_))*3*(&v. - k(_j_))**2
                	- 3*(&v. > kmin)*_phi_*(&v. - kmin)**2
                	- 3*(&v. > kmax)*(1-_phi_)*(&v. - kmax)**2;
        	end;
    	end;
	end;

*	derive spline variables in reverse order;
	else if &reverse. then do;
    	if &var. ^= . then do;
			z(&n_spl.) = &v.;
			zp(&n_spl.) = 1;
			interior = &m. - 2;

        	do _j_= 1 to interior;
				_h_ = &m. - _j_;
            	_phi_ = (k(_h_)-kmin)/(kmax-kmin);
            	z(_j_) = (k(_h_) > &v.)*(k(_h_) - &v.)**3
               	 	- (kmax > &v.)*_phi_*(kmax - &v.)**3
                	- (kmin > &v.)*(1-_phi_)*(kmin - &v.)**3;

            	zp(_j_) =  (k(_h_) >  &v.)*(-3)*(k(_h_) - &v.)**2
               		- (-3)*(kmax > &v.)*_phi_*(kmax - &v.)**2
                	- (-3)*(kmin > &v.)*(1-_phi_)*(kmin - &v.)**2;
        	end;
    	end;
		drop  _h_ interior ;
	end;

    drop _phi_ _j_ _p1 - _p&m. kmin kmax;

run;
%end;

/*	if no spline derivatives required		*/
%if &dgen. eq %then %do;
data &set.;
    set &set. ;

*   arrays for computation of spline functions;

    array z(*)  &gen.1 - &gen.&n_spl. ;     *   spline variables;
    array k(*) _p1 - _p&m. (&knots.);       *   knot points;

	kmin = _p1;
	kmax = _p&m.;

	if not &reverse. then do;
    	if &var ^= . then do;
        	z(1)    = &v.;

        	do _j_=2 to &n_spl.;
            	_phi_ = (kmax-k(_j_))/(kmax-kmin);
            	z(_j_) = (&v. > k(_j_))*(&v. - k(_j_))**3
                	- (&v. > kmin)*_phi_*(&v. - kmin)**3
                	- (&v. > kmax)*(1-_phi_)*(&v. - kmax)**3;
			end;
		end;
	end;

*	derive spline variables in reverse order;
	else if &reverse. then do;
    	if &var. ^= . then do;
			z(&n_spl.) = &v.;
			interior = &m. - 2;

        	do _j_= 1 to interior;
				_h_ = &m. - _j_;
            	_phi_ = (k(_h_)-kmin)/(kmax-kmin);
            	z(_j_) = (k(_h_) > &v.)*(k(_h_) - &v.)**3
               	 	- (kmax > &v.)*_phi_*(kmax - &v.)**3
                	- (kmin > &v.)*(1-_phi_)*(kmin - &v.)**3;
        	end;
    	end;
		drop  _h_ interior ;
	end;

    drop _phi_ _j_ _p1 - _p&m. kmin kmax;
run;
%end;

*   orthogonalisation required;
%if &orthog.  or &tmatrix ne %then %do;
	

proc iml;
    use &set. where (&if2.);

    read all var("&gen.1" : "&gen.&n_spl." ) into z;
    if "&dgen." ^= "" then read all var("&dgen.1" : "&dgen.&n_spl.") into zp;

    varnames = ("&gen.1" : "&gen.&n_spl.");

	&iml_code.;		

    if "&dgen." ^= "" then do;
        varnames2 = ("&dgen.1" : "&dgen.&n_spl.");
        create _orth_d_ from dgen  [colname = varnames2] ;
        append from dgen;
    end;
quit;

%if &dgen. ne %then %let m_l = _orth_d_;
%if &dgen. eq %then %let m_l =;

data &set.;
    merge &set.
        _ortho_
        &m_l.;
run;

proc datasets library = work nolist force nowarn;
    delete _ortho_ _ortho_id_ &m_l. _knots_ _t_knots_ _t_knots2_;
quit;
run;

%end;

*   error conditions branch to here;
%rcs_fin:
%mend rcsgen;

*   macro stset;
*   used to create a standard dataset for time-to-event analysis.  So we dont have to do it
    within the stpm2() macro

    set     = dataset to be used
    death   = variable (and value representing death)
    time    = variable holding time to event data
    id      = study id
    options = noprint .  turn off printing of report.  Default is to print report
	enter	= <variable>  for delayed-entry (period) analysis
				where <variable> holds entry time in same scale as event time
				
    

    this macro will create the dataset _events_, log transform the time variable and
    add a variable for the constant term
    
    Study ID must be unique for each row
    rows with time variable <= 0 or missing will be excluded at this point

	as of 22 May 2020, sort the _events_ file by _study_id_ to facilitate merging
	prediction results

	29 June 2020  added entry time variable
    ;
%macro stset(set, death, time, id, enter =  , options =);

    %let d_var = %scan(&death.,1,'()');
    %let d_val = %scan(&death.,2,'()');
	%let enter_fail = 0;
    
    %let printopt = 1;
    %if %sysfunc(findw(&options., noprint)) ^= 0 %then %let printopt =0; 

	%let m_enter =;
	%if &enter. ne %then 
		%let m_enter = max(&enter.) = max_enter
			sum(enter_fail) = enter_fail
			min(pos_entry) = pos_entry;


data _events_;
    set &set.;


    _death_ = .;
    if &d_var. ^= . then do;
        if &d_var. in( &d_val.) then _death_ = 1;
        else if &d_var. not in( &d_val.) then _death_ = 0;
    end;
    _censor_ = 1-_death_;

    cons = 1;

    _t_ = &time.;
    if _t_ > 0 and _t_ ^= . then _ln_t = log(&time.);
    else do;
        exclude = 1;
        _t_ = .;
		_death_ 	= .;
		_censor_ 	=.;
    end;
    _study_id_ = &id.;

	enter_fail = 0;
	%if &enter. ne  %then %str(
			if &enter. ^= . then do;
				_t0_ = &enter.;
				if &enter. > 0 then pos_entry = &enter.;
			end;
			enter_fail = _t_ < &enter.;

			if enter_fail then do;
				_death_ 	= .;
				_censor_ 	=.;
    		end;	
			);

	remaining = _death_ ^= .; 
run;

%let id =;
%let dup = ;
proc sql noprint;
    select _study_id_, count(*) as  count
        into : id, : dup
        from _events_ 
        group by _study_id_
        having count > 1;
quit;

*   check for duplicated study ID;
%if &dup. ne %then %do;
    %err_mess(duplicate study IDs found:  each row must have a unique identifier);
    %put *** the following IDs were duplicated   ***;
    %put &=id.;
    %go to fin;
%end;

proc means data = _events_ noprint;
    output out = stset_out
        max(_t_) = max_time
        min(_t_) = min_time
		&m_enter.        
		sum(exclude) = excluded
		sum(remaining) = remaining
        sum(_death_) = deaths
        sum(_censor_) = censored;
run;

%if &printopt. %then %do;
%let file = %scan(&set.,1,'(');
title "Summary of events read from &file. dataset, saved in _events_ dataset";
proc print data = stset_out noobs label split = '*';
    var _freq_ 
		excluded 
		%if &enter. ne  %then enter_fail; 
		remaining
		deaths censored max_time min_time 
		%if &enter. ne  %then max_enter pos_entry; ; 

    label excluded  = 'survival*<= 0'
        deaths      = 'Deaths'
        censored    = 'Censored'
		remaining	= 'written* to file'
        max_time    = 'Maximum*survival*time'
        min_time    = 'Minimum*survival*time'
		%if &enter. ne %then %str(max_enter = 'maximum*entry*time'
			enter_fail = 'exit before*entry time'
			pos_entry = 'earliest*entry');
        _freq_      = 'cases*read';
    format max_time min_time %if &enter. ne  %then max_enter pos_entry; 6.3
		_freq_ remaining deaths censored
		%if &enter. ne  %then  enter_fail; comma9.0;
run;

%end;

*   final sort and drop excluded records;
proc sort data = _events_ (where = (exclude = . and enter_fail = 0))
	out = _events_ (drop = exclude _censor_ remaining %if &enter. ne  %then  pos_entry enter_fail; ) ;
  	by _study_id_;
run;

proc datasets library = work nolist nowarn force;
    delete stset_out;
quit;
run;

%fin:

%mend stset;

*   prepend an underscore to a list of variables names;
%macro add_h(set, vars);
    %local ren i;
    %let ren = ;
    %do i = 1 %to %sysfunc(countw(&vars.));
        %let ren = &ren. %scan(&vars.,&i.) = _%scan(&vars.,&i.);
    %end;
    %if &verbose. %then %put &ren.;
        proc datasets library = work nolist nowarn force;
            modify &set.;
            rename &ren.;
        quit;
        run;
%mend add_h;

*   remove an underscore from list of variables names;
%macro rem_h(set, vars);
    %local ren i;
    %let ren = ;
    %do i = 1 %to %sysfunc(countw(&vars.));
        %let ren = &ren. _%scan(&vars.,&i.) = %scan(&vars.,&i.);
    %end;
    %if &verbose. %then %put &ren.;
        proc datasets library = work  nolist nowarn force;
            modify &set.;
            rename &ren.;
        quit;
        run;
%mend rem_h;


/*  MacroStore

    save local macro variables as strings in model description file
    define the dataset, with string lengths if it does not exist, otherwise just append

    parameter   contents
    --------------------
    comp        name of local macro variable
    comm        description of contents of macro variable
*/

%macro ms(comp, comm);
    %if not %sysfunc(exist(_model_)) %then %let def = length comp $ 32 desc $ 1024 comments $64;
    %else %let def = modify _model_;;
    data _model_;
        &def.;
        if _n_ = 1 then do;
            comp = "&comp.";
            desc = symget("&comp.");
            comments = "&comm.";
            output;
        end;
    run;
%mend;


*   use this macro to add an alternate time variable to the standard dataset
    this is a suggestion for the use of the range program in stata
    paramters (all must be supplied)
    var     variable name to add
    min     minimum value
    max     maximum value
    n       number of intervals;
%macro range(var, min, max, n);
data _events_;
    set _events_;
    if _n_ <= &n. then &var. = &min. + (_n_-1)*((&max. - &min.)/(&n.-1));
run;
%mend;

*   macro to turn off all output ;

%macro ODSOff();
ods graphics off;
ods exclude all;
ods noresults;
%mend;


*   macro to turn output back on;
%macro ODSOn();
ods graphics on;
ods exclude none;
ods results;
%mend;

/*  macros to compute natural logs and to exponentiate
    the numbers in a macro strin
*/
/*  exponentiation  */
%macro exp(k);
%local k i ret;
%let ret = ;
%do i = 1 %to %sysfunc(countw(&k., ' '));
%let ret = &ret. %sysfunc(exp(%scan(&k.,&i.,' ')));
%end;
&ret
%mend;

/*  natural log */
%macro ln(k);
%local k i ret;
%let ret = ;
%do i = 1 %to %sysfunc(countw(&k., ' '));
    %let ret = &ret. %sysfunc(log(%scan(&k.,&i.,' ')));
%end;
&ret
%mend;


*   report error message and set an error flag that could be used;
%macro err_mess(message);
    %put ;
    %put ERROR: *** Error in parameter specification ***;
    %put ERROR: *** &message.;
    %put ERROR: ***;
    %let err = 1;
%mend err_mess;

/*	compute centile of a survival curve	*/
/*	added Feb, 2020						*/
%macro centile(alpha);

*	it all happens in IML;  
proc iml;
*	housekeeping before the real work begins;

*	covariance matrix;
	use _cov_;
	read all into temp;
	cov = temp[,2:ncol(temp)];

*	if this is a cure model, then there are missing entries in the covariance matrix to set to zero;
	if &cure_model. then do;
		idx = loc( cov = . );
		cov[ idx ] = 0;
	end;

*	parameter estimates;
	use _parms_;
    read all var {estimate} into est;

/*	e = e`;				*	transpose parameter estimates;*/

*	covariate data and study ID for each subject of interest; 
	if "&timevar." ^= "" then use _events_ where (&timevar. ^=.);
	else use _events_;

	read all var {_t_ cons &covar.} into obs;

	read all var {_study_id_} into ids;
	
*	T matrices for baselines splines and any tvc variables;
	if &orthog. then do;
		if ^&rcsbaseoff. then do;
			use _t_bh_;
			read all  into t_bh;
		end;
		if "&tvc." ^= "" then do;
			tvc_mats = {&tvc_mat_names.};
			do i = 1 to &n_tvc.;
				use (tvc_mats[i]);
				read all into temp;
				call valset(tvc_mats[i], temp);
			end;
		end;
	end;
				
*	build matrix of tvc knots;
	if "&tvc." ^= "" then do;
		use _model_ (where = (substr(comp,1,7) = 'lnknots');
		read all var {comp desc} into tvc_knots;
		k_tvc = j(&n_tvc.,11,0);		*	col 1 is number of knots, 2 - 11 are knots (log scale);
		do i = 1 to &n_tvc.;
			n = countw(tvc_knots[i,2], ' ');
			k_tvc[i,1] = n;
			do j=1 to n;
				k_tvc[i,j+1] = num(scan(tvc_knots[i,2], j, ' '));
			end;
		end;
	end;

*	if this is a cure model, we need the last knot for any spline variable;
	if &cure_model. & &orthog. then do;
		k = {&ln_bhknots.};
		nk = ncol(k);
		kmin = k[1];
		kmax = k[nk];
		rcs = j(1,nk,1);

*	last spline value;
		do j= 1 to nk-2;
			_h_ = nk - j;
         	_phi_ = (k[_h_]-kmin)/(kmax-kmin);
          	rcs[j] = (k[_h_] > kmax)*(k[_h_] - kmax)**3
             - (kmax > kmax)*_phi_*(kmax - kmax)**3
             - (kmin > kmax)*(1-_phi_)*(kmin - kmax)**3;
    	end;
		rcs[nk-1] = kmax;
		rcs = rcs*inv(t_bh);

*	sim. for any tvc variables;
		if "&tvc." ^= "" then do;
		rcs_tvc = j(&n_tvc., 10, 0);
			do i = 1 to &n_tvc.;
				n = k_tvc[i,1];
				k = j(1,n,1);
				k = k_tvc[i,2:n+1];
				kmin = k[1];
				kmax = k[n];
				rcs_t = j(1,n,1);
				do j= 1 to n-2;
					_h_ = n - j;
         			_phi_ = (k[_h_]-kmin)/(kmax-kmin);
          			rcs_t[j] = (k[_h_] > kmax)*(k[_h_] - kmax)**3
             			- (kmax > kmax)*_phi_*(kmax - kmax)**3
             			- (kmin > kmax)*(1-_phi_)*(kmin - kmax)**3;
    			end;
				rcs_t[n-1] = kmax;
				T = value(tvc_mats[i]);
				rcs_t = rcs_t*inv(T);
				rcs_tvc[i,1:n] = rcs_t[1,1:n];
			end;
		end;
	end;

*	prediction function;
	start p_centile(time) global(t_bh, rcs, k_tvc,  rcs_tvc, tvc_mats &tvc_mat_list.,  est, xobs) ;

    	ln_t = log(time);
		x  = xobs;
		e = est;

*   generate spline variables;
    if ^&rcsbaseoff. then do;
		k 	= {&ln_bhknots.};			*	knots;
		z 	= j(1,ncol(k)-1,0);			*	splines;
    	zp 	= j(1,ncol(k)-1,0);     	*	1st derivatives;

		kmin = k[1];
		kmax = k[ncol(k)];

		if ^&reverse. then do;
        	z[1]    = ln_t;
        	zp[1]   = 1;

        	do j=2 to ncol(k)-1;
            	_phi_ = (kmax-k[j])/(kmax-kmin);
            	z[j] = (ln_t > k[j])*(ln_t - k[j])**3
                	- (ln_t > kmin)*_phi_*(ln_t - kmin)**3
                	- (ln_t > kmax)*(1-_phi_)*(ln_t - kmax)**3;

*       first derivatives of above spline functions;
           	 	zp[j] =  (ln_t > k[j])*3*(ln_t - k[j])**2
                	- 3*(ln_t > kmin)*_phi_*(ln_t - kmin)**2
                	- 3*(ln_t > kmax)*(1-_phi_)*(ln_t - kmax)**2;
        	end;
			if &orthog. then do;
				zp = zp*inv(t_bh[1:ncol(z),1:ncol(z)]);
    			z = z||J(1,1,1);
    			z = z*inv(t_bh);
				z = z[1,1:ncol(z)-1];
			end;
			y = j(1,2*ncol(z),0);
			do j = 1 to ncol(z);
				y[1,1+2*(j-1)] = z[1,j];
				y[1,2*j] = zp[1,j];
			end;
			x = x||y;
		end;

		if &reverse. then do;
        	z[ncol(k)-1]    = ln_t;
        	zp[ncol(k)-1]   = 1;
			interior = ncol(k) - 2;

        	do j= 1 to interior;
				h = ncol(k) - j;
            	_phi_ = (k[h]-kmin)/(kmax-kmin);
            	z[j] = (k[h] > ln_t)*(k[h] - ln_t)**3
               	 	- (kmax > ln_t)*_phi_*(kmax - ln_t)**3
                	- (kmin > ln_t)*(1-_phi_)*(kmin - ln_t)**3;

            	zp[j] =  (k[h] >  ln_t)*(-3)*(k[h] - ln_t)**2
               		- (-3)*(kmax > ln_t)*_phi_*(kmax - ln_t)**2
                	- (-3)*(kmin > ln_t)*(1-_phi_)*(kmin - ln_t)**2;
        	end;
			if &orthog. then do;
	    		zp = zp*inv(t_bh[1:ncol(z),1:ncol(z)]);
				z = z||J(1,1,1);
				z = z*inv(t_bh);
				z = z-rcs;
				z = z[1,1:ncol(z)-1];
			end;
			y = j(1,2*ncol(z),0);
			do j = 1 to ncol(z);
				y[1,1+2*(j-1)]	= z[1,j];
				y[1,2*j] 		= zp[1,j];
			end;
			x = x||y;
		end;
	end;

*	generate TVC variables;
	if "&tvc." ^= "" then do;
		do i = 1 to &n_tvc.;
			n = k_tvc[i,1];
/*			k = j(1,n,0);*/
			k = k_tvc[i,2:n+1];
			z 	= j(1,n-1,0);		*	splines;
    		zp 	= j(1,n-1,0);     	*	1st derivatives;

			kmin = k[1];
			kmax = k[n];

			if ^&reverse. then do;
        		z[1]    = ln_t;
        		zp[1]   = 1;

        		do j= 2 to n-1;
            		_phi_ = (kmax-k[j])/(kmax-kmin);
            		z[j] = (ln_t > k[j])*(ln_t - k[j])**3
                		- (ln_t > kmin)*_phi_*(ln_t - kmin)**3
                		- (ln_t > kmax)*(1-_phi_)*(ln_t - kmax)**3;

           	 		zp[j] =  (ln_t > k[j])*3*(ln_t - k[j])**2
                		- 3*(ln_t > kmin)*_phi_*(ln_t - kmin)**2
                		- 3*(ln_t > kmax)*(1-_phi_)*(ln_t - kmax)**2;
				end;
				if &orthog.  then do;
					T = value(tvc_mats[i]);
	    			zp = zp*inv(T[1:ncol(z),1:ncol(z)]);
					z = z||j(1,1,1);
					z = z*inv(T);
					z = z[1,1:ncol(z)-1];
				end;
			end;

			if &reverse. then do;
        		z[n-1]    = ln_t;
        		zp[n-1]   = 1;
				interior = n - 2;

        		do j = 1 to interior;
					h = n - j;
            		_phi_ = (k[h]-kmin)/(kmax-kmin);
            		z[j] = (k[h] > ln_t)*(k[h] - ln_t)**3
               	 		- (kmax > ln_t)*_phi_*(kmax - ln_t)**3
                		- (kmin > ln_t)*(1-_phi_)*(kmin - ln_t)**3;

            		zp[j] =  (k[h] >  ln_t)*(-3)*(k[h] - ln_t)**2
               			- (-3)*(kmax > ln_t)*_phi_*(kmax - ln_t)**2
                		- (-3)*(kmin > ln_t)*(1-_phi_)*(kmin - ln_t)**2;
        		end;
				if &orthog.  then do;				
					T = value(tvc_mats[i]);
	    			zp = zp*inv(T[1:ncol(z),1:ncol(z)]);
					z = z||j(1,1,1);
					z = z*inv(T);
					z = z-rcs_tvc[i,1:n];
					z = z[1,1:ncol(z)-1];
				end;
			end;
			y = j(1,2*ncol(z),0);
				do j = 1 to ncol(z);
				y[1,1+2*(j-1)] = z[1,j];
				y[1,2*j] = zp[1,j];
			end;
			x = x||y;
		end;

	end;

*	account for the scale in use (hazard, odds, theta, normal);
	if "&scale." 		= 'hazard' 	then p = exp(-exp((&pred_p.))) 		- (&alpha/100);
	else if "&scale." 	= 'odds' 	then p = 1/(1+exp(&pred_p.)) 		- (&alpha/100);
	else if "&scale." 	= 'theta' 	then p = exp(-exp((&pred_p.))) 		- (&alpha/100);
	else if "&scale."	= 'normal'	then p = cdf('normal',-(&pred_p.))	- (&alpha/100);

	return (p);
	finish p_centile;
	
	
*	express the centile as a function of the paramter estimates;
	start ci_centile(e) global(est);
		est = e;
		p = froot("p_centile", {&bounds.});
		return (p);
	finish ci_centile;

*	call estimation functions;
	if &ci. then cent = j(nrow(obs),4,0);
	else cent = j(nrow(obs),2,0);
	
	xobs = obs[1,];

    do i= 1 to nrow(obs);
		xobs = obs[i,] ;
		if ^&ci. then cent[i,1] = froot("p_centile", {&bounds.});		*	centile estimate;
		if &ci. then do;
 			call nlpfdd(p, grd, h, "ci_centile", est);
			cent[i,1] = p;					*	centile estimate;
            SE_p =  sqrt(grd*cov*t(grd));   *   delta method to estimate SE;
         	cent[i,3] = p-1.96*se_p;        *   Lower 95% conf. limit;
            cent[i,4] = p+1.96*se_p;        *   Upper 95% conf. limit;
		end;
		cent[i,2] = ids[i,1];   *   study id;
	end;

    if &ci. then create _pred_ from cent[colname = {"&var."  '_study_id_' "&var._lci" "&var._uci"}];
    else create _pred_ from cent[colname = {"&var."  '_study_id_'}];
    append from cent;
quit;

%mend ;

    
/*	sas macro sample 25082:  Determine if a particular variable is present within a data set	*/
/* First parameter is the variable name you want to check. */
/* Second parameter is the name of the dataset to search.  */
%macro varcheck(varname,dsname);

   %let dsid = %sysfunc(open(&dsname));
   %let val = %sysfunc(varnum(&dsid,&varname));
   %let rc = %sysfunc(close(&dsid));

   &val

%mend varcheck;
