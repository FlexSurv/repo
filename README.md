# repo
 SAS macros for life table and flexible parametric survival analysis of cancer patient survival
 
 The SAS© macros presented here are based on the Stata© user-contributed programs 
 stpm2 (by Paul Lambert, 2008 and following) and its ancillary programs, in particular the StLifeLost 
 module (by Therese Andersson, 2012). This repository contains the macro documentation 
 and annotated examples of the use of the macros to explore the flexible parametric 
 survival model, supporting datasets and background material. 
 
 Also present in this repository is %rel_surv, a SAS program developed with the support 
 of the Canadian Partnership Against Cancer. This program has been re-written and enhanced
 from the original (due to Paul Dickman) and is designed for life table analysis of cancer patient 
 survival. It includes Pohar Perme estimator of net survival as one of its potential outputs.

 Note that for regression analysis, the SAS module IML must be present in the local SAS 
 installation as it is required for evaluating the prediction functions and generating 
 the restricted cubic splines used in model fitting.

folder structure present:

 exercises:  		using the macros in this repo for the analysis of cancer survival data
 data:  			public-use datasets used in the excercises
 documentation:  	annotations for the exercises and macro documentation
 macros:  			all macros and other sas code needed for the exercises
 other programs:  	the example program from the Appendix to LEL paper (2023)
  
Although the programs and macros presented here have been tested and compared with results from 
other similar programs, they cannot be guaranteed to be free of errors or omissions. 
 30 April 2023
