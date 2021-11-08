/* C-SPAN Formats 
Cancer site definitions follow the 2009 Canadian Cancer Statistics */

*	local format to specify exactly what is considered to be a higher stage.  
	For example, raw text sort order puts stage 'IINOS' to be higher stage than any stage 'III',
	which can be avoided by using an order mapping such as is created here;
proc format;
	invalue cs_order
		'I' =	1
		'IA'=	2
		'IB'=	3
		'II'=	4
		'IIA'=	5
		'IIB'=	6
		'IIC'=	7
		'III'=	8
		'IINOS'=	9
		'IIIA'=	10
		'IIIB'=	11
		'IIIC'=	12
		'IIINOS'=	13
		'IV'=	14
		'IVA'=	15
		'IVB'=	16
		'IVC'=	17
		'IVNOS'=	18
		other = 0;


*	format to group TNM 'best' stage to stage groups I - IV;
	value $ stg
	'I','IA','IB' 						= 'Stage I'
	'II','IIA','IIB','IIC','IINOS' 		= 'Stage II'
	'III','IIIA','IIIB','IIIC','IIINOS' = 'Stage III'
	'IV','IVA','IVB','IVC','IVNOS' 		= 'Stage IV'
	'UNK' 								= 'Stage UNK';

	value $ stg_G
	'I','IA','IB' ,
	'II','IIA','IIB','IIC','IINOS' 		= 'early'
	'III','IIIA','IIIB','IIIC','IIINOS',
	'IV','IVA','IVB','IVC','IVNOS' 		= 'later'
	'UNK' 								= 'Stage UNK';

  /* These age group formats group age into the appropriate categories, depending on cancer site
     so that they can be used with the weights for age-standardization that Larry Ellison from 
     Statistics Canada shared */

  value agegrpn
    0                             = '0-99'
    1                             = '0-39'
    2                             = '40-49'
    3                             = '50-59'
    4                             = '60-69'
    5                             = '70-79'
    6                             = '80-99'
    7                             = '0-34'
    8                             = '0-44'
    9                             = '0-54'
    10                            = '20-44'
    11                            = '35-44'
    12                            = '45-54'
    13                            = '55-64'
    14                            = '65-74'
    15                            = '75-84'
    16                            = '75-99'
    17                            = '85-99'
    18                      	  = '0-44'
    19                     		  = '45-59'
    20                     		  = '60-74'
    21                      	  = '75-99';
;  
run;

proc format;
  ***** Breast (Canadian weights only);
  invalue agegrp_breast
    0 -< 40                      = 1
    40 -< 50                      = 2
    50 -< 60                      = 3
    60 -< 70                      = 4
    70 -< 80                      = 5
    80-high                       = 6;

  ***** Prostate (both Canadian and ICSS weights);
  invalue agegrp_prostate
    0 -<55                       = 9
    55 -<65                       = 13
    65 -<75                       = 14
    75 -<85                       = 15
    85-high                       = 17;

  ***** All other sites;

  invalue agegrp_o
    0 -< 45                      = 8
    45 -< 55                      = 12
    55 -< 65                      = 13
    65 -< 75                      = 14
    75-high                       = 16;
    
  ***** Breast (Canadian weights only);

  value agegrp_breast
    15 -< 40                      = '15-39'
    40 -< 50                      = '40-49'
    50 -< 60                      = '50-59'
    60 -< 70                      = '60-69'
    70 -< 80                      = '70-79'
    80-high                       = '80-99';

  ***** Prostate (both Canadian and ICSS weights);
  value agegrp_prostate
    15 -<55                       = '15-54'
    55 -<65                       = '55-64'
    65 -<75                       = '65-74'
    75 -<85                       = '75-84'
    85-high                       = '85-99';

  ***** All other sites;

  value agegrp_o
    15 -< 45                      = '15-44'
    45 -< 55                      = '45-54'
    55 -< 65                      = '55-64'
    65 -< 75                      = '65-74'
    75-high                       = '75-99';

  
  /* This ICDO-3 format is used with the weights for age-standardization 
	 that Larry Ellison from Statistics Canada shared */

  value $icdoiii
    '0000'='All Cancers'
    '0009'='All nonsex Cancers'
    '0101'='Lip'
    '0102'='Tongue'
    '0103'='Salivary Gland'
    '0104'='Floor of Mouth'
    '0105'='Gum & Oth Mouth'
    '0106'='Nasopharynx'
    '0107'='Oropharynx'
    '0108'='Hypopharynx'
    '0109'='Other Buccal & Phar'
    '0110'='Oral'
    '0199'='Pharynx'
    '0201'='Esophagus'
    '0202'='Stomach'
    '0203'='Small Intestine'
    '0204'='Colon'
    '0205'='Rectum'
    '0206'='Anus'
    '0207'='Liver'
    '0208'='Gallbladder'
    '0209'='Pancreas'
    '0210'='Other Digest Sys'
    '0299'='Colorectal'
    '0301'='Larynx'
    '0302'='Lung and Bronchus'
    '0303'='Other Resp Sys'
    '0400'='Bones & Joints'
    '0500'='Soft tissue'
    '0601'='Skin Melanoma'
    '0602'='Other Skin'
    '0700'='Breast'
    '0801'='Cervix Uteri'
    '0802'='Corpus Uteri'
    '0803'='Uterus, NOS'
    '0804'='Ovary'
    '0805'='Other Female'
    '0901'='Prostate'
    '0902'='Testis'
    '0903'='Penis'
    '0904'='Other Male'
    '1001'='Bladder'
    '1002'='Kidney'
    '1003'='Ureter'
    '1004'='Oth Urinary Sys'
    '1099'='Oth Urin Sys(Ureter)'
    '1100'='Eye'
    '1201'='Brain'
    '1202'='Other Nervous Sys'
    '1301'='Thyroid'
    '1302'='Other Endocrine'
    '1401'='Hodgkin Lymphoma'
    '1402'='N.H.L.'
    '1500'='Multiple Myeloma'
    '1601'='AL Leukemia'
    '1602'='CL Leukemia'
    '1603'='AM Leukemia'
    '1604'='CM Leukemia'
    '1605'='Other Leukemia'
    '1699'='Leukemias'
    '1700'='Unknown'
    '1800'='Mesothelioma'
    '1900'='Kaposi Sarcoma'
    '2000'='Childhood'
    '2199'='All other'
	'5000'='All Cancers';
 

   value CCRgrpf
            .=' '
      100= 'Oral'
      101= 'Lip'
      102= 'Tongue'
      103= 'Major salivary gland'
      104= 'Floor of mouth'
      105= 'Gum and other mouth'
      106= 'Nasopharynx'
      107= 'Oropharynx'
      108= 'Hypopharynx'
      109= 'Other buccal cavity and pharynx'

      201= 'Esophagus' 
      202= 'Stomach'
      203= 'Small intestine'
      200= 'Colorectal'
      204= 'Colon excluding rectum'
      205= 'Rectum and rectosigmoid'
      206= 'Anus'
      207= 'Liver'
      208= 'Gallbladder'
      209= 'Pancreas' 
      210= 'Other digestive system'
            
      301= 'Larynx'
      302= 'Lung and bronchus'
      303= 'Other respiratory system'

      400= 'Bones and joints'
   
      500= 'Soft tissue (including heart)'

      601= 'Melanomas of the skin'
      602= 'Other skin'

      700= 'Breast'

      801=  'Cervix uteri'
      800=  'Body of Uterus'
      802=  'Corpus uteri'
      803=  'Uterus, NOS'
      804=  'Ovary'
      805=  'Other female genital system'

      901= 'Prostate'
      902= 'Testis'
      903= 'Penis'
      904= 'Other male genital system'
         
      1001=  "Bladder (including in situ)" 
      1002=  'Kidney and renal pelvis'
      1003=  'Ureter'
      1004=  'Other urinary system'
               
      1100= 'Eye'

	  1200 = 'Brain and CNS'
      1201=  'Brain'
      1202=  'Other nervous system'

      1301=  'Thyroid'
      1302=  'Other endocrine'

      1401=  "Hodgkin Lymphoma"
      1402=  "Non-Hodgkin lymphomas"
            
      1500 = 'Multiple myeloma'
   
      1600= 'Leukemia'
      1601= 'Acute lymphocytic Leukemia'
      1602= 'Chronic lymphocytic Leukemia'
      1603= 'Acute myeloid leukemia'
      1604= 'Chronic myeloid leukemia'
      1605= 'Other Leukemias'

      1800= 'Mesothelioma'

      1900= 'Kaposi Sarcoma'

      1700= "Other, ill-defined, and unknown sites" 

	  5000 = "All Cancers"

      9999= 'All Other Cancers';

   Value CCRaggf  
      101-109 = '100'
      201 = '201'
      202 = '202'
      204-205 = '200'       
      207 = '207'
      209 = '209'
      301 = '301'
      302 = '302'
      601 = '601'
      700 = '700'
      801 = '801'
      802-803 = '800'
      804 = '804'
      901 = '901'
      902 = '902'
      1001 = '1001'
      1002 = '1002'
      1201 - 1202 = '1201' 
      1301 = '1301'
      1401 = '1401'
      1402 = '1402'
      1500 = '1500'
      1601-1605 = '1600'
      other = '9999'
                  ;
run;

