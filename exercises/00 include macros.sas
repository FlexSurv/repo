/*	00 include macros		

	set these locations to load all the macros needed
	for the exercises

	March 12, 2017	
	
	Use this version for local SAS installs.  Folders use Windows-syle syntax
*/

*	location of directory tree with the following structure:
	root
		- data		location of datasets used in exercises
		- exercises	sas versions of the stata .do exercises from workshop
		- macros	sas macros for relative survival and flexible parametric regression;

* WHERE fpsaus == flexible parametric survival analysis using SAS;
	
%let fpsaus = <path to the root of the survival resources tree>;

/*	for example:	*/
/*	%let fpsaus = C:\work\survival resources;  */

*	macros used in exercises;
%include "&fpsaus.\macros\other.sas";
%include "&fpsaus.\macros\regression methods.sas";
%include "&fpsaus.\macros\life table methods.sas";
%include "&fpsaus.\macros\formats.sas";

libname data "&fpsaus.\data";

