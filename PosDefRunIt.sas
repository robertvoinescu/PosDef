%macro PosDefRunIt(InputTable=SCORR,
			OutputTable=SCORRPD,
			NameCol=_name_,
			eigenLT=0.001,
			eigenReplace=0.01,
			maxIter=12,
			method = 1);
			
	%if %symexist(PythonPosDefCounter)=0 %then %do;
		%global PythonPosDefCounter;
		%let PythonPosDefCounter=0;
	%end;
	%let PythonPosDefCounter=%sysevalf(&PythonPosDefCounter.+1);

	%if %symexist(overrideLT)=0 %then %do;
		%let eigenLT=overrideLT;
	%end;

	%if %symexist(algo)=0 %then %do;
		%let algo="original";
	%end;

	%let finalReb=0;
	%if %symexist(overrideReb)=0 %then %do;
		%let finalReb=overrideReb;
	%end;
	
	%let PosDefSubDir = Powersimm\sasmacro\ForwardPriceSim;
	%let WorkDirectory = %sysfunc(getoption(work));

    data PosDefParms;
		OutputTable     = "&OutputTable.&PythonPosDefCounter."    ;
		NameCol         = "&NameCol."        ;
		eigenLT         =  &eigenLT.         ;
		eigenReplace    =  &eigenReplace.    ;
		maxIter         =  &maxIter.         ;
		method          =  &method.          ;
		finalReb        =  &FinalReb.        ;
    run;
    
   	proc export data=PosDefParms outfile= "&WorkDirectory.\PosDefParms&PythonPosDefCounter..csv" dbms=csv replace;
	run;

	proc export data=&InputTable. outfile= "&WorkDirectory.\&InputTable.&PythonPosDefCounter..csv" dbms=csv replace;
	run;

	filename pos "&WorkDirectory.\RunPythonPosDef&PythonPosDefCounter..bat";

	%if algo = "original" %then %do;
		data _null_;
			file pos;
			pythonpath = %sysfunc(quote("C:\Program Files\Python37\python.exe"));
			msgline = pythonpath || " &BookMacroCodeBase.\&PosDefSubDir.\PosDefRunIt.py &InputTable.&PythonPosDefCounter. PosDefParms&PythonPosDefCounter. &WorkDirectory.";
			put msgline; 
		run;
	%else %if algo = "matlab";

	/* RV: call python executable to do find nearest covariance matrix requires numpy and pandas installed*/
	options noxwait xsync;
	x "&WorkDirectory.\RunPythonPosDef&PythonPosDefCounter..bat";
	
	/* RV: Read out the python data */
	proc import datafile="&WorkDirectory.\&OutputTable.&PythonPosDefCounter..csv" out=&OutputTable. dbms=csv replace;
	run;
%mend PosDefRunIt;
