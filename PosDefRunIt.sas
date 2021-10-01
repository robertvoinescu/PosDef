%macro PosDefRunIt(InputTable=SCORR,
			OutputTable=SCORRPD,
			NameCol=_name_,
			eigenLT	    =0.000000001,
			eigenReplace=0.00000001,
			maxIter=1000,
			method=1,
			final_reb=0);
			
	%if %symexist(PythonPosDefCounter)=0 %then %do;
		%global PythonPosDefCounter;
		%let PythonPosDefCounter=0;
	%end;
	%let PythonPosDefCounter=%sysevalf(&PythonPosDefCounter.+1);
	

	%let PosDefSubDir  = Powersimm\sasmacro\ForwardPriceSim;
	%let WorkDirectory = %sysfunc(getoption(work));
	%let LogFile 	   = "&OutputLogPath\JobId"||trim(left("&JobId."))||"_"||compress(put(datetime(),datetime18.),' :')||"_posdef.log" 

    data PosDefParams;
		OutputTable     = "&OutputTable.&PythonPosDefCounter."    ;
		NameCol         = "&NameCol."        ;
		eigenLT         =  &eigenLT.         ;
		eigenReplace    =  &eigenReplace.    ;
		maxIter         =  &maxIter.         ;
		method          =  &method.          ;
		final_reb	=  &final_reb.	     ;
    run;

    
   	proc export data=PosDefParams outfile= "&WorkDirectory.\PosDefParams&PythonPosDefCounter..csv" dbms=csv replace;
	run;

	proc export data=&InputTable. outfile= "&WorkDirectory.\&InputTable.&PythonPosDefCounter..csv" dbms=csv replace;
	run;

	filename pos "&WorkDirectory.\RunPythonPosDef&PythonPosDefCounter..bat";

	data _null_;
		file pos;
		pythonpath = %sysfunc(quote("C:\Program Files\Python37\python.exe"));
		msgline = pythonpath || " &BookMacroCodeBase.\&PosDefSubDir.\PosDefRunIt.py &InputTable.&PythonPosDefCounter. PosDefParams&PythonPosDefCounter. &WorkDirectory. &LogFile.";
		put msgline; 
	run;


	/* RV: call python executable to do find nearest covariance matrix requires numpy and pandas installed*/
	options noxwait xsync;
	x "&WorkDirectory.\RunPythonPosDef&PythonPosDefCounter..bat";

	
	/* RV: Read out the python data */
	proc import datafile="&WorkDirectory.\&OutputTable.&PythonPosDefCounter..csv" out=&OutputTable. dbms=csv replace;
		guessingrows=max;
	run;
%mend PosDefRunIt;
