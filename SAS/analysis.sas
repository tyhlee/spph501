/*  Import the prepped data */


FILENAME REFFILE '/folders/myfolders/sasuser.v94/copd_sas_corrected.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV
	OUT=WORK.IMPORT;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.IMPORT; RUN;


%web_open_table(WORK.IMPORT);
PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV REPLACE
	OUT=WORK.copd;
	GETNAMES=YES;
	
RUN;

PROC CONTENTS DATA=WORK.copd;
RUN;
***Analysis;


/* TITLE "RE Frailty, Gap time, Weibull"; */
/* PROC NLMIXED DATA=WORK.copd gCONV=0 QPOINTS=5; */
/* 	PARMS ln_gamma=0.1180 */
/* 	b_0= -0.8 */
/* 	b_trtSal=0.3127 */
/* 	b_trtSal_Flu=0.2507 */
/* 	b_age=-0.5213	 */
/* 	b_BMI=-1.0989	 */
/* 	b_fev1pp=-1.4856 */
/* 	b_fev1fvcratio=0.2384	 */
/* 	b_rand_season2=-0.4905 */
/* 	b_rand_season3=-0.2933 */
/* 	b_rand_season4=-0.1329 */
/* 	b_num_re1=-0.4762	 */
/* 	b_num_re2=-0.3804 */
/* 	b_num_re3=-0.9259	 */
/* 	b_event_season2=-0.1268	 */
/* 	b_event_season3=-0.3455	 */
/* 	b_event_season4=-0.3801 */
/* 	b_male=0.0057 */
/* 	b_nowsmk=0.0614 */
/* 	b_oxygen=0.4871 */
/* 	b_OPTIMAL=0.1067 */
/* 	ln_v=0.7565; */
/*  */
/* 	lin=b_0+ */
/* 		b_trtSal*trtSal+ */
/* 		b_trtSal_Flu*trtSal_Flu+ */
/* 		b_age*age+ */
/* 		b_BMI*BMI+ */
/* 		b_fev1pp*fev1pp+ */
/* 		b_fev1fvcratio*fev1fvcratio+ */
/* 		b_rand_season2*rand_season2+ */
/* 		b_rand_season3*rand_season3+ */
/* 		b_rand_season4*rand_season4+ */
/* 		b_num_re1*num_re1+ */
/* 		b_num_re2*num_re2+ */
/* 		b_num_re3*num_re3+ */
/* 		b_event_season2*event_season2+ */
/* 		b_event_season3*event_season3+ */
/* 		b_event_season4*event_season4+ */
/* 		b_male*male+ */
/* 		b_nowsmk*nowsmk+ */
/* 		b_oxygen*oxygen+ */
/* 		b_OPTIMAL*OPTIMAL+ */
/* 		z; */
/* 	alpha=EXP(lin); */
/*  */
/* 	S_t=exp(-(alpha*gaptime)**EXP(ln_gamma));  */
/* 	h=EXP(ln_gamma)*alpha*((gaptime*alpha)**(EXP(ln_gamma)-1));  */
/*   */
/* 	ll=(recurrent_event~=1)*log(S_t)+(recurrent_event=1)*log(h*S_t);  */
/*   */
/* 	MODEL gaptime~general(ll);  */
/* 	RANDOM z~NORMAL(0,EXP(ln_v)) SUBJECT=ID OUT=__z_re;  */
/*   */
/* 	PREDICT lin OUT=__lin_re;  */
/* 	 */
/* 	ESTIMATE 'v' EXP(ln_v); */
/* 	ESTIMATE 'gamma' EXP(ln_gamma); */
/* 	 */
/* 	ESTIMATE 're_trtSal' b_trtSal*EXP(ln_gamma); */
/* 	ESTIMATE 're_trtSal_Flu' b_trtSal_Flu*EXP(ln_gamma); */
/* 	ESTIMATE 're_age' b_age*EXP(ln_gamma); */
/* 	ESTIMATE 're_BMI' b_BMI*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1pp' b_fev1pp*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1fvcratio' b_fev1fvcratio*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season2' b_rand_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season3' b_rand_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season4' b_rand_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re1' b_num_re1*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re2' b_num_re2*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re3' b_num_re3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season2' b_event_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season3' b_event_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season4' b_event_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_male' b_male*EXP(ln_gamma); */
/* 	ESTIMATE 're_nowsmk' b_nowsmk*EXP(ln_gamma); */
/* 	ESTIMATE 're_oxygen' b_oxygen*EXP(ln_gamma); */
/* 	ESTIMATE 're_OPTIMAL' b_OPTIMAL*EXP(ln_gamma); */
/*  */
/* 	ODS OUTPUT ParameterEstimates=__param_re;  */
/* 	ODS OUTPUT AdditionalEstimates=__add_param_re; */
/* RUN;  */
 
/* TITLE "RE Gamma Frailty, Gap time, Weibull"; */
/* PROC NLMIXED DATA=WORK.copd gCONV=0 QPOINTS=5; */
/* 	PARMS ln_gamma=0.1180 */
/* 	b_0= -0.8 */
/* 	b_trtSal=0.3127 */
/* 	b_trtSal_Flu=0.2507 */
/* 	b_age=-0.5213	 */
/* 	b_BMI=-1.0989	 */
/* 	b_fev1pp=-1.4856 */
/* 	b_fev1fvcratio=0.2384	 */
/* 	b_rand_season2=-0.4905 */
/* 	b_rand_season3=-0.2933 */
/* 	b_rand_season4=-0.1329 */
/* 	b_num_re1=-0.4762	 */
/* 	b_num_re2=-0.3804 */
/* 	b_num_re3=-0.9259	 */
/* 	b_event_season2=-0.1268	 */
/* 	b_event_season3=-0.3455	 */
/* 	b_event_season4=-0.3801 */
/* 	b_male=0.0057 */
/* 	b_nowsmk=0.0614 */
/* 	b_oxygen=0.4871 */
/* 	b_OPTIMAL=0.1067 */
/* 	theta=1; */
/* 		 */
/* 	p=cdf("NORMAL",z); */
/* 	if p > .99999 then p=0.99999; */
/* 	g2=quantile('GAMMA',p,1/theta); */
/* 	g=g2 * theta; */
/*  */
/* 	lin=b_0+ */
/* 		b_trtSal*trtSal+ */
/* 		b_trtSal_Flu*trtSal_Flu+ */
/* 		b_age*age+ */
/* 		b_BMI*BMI+ */
/* 		b_fev1pp*fev1pp+ */
/* 		b_fev1fvcratio*fev1fvcratio+ */
/* 		b_rand_season2*rand_season2+ */
/* 		b_rand_season3*rand_season3+ */
/* 		b_rand_season4*rand_season4+ */
/* 		b_num_re1*num_re1+ */
/* 		b_num_re2*num_re2+ */
/* 		b_num_re3*num_re3+ */
/* 		b_event_season2*event_season2+ */
/* 		b_event_season3*event_season3+ */
/* 		b_event_season4*event_season4+ */
/* 		b_male*male+ */
/* 		b_nowsmk*nowsmk+ */
/* 		b_oxygen*oxygen+ */
/* 		b_OPTIMAL*OPTIMAL+ */
/* 		log(g); */
/* 	alpha=EXP(lin); */
/*  */
/* 	S_t=exp(-(alpha*gaptime)**EXP(ln_gamma));  */
/* 	h=EXP(ln_gamma)*alpha*((gaptime*alpha)**(EXP(ln_gamma)-1));  */
/*   */
/* 	ll=(recurrent_event~=1)*log(S_t)+(recurrent_event=1)*log(h*S_t);  */
/*   */
/* 	MODEL gaptime~general(ll);  */
/* 	RANDOM z~NORMAL(0,1) SUBJECT=ID OUT=__z_re_gamma;  */
/*   */
/* 	PREDICT lin OUT=__lin_re_gamma;  */
/* 	 */
/* 	ESTIMATE 'theta_inv' 1/theta; */
/* 	ESTIMATE 'gamma' EXP(ln_gamma); */
/* 	 */
/* 	ESTIMATE 're_trtSal' b_trtSal*EXP(ln_gamma); */
/* 	ESTIMATE 're_trtSal_Flu' b_trtSal_Flu*EXP(ln_gamma); */
/* 	ESTIMATE 're_age' b_age*EXP(ln_gamma); */
/* 	ESTIMATE 're_BMI' b_BMI*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1pp' b_fev1pp*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1fvcratio' b_fev1fvcratio*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season2' b_rand_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season3' b_rand_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season4' b_rand_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re1' b_num_re1*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re2' b_num_re2*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re3' b_num_re3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season2' b_event_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season3' b_event_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season4' b_event_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_male' b_male*EXP(ln_gamma); */
/* 	ESTIMATE 're_nowsmk' b_nowsmk*EXP(ln_gamma); */
/* 	ESTIMATE 're_oxygen' b_oxygen*EXP(ln_gamma); */
/* 	ESTIMATE 're_OPTIMAL' b_OPTIMAL*EXP(ln_gamma); */
/*  */
/* 	ODS OUTPUT ParameterEstimates=__param_re_gamma;  */
/* 	ODS OUTPUT AdditionalEstimates=__add_param_re_gamma; */
/* RUN;   */
 
/* TITLE "Death Frailty, Gap time, Weibull";  */
/* PROC NLMIXED DATA=WORK.copd gCONV=0 QPOINTS=5 MAXTIME=500 MAXIT=80 MAXITER=80 GTOL=1E-7;  */
/* 	PARMS ln_gamma_d=0.2295 */
/* 	b_d_0=-8	  */
/* 	b_d_trtSal=-0.15 */
/* 	b_d_trtSal_Flu=-0.50 */
/* 	b_d_age=0.08  */
/* 	b_d_packyears=0.66	  */
/* 	b_d_BMI=-2.7 */
/* 	b_d_fev1pp=-3.3  */
/* 	b_d_fev1fvcratio=0.61	 */
/* 	b_d_num_re1=0.55	  */
/* 	b_d_num_re2=1.179  */
/* 	b_d_num_re3=1.413 */
/* 	b_d_male=0.49 */
/* 	b_d_nowsmk=0.61 */
/* 	b_d_oxygen=0.47 */
/* 	ln_v=0; */
/* 		 */
/* 	 */
/* 	lin_d=b_d_0 +	 */
/* 		b_d_trtSal*trtSal+ */
/* 		b_d_trtSal_Flu*trtSal_Flu+ */
/* 		b_d_age*age+ */
/* 		b_d_packyears*packyears+ */
/* 		b_d_BMI*BMI+ */
/* 		b_d_fev1pp*fev1pp+ */
/* 		b_d_fev1fvcratio*fev1fvcratio+ */
/* 		b_d_num_re1*num_re1+ */
/* 		b_d_num_re2*num_re2+ */
/* 		b_d_num_re3*num_re3+ */
/* 		b_d_male*male+ */
/* 		b_d_nowsmk*nowsmk+ */
/* 		b_d_oxygen*oxygen+ */
/* 		z; */
/*  */
/* 	alpha_d=EXP(lin_d); */
/*  */
/* 	S_t_d=exp(-(alpha_d*gaptime)**EXP(ln_gamma_d)); */
/* 	h_d=EXP(ln_gamma_d)*alpha_d*((gaptime*alpha_d)**(EXP(ln_gamma_d)-1)); */
/*  */
/* 	ll=(death~=1)*log(S_t_d)+(death=1)*log(h_d*S_t_d); */
/*  */
/* 	MODEL gaptime~general(ll); */
/* 	RANDOM z~NORMAL(0,exp(ln_v)) SUBJECT=ID OUT=__z_death_normal; 	 */
/*  */
/* 	PREDICT lin_d OUT=__lin_death_normal; */
/* 	 */
/* 	ESTIMATE 'v_d' exp(ln_v); */
/* 	ESTIMATE 'gamma_d' exp(ln_gamma_d); */
/* 	 */
/* 	ESTIMATE 'd_trtSal' b_d_trtSal*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_trtSal_Flu' b_d_trtSal_Flu*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_age' b_d_age*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_packyears' b_d_packyears*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_BMI' b_d_BMI*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_fev1pp' b_d_fev1pp*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_fev1fvcratio' b_d_fev1fvcratio*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re1' b_d_num_re1*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re2' b_d_num_re2*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re3' b_d_num_re3*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_male' b_d_male*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_nowsmk' b_d_nowsmk*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_oxygen' b_d_oxygen*EXP(ln_gamma_d); */
/* 	 */
/* 	ODS OUTPUT ParameterEstimates=__param_death_normal;  */
/* 	ODS OUTPUT ParameterEstimates=__add_param_death_normal; */
/* RUN; */

 
TITLE "Death Gamma Frailty, Gap time, Weibull"; 
PROC NLMIXED DATA=WORK.copd gCONV=0 QPOINTS=10 MAXIT=80 MAXITER=80 GTOL=1E-5 corr cov; 
	PARMS ln_gamma_d=0.1548	
	b_d_0=-2.7168
	b_d_trtSal=-0.11
	b_d_trtSal_Flu=-0.2866	
	b_d_age=0.4341
	b_d_packyears=0.8038
	b_d_BMI=-0.3106	
	b_d_fev1pp=-3.72
	b_d_fev1fvcratio=-0.2173
	b_d_num_re1=1.18
	b_d_num_re2=1.6245
	b_d_num_re3=1.8640
	b_d_male=0.56
	b_d_nowsmk=0.6185
	b_d_oxygen=0.5416
	theta=1;
		
	p=cdf("NORMAL",z);
	if p > .99999 then p=0.99999;
	g2=quantile('GAMMA',p,1/theta);
	g=g2 * theta;
	
	lin_d=b_d_0 +	
		b_d_trtSal*trtSal+
		b_d_trtSal_Flu*trtSal_Flu+
		b_d_age*age+
		b_d_packyears*packyears+
		b_d_BMI*BMI+
		b_d_fev1pp*fev1pp+
		b_d_fev1fvcratio*fev1fvcratio+
		b_d_num_re1*num_re1+
		b_d_num_re2*num_re2+
		b_d_num_re3*num_re3+
		b_d_male*male+
		b_d_nowsmk*nowsmk+
		b_d_oxygen*oxygen+
		log(g);

	alpha_d=EXP(lin_d);

	S_t_d=exp(-(alpha_d*gaptime)**EXP(ln_gamma_d));
	h_d=EXP(ln_gamma_d)*alpha_d*((gaptime*alpha_d)**(EXP(ln_gamma_d)-1));

	ll=(death~=1)*log(S_t_d)+(death=1)*log(h_d*S_t_d);

	MODEL gaptime~general(ll);
	RANDOM z~NORMAL(0,1) SUBJECT=ID OUT=__z_death_gamma; 	

	PREDICT lin_d OUT=__lin_death_gamma;
	
	ESTIMATE 'theta_inv' 1/theta;
	ESTIMATE 'gamma_d' exp(ln_gamma_d);
	
	ESTIMATE 'd_trtSal' b_d_trtSal*EXP(ln_gamma_d);
	ESTIMATE 'd_trtSal_Flu' b_d_trtSal_Flu*EXP(ln_gamma_d);
	ESTIMATE 'd_age' b_d_age*EXP(ln_gamma_d);
	ESTIMATE 'd_packyears' b_d_packyears*EXP(ln_gamma_d);
	ESTIMATE 'd_BMI' b_d_BMI*EXP(ln_gamma_d);
	ESTIMATE 'd_fev1pp' b_d_fev1pp*EXP(ln_gamma_d);
	ESTIMATE 'd_fev1fvcratio' b_d_fev1fvcratio*EXP(ln_gamma_d);
	ESTIMATE 'd_num_re1' b_d_num_re1*EXP(ln_gamma_d);
	ESTIMATE 'd_num_re2' b_d_num_re2*EXP(ln_gamma_d);
	ESTIMATE 'd_num_re3' b_d_num_re3*EXP(ln_gamma_d);
	ESTIMATE 'd_male' b_d_male*EXP(ln_gamma_d);
	ESTIMATE 'd_nowsmk' b_d_nowsmk*EXP(ln_gamma_d);
	ESTIMATE 'd_oxygen' b_d_oxygen*EXP(ln_gamma_d);
	
	ODS OUTPUT ParameterEstimates=__param_death_gamma; 
	ODS OUTPUT ParameterEstimates=__add_param_death_gamma;
RUN;

/* TITLE "Shared Frailty, Gap time, Weibull"; */
/* PROC NLMIXED DATA=WORK.copd GCONV=0 QPOINTS=6 GTOL=1E-8 COV CORR; */
/* 	PARMS ln_gamma=0.1183 ln_gamma_d=0.1417	 */
/* 	b_0=-0.7753	 */
/* 	b_trtSal=0.32 */
/* 	b_trtSal_Flu=0.244 */
/* 	b_age=-0.513 */
/* 	b_BMI=-1.14 */
/* 	b_fev1pp=-1.515 */
/* 	b_fev1fvcratio=0.25 */
/* 	b_rand_season2=-0.49 */
/* 	b_rand_season3=-0.30 */
/* 	b_rand_season4=-0.14 */
/* 	b_num_re1=-0.44 */
/* 	b_num_re2=-0.35 */
/* 	b_num_re3=-0.776 */
/* 	b_event_season2=-0.125 */
/* 	b_event_season3=-0.3564	 */
/* 	b_event_season4=-0.37 */
/* 	b_male=0.066 */
/* 	b_nowsmk=0.07 */
/* 	b_oxygen=0.516 */
/* 	b_OPTIMAL=0.126 */
/* 	b_d_0=-3.2636	  */
/* 	b_d_trtSal=-0.10 */
/* 	b_d_trtSal_Flu=-0.2645 */
/* 	b_d_age=1.87 */
/* 	b_d_packyears=0.5176 */
/* 	b_d_BMI=-1.2 */
/* 	b_d_fev1pp=-4.376 */
/* 	b_d_fev1fvcratio=-0.05880		 */
/* 	b_d_num_re1=1.16  */
/* 	b_d_num_re2=1.53  */
/* 	b_d_num_re3=1.87 */
/* 	b_d_male=0.578 */
/* 	b_d_nowsmk=0.65 */
/* 	b_d_oxygen=0.57 */
/* 	kinv=50 */
/* 	ln_v=0.7128; */
/* 	 */
/* 	lin=b_0+ */
/* 		b_trtSal*trtSal+ */
/* 		b_trtSal_Flu*trtSal_Flu+ */
/* 		b_age*age+ */
/* 		b_BMI*BMI+ */
/* 		b_fev1pp*fev1pp+ */
/* 		b_fev1fvcratio*fev1fvcratio+ */
/* 		b_rand_season2*rand_season2+ */
/* 		b_rand_season3*rand_season3+ */
/* 		b_rand_season4*rand_season4+ */
/* 		b_num_re1*num_re1+ */
/* 		b_num_re2*num_re2+ */
/* 		b_num_re3*num_re3+ */
/* 		b_event_season2*event_season2+ */
/* 		b_event_season3*event_season3+ */
/* 		b_event_season4*event_season4+ */
/* 		b_male*male+ */
/* 		b_nowsmk*nowsmk+ */
/* 		b_oxygen*oxygen+ */
/* 		b_OPTIMAL*OPTIMAL+ */
/* 		z; */
/* 		 */
/* 	alpha=EXP(lin); */
/*  */
/* 	S_t=exp(-(alpha*gaptime)**EXP(ln_gamma)); */
/* 	h=EXP(ln_gamma)*alpha*((gaptime*alpha)**(EXP(ln_gamma)-1)); */
/*  */
/* 	ll1=(recurrent_event~=1)*log(S_t)+(recurrent_event=1)*log(h*S_t); */
/*  */
/* 	lin_d=b_d_0 +	 */
/* 		b_d_trtSal*trtSal+ */
/* 		b_d_trtSal_Flu*trtSal_Flu+ */
/* 		b_d_age*age+ */
/* 		b_d_packyears*packyears+ */
/* 		b_d_BMI*BMI+ */
/* 		b_d_fev1pp*fev1pp+ */
/* 		b_d_fev1fvcratio*fev1fvcratio+ */
/* 		b_d_nowsmk*nowsmk+ */
/* 		b_d_oxygen*oxygen+ */
/* 		b_d_num_re1*num_re1+ */
/* 		b_d_num_re2*num_re2+ */
/* 		b_d_num_re3*num_re3+ */
/* 		b_d_male*male+ */
/* 		z/kinv; */
/*  */
/* 	alpha_d=EXP(lin_d); */
/*  */
/* 	S_t_d=exp(-(alpha_d*gaptime)**EXP(ln_gamma_d)); */
/* 	h_d=EXP(ln_gamma_d)*alpha_d*((gaptime*alpha_d)**(EXP(ln_gamma_d)-1)); */
/*  */
/* 	ll2=(death~=1)*log(S_t_d)+(death=1)*log(h_d*S_t_d); */
/* 	 */
/* 	ll=ll1+ll2; */
/*  */
/* 	MODEL gaptime~general(ll); */
/* 	RANDOM z~NORMAL(0,EXP(ln_v)) SUBJECT=ID OUT=__joint_Z; */
/* 	 */
/* 	PREDICT lin OUT=__joint_lin; */
/* 	PREDICT lin_d OUT=__joint_lin_d; */
/* 	 */
/* 	Estimate 'k' 1/kinv; */
/* 	Estimate 'gamma' EXP(ln_gamma); */
/* 	Estimate 'gamma_d' EXP(ln_gamma_d); */
/* 	Estimate 'v' EXP(ln_v); */
/* 	 */
/* 	ESTIMATE 're_trtSal' b_trtSal*EXP(ln_gamma); */
/* 	ESTIMATE 're_trtSal_Flu' b_trtSal_Flu*EXP(ln_gamma); */
/* 	ESTIMATE 're_age' b_age*EXP(ln_gamma); */
/* 	ESTIMATE 're_BMI' b_BMI*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1pp' b_fev1pp*EXP(ln_gamma); */
/* 	ESTIMATE 're_fev1fvcratio' b_fev1fvcratio*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season2' b_rand_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season3' b_rand_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_rand_season4' b_rand_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re1' b_num_re1*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re2' b_num_re2*EXP(ln_gamma); */
/* 	ESTIMATE 're_num_re3' b_num_re3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season2' b_event_season2*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season3' b_event_season3*EXP(ln_gamma); */
/* 	ESTIMATE 're_event_season4' b_event_season4*EXP(ln_gamma); */
/* 	ESTIMATE 're_male' b_male*EXP(ln_gamma); */
/* 	ESTIMATE 're_nowsmk' b_nowsmk*EXP(ln_gamma); */
/* 	ESTIMATE 're_oxygen' b_oxygen*EXP(ln_gamma); */
/* 	ESTIMATE 're_OPTIMAL' b_OPTIMAL*EXP(ln_gamma); */
/* 	 */
/* 	ESTIMATE 'd_trtSal' b_d_trtSal*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_trtSal_Flu' b_d_trtSal_Flu*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_age' b_d_age*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_packyears' b_d_packyears*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_BMI' b_d_BMI*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_fev1pp' b_d_fev1pp*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_fev1fvcratio' b_d_fev1fvcratio*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re1' b_d_num_re1*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re2' b_d_num_re2*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_num_re3' b_d_num_re3*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_male' b_d_male*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_nowsmk' b_d_nowsmk*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_oxygen' b_d_oxygen*EXP(ln_gamma_d); */
/*  */
/*  */
/* 	ESTIMATE 'CON1_0' b_num_re1*EXP(ln_gamma); */
/* 	ESTIMATE 'CON2_1' b_num_re2*EXP(ln_gamma)-b_num_re1*EXP(ln_gamma); */
/* 	ESTIMATE 'CON3_2' b_num_re3*EXP(ln_gamma)-b_num_re2*EXP(ln_gamma); */
/*  */
/* 	ESTIMATE 'd_CON1_0' b_d_num_re1*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_CON2_1' b_d_num_re2*EXP(ln_gamma_d)-b_d_num_re1*EXP(ln_gamma_d); */
/* 	ESTIMATE 'd_CON3_2' b_d_num_re3*EXP(ln_gamma_d)-b_d_num_re2*EXP(ln_gamma_d); */
/* 	 */
/* 	ODS OUTPUT ParameterEstimates=__joint_estimates; */
/* 	ODS OUTPUT AdditionalEstimates=__joint_additional; */
/* RUN; */