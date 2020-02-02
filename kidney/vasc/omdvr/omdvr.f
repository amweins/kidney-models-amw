	PROGRAM OMDVR
C CAPILLARY MODEL FROM ATHERTON ET.AL. AJP 247:F61-72, 1984.
C   FORMULATED AS AN INITIAL VALUE PROBLEM.
C THERE ARE 6 IONIC SPECIES: SODIUM,POTASSIUM,CHLORIDE, AND 
C BICARBONATE, AS WELL AS THE TWO PHOSPHATE IONS, AND AN
C IMPERMEANT SPECIES. IN ADDITION, TWO NONELECTROLYTES,
C GLUCOSE AND UREA, ARE INCLUDED.
C WITHIN THE CAPILLARY THERE ARE 14 HEMOGLOBIN SPECIES.
C  SUBSCRIPTS: C-CAPILLARY  S-INTERSTITIUM
C THE JACOBIAN IS COMPUTED IN DOUBLE PRECISION AND
C  THE MATRIX IS INVERTED IN DOUBLE PRECISION BY LES.
C
C INCLUSION OF NEW SOLUTES FOR COMPATIBILITY WITH NEPHRON MODEL (4/15)
C SPATIAL RESOLUTION OF CS VARIABLES.
C
C  SOLUTE INDEX:
C	1-          NA+
C	2-          K+
C	3-          CL-
C	4-          HCO3-
C	5-          H2CO3
C	6-          CO2
C	7<-5        HPO4--
C	8<-6        H2PO4-
C	9<-8        UREA
C	10-         NH3
C	11-         NH4+
C	12<-9       H+
C	13-         HCO2-
C	14-         H2CO2
C	15<-7       GLUC
C	16<-10      PROT
C
C	17<-11      ONH2 
C       18<-12      NH2  
C	19<-13      ONH3+
C	20<-14      NH3+ 
C	21<-15      ONHCO2-
C	22<-16      NHCO2-
C	23<-17      OXH 
C	24<-18      XH   
C	25<-19      OX-  
C	26<-20      X-   
C	27<-21      OYH  
C	28<-22      YH   
C	29<-23      OY-  
C	30<-24      Y-   
C
	INTEGER SOLS,X,TAU,T,TLIM,COUNT,EXP,CHOP
C
C       SOLS-   NUMBER OF SOLUTES
C	X-	INDICATES INITIAL (=1) OR FORWARD (=2) SPATIAL STEP
C       TAU-    INDICATES STEADY STATE (=0) OR TIMED (=1) EXPT.
C       T-      INDEX FOR OLD (=1) OR NEW (=2) TIME STEP
C       COUNT-  NUMBER OF EXPT. BEING SIMULATED
C       EXP-    TOTAL NUMBER OF EXPTS. IN THE RUN
C	CHOP-	SPATIAL CHOP OF CAPILLARY
C	SWR-	RESETTING SWITCH FOR NEXT SPATIAL (=1) OR TIME (=2) STEP
C
C PARAMETERS
	DOUBLE PRECISION
     1   Z(15),ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15),
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
C 
C VARIABLES
	DOUBLE PRECISION
     1   VC(901,2),PC(901,2),CC(31,901,2),IMPC(901,2),
     1   LCHC(901),XC(15,901),
     1   VS(901,2),PS(901,2),CS(16,901,2),
     1   IMPS(901,2),LCHS(901),XS(15,901),
     1   CCS(15),SAT(901,2),CCN(901,2),CCX(901,2),CCY(901,2),
     1   SC(901,2),AC(901,2),FVC(901,2),FKC(31,901,2),
     1   QV,QC(31),FVCS(901,2),FKCS(16,901,2),
     1   FVCSE(901,2),FKCSE(16,901,2),FVCSI(901,2),FKCSI(16,901,2),
     1   HCTC(901,2),FBC(901,2),POX(901,2)
C
	DOUBLE PRECISION RJ,TPS,TNS,TFS,QPS,QNS,QFS
C
C       Z-      VALENCE OF I'TH SOLUTE
C       ZIMPC-   VALENCE OF SERUM PROTEINS
C       RT-     GAS CONST. TIMES TEMP  (MMHG.ML/MMOL)
C       RTE-    GAS CONST. TIMES TEMP  (JOULE/MMOL)
C       F-      FARADAY (COUL/MOL)
C
C       EPSI-   TOLERANCE FOR THE ERROR VECTOR
C       SCALE-  NUMERICAL SCALING FACTOR
C       DT-     TIME STEP
C       RTAU-   REAL REPRESENTATION OF THE INTEGER TAU
C	TIME-	ELAPSED TIME OF THE SIMULATED EXPT.
C
C	CL-	CAPILLARY LENGTH
C	DX-	SPATIAL STEP OF INTEGRATION
C	RC0-	CAPILLARY RADIUS AT ZERO TRANSCAPILLARY PRESSURE
C	MUC-	CAPILLARY COMPLIANCE
C	ETA-	CAPILLARY FLUID VISCOSITY
C	SC-	LUMINAL CIRCUMFERENCE
C	AC-	LUMINAL AREA
C
C	PK.-	.=C,P,N,F  FOR PK OF BUFFERS, HCO3, HPO4, NH3, OR HCO2
C	KX.-	.=O,D EQUILIBRIUM CONSTANT FOR HGB X-BUFFER
C	KY.-	.=O,D EQUILIBRIUM CONSTANT FOR HGB Y-BUFFER
C	KC.-	.=O,D EQUILIBRIUM CONSTANT FOR NHCO2-HGB BUFFER
C	KD.-	.=O,D EQUILIBRIUM CONSTANT FOR NH2-HGB BUFFER
C	SAT-	HGB OXYGEN SATURATION
C	CCN-	SUM OF HGB-AMINO COMPOUNDS
C	CCX-	SUM OF HGB X-BUFFER
C	CCY-	SUM OF HGB Y-BUFFER
C	A.-     .=0,1,2  LANDIS-PAPPENHEIMER COEFFICIENTS
C
C       LPCS(E,I)-   WATER PERMEABILITY (CM/SEC.MMHG) FOR TJ OR CELLULAR PATH
C       SCS(E,I)-    REFLECTION COEFFICIENT FOR TJ OR CELLULAR PATH
C       HCS(E,I)-    SOLUTE PERMEABILITY (CM/SEC) FOR TJ OR CELLULAR PATH
C       CCS-         MEAN SOLUTE CONCENTRATION
C       LCS(E,I)-    ONSAGER COEFFICIENTS FOR TJ OR CELLULAR PATH
C
C       V.-     VOLTAGE (MVOLT)
C       P.-     PRESSURE (MMHG)
C       C.-     CONCENTRATION (MMOL/ML)
C       X.-     ELECTROCHEMICAL POTENTIAL
C	LCH.-	LOG CONCENTRATION HYDROGEN (PH)
C               . = C,S
C       IMP.-   IMPERMEANT SPECIES CONCENTRATION
C		. = C,S
C
C	FVC-	LUMINAL VOLUME FLOW RATE
C	FKC-	LUMINAL SOLUTE FLUX RATE
C	QC-	NET GENERATION OF SPECIES WITHIN THE LUMEN
C	QV-	NET GENERATION OF VOLUME WITHIN THE LUMEN
C	FVCS(E,I)-  VOLUME REABSORPTIVE RATE (ML/S.CM**2) FOR TJ OR CELLULAR PATH
C	FKCS(E,I)-  SOLUTE REABSORPTIVE RATE (MMOL/S.CM**2) FOR TJ OR CELLULAR PATH
C	HCTC-	CAPILLARY HEMATOCRIT
C	FBC-	CAPILLARY BLOOD FLOW
C
	COMMON/PARVR/ SOLS,T,X,CHOP,
     1   Z,ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS,LCS,HCS,
     1   LPCSE,SCSE,LCSE,HCSE,
     1   LPCSI,SCSI,LCSI,HCSI,
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
	COMMON/VARVR/
     1   VC,PC,CC,VS,PS,CS,IMPC,IMPS,
     1   LCHC,LCHS,XC,XS,CCS,
     1   SAT,CCN,CCX,CCY,HCTC,FBC,POX,
     1   SC,AC,FVC,FKC,QC,QV,FVCS,FKCS,
     1   FVCSE,FKCSE,FVCSI,FKCSI
C
C
	OPEN (UNIT=22,FILE='omdbound.dat',STATUS='old')
	OPEN (UNIT=21,FILE='omdresult.dat')
	OPEN (UNIT=20,FILE='omdparam.dat',STATUS='old')
	OPEN (UNIT=23,FILE='omdend.dat')
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
C
	SOLS=15
	PRINT 27
   27   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
	READ *, EXP
	DO 600 COUNT=1,EXP
C
        CALL OMDNEWT(COUNT,1,1)
C
  600   CONTINUE
C
	STOP
	END
