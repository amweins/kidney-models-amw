	PROGRAM AHLM
C
C BASED ON PROGRAM COMP
C COMP:COMPARTMENT MODEL FOR EPITHELIAL FLUXES
C REVISED 12/12/94 (STARTING WITH PROXIMAL TUBULE MODEL IN NAH)
C TO REPRESENT THE CORTICAL COLLECTING DUCT EPITHELIUM
C WILL SUPERCEDE THE STRIETER MODEL BY ADDING CO2, H2CO3, NH3, NH4
C
	INTEGER SOLS,TAU,T,TLIM,COUNT,EXP,NPAR,X,CHOP
C
C       SOLS-   NUMBER OF SOLUTES
C       TAU-    INDICATES STEADY STATE (=0) OR TIMED (=1) EXPT.
C       T-      INDEX FOR OLD (=1) OR NEW (=2) TIME STEP
C       COUNT-  NUMBER OF EXPT. BEING SIMULATED
C       EXP-    TOTAL NUMBER OF EXPTS. IN THE RUN
C	NPAR-	NUMBER OF PARAMETERS ICREMENTED IN A TRANSIENT PROBLEM.
C	X-	INDICATES SPATIAL STEP
C	CHOP-	SPATIAL CHOP OF TUBULE
C
	CHARACTER*1 ISOFM
C
C	ISOFM-	ISOFORM OF NKCC AND KCC TRANSPORTERS (=B, A, OR F)
C
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,DIST,
     1   EPSI,DT,RTAU,TIME,L0,L(1601,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5),BCO2(3)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(1601,2),PM(1601,2),CM(15,1601,2),
     1   IMPM(1601,2),LCHM(1601),XM(15,1601),
     1   VS(1601,2),PS(1601,2),CS(15,1601,2),
     1   IMPS(1601,2),LCHS(1601),XS(15,1601),
     1   TL,DX,RM0,MUM,ETA,ZIMPS,
     1   SM(1601,2),AM(1601,2),FVM(1601,2),FKM(16,1601,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(1601,2),MUA,CHVL0,CHVL(1601,2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(1601,2),PE(1601,2),CE(15,1601,2),LCHE(1601),XE(15,1601),
     1   FEVM(1601,2),FEKM(15,1601,2),FEVS(1601,2),FEKS(15,1601,2),
     1   CURE(1601)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,1601,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,1601,2),HCBUF(3,1601,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),
     1   ATMI(3,15,1601),ATIS(3,15,1601),ATIE(3,15,1601),
     1   VI(3,1601,2),PI(3,1601,2),CI(3,15,1601,2),IMP(3,1601),
     1   LCHI(3,1601),XI(3,15,1601),FIVM(3,1601,2),FIKM(3,15,1601,2),
     1   FIVS(3,1601,2),FIKS(3,15,1601,2),CURI(3,1601),JV(3,1601,2),
     1   JK(3,15,1601,2)
C SPECIAL TRANSPORTERS
	CHARACTER*1 ISOFM
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,1601,2),JHK(3,1601,2),JHP(3,1601,2),
     1   JAE1(3,1601,2),JTSC(3,1601,2),JNHE3(3,3,1601,2),QIAMM,
     1   NNKCC(3),NKCL(3),JNKCC(3,4,1601,2),JKCC(3,3,1601,2)
C     INDICES FOR THE CELL TYPE: 1=PRINCIPAL; 2=ALPHA; 3=BETA
C
C       Z-      VALENCE OF I'TH SOLUTE
C       ZIMP-   VALENCE OF IMPERMEANT SPECIES
C       RT-     GAS CONST. TIMES TEMP  (MMHG.ML/MMOL)
C       RTE-    GAS CONST. TIMES TEMP  (JOULE/MMOL)
C       F-      FARADAY (COUL/MOL)
C
C       EPSI-   TOLERANCE FOR THE ERROR VECTOR
C       DT-     TIME STEP
C       RTAU-   REAL REPRESENTATION OF THE INTEGER TAU
C	TIME-	ELAPSED TIME OF THE SIMULATED EXPT.
C
C       L-      CELL HEIGHT
C       L0-     CELL MEMBRANE HEIGHT
C       MUA-    CHANNEL AREA COMPLIANCE
C       MUV-    CHANNEL VOLUME COMPLIANCE
C
C       AI0-    CROSS-SECTION AREA OF CELL AT MUCOSAL SURFACE
C       AME-    AREA OF TIGHT JUNCTION
C       AE-     CHANNEL AREA AT BASE
C       AIE-    LATERAL CELL AREA
C
C       CHVL-   CHANNEL VOLUME
C       CLVL-   CELL VOLUME
C       CHVL0-  REFERENCE CELL VOLUME
C       CLVL0-  REFERENCE CELL VOLUME
C       IMP0-  CELL IMPERMEANT SPECIES CONC. AT REFERENCE VOLUME
C       IMP-    CELL IMPERMEANT SPECIES CONC. AT CLVL
C
C	CBUF-	CELL BUFFER CONCENTRATION
C	HCBUF-	PROTONATED CELL BUFFER
C	TBUF-	TOTAL CELL BUFFER CONCENTRATION FOR THE REFERENCE VOLUME
C	PK.-	.=B,C,P  FOR PK OF CELL BUFFERS, HCO3, OR HPO4
C	KHY.,KDHY.-  .=P,A,B,E HYDRATION/DEHYDRATION CONSTANT FOR CO2
C
C       LP..-   WATER PERMEABILITY (CM/SEC.MMHG)
C       S..-    REFLECTION COEFFICIENT
C       H..-    SOLUTE PERMEABILITY (CM/SEC)
C       C..-    MEAN SOLUTE CONCENTRATION
C	        ME-     TIGHT JUNCTION
C               MI-     CELL MUCOSA
C               IE-     LATERAL CELL
C               ES-     CHANNEL BASEMENT MEMBRANE
C               IS-     CELL SEROSA
C
C       LM,LS-  MATRICES OF ONSAGER COEFFICIENTS FOR MUCOSAL
C               AND BASOLATERAL CELL MEMBRANES
C
C       V.-     VOLTAGE (MVOLT)
C       P.-     PRESSURE (MMHG)
C       C.-     CONCENTRATION (MMOL/ML)
C       X.-     ELECTROCHEMICAL POTENTIAL
C	LCH.-	LOG CONCENTRATION HYDROGEN (PH)
C               . = M,S,I,E
C
C       IMP.-   IMPERMEANT SPECIES CONCENTRATION
C		. = M,S
C
C       NP-     NA+ EXTRUSION OUT THE LATERAL CELL MEMBRANE
C	KNPN-	NA AFFINITY FOR THE PUMP
C	KNPK-	K AFFINITY FOR THE PUMP
C       KNH4-   RELATIVE AFFINITY FOR NH4+ FOR THE PUMP
C	ATMI-	ACTIVE TRANSPORT ACROSS THE LUMINAL MEMBRANE
C	ATIS-	ACTIVE TRANSPORT ACROSS THE PERITUBULAR MEMBRANE
C	NPHK-	ABUNDANCE OF H-K-ATPASE
C	NAE1-	ABUNDANCE OF AE1
C	NTSC-	ABUNDANCE OF TSC
C	NNHE3-	ABUNDANCE OF NHE3
C	NNKCC-	ABUNDANCE OF NKCC
C	NKCL-	ABUNDANCE OF KCC
C
C	JNAK-	FLUX (MMOL/S) ACROSS THE PERITUBULAR NA,K-ATPASE
C	JHK-	FLUX ACROSS THE LUMINAL H,K-ATPASE
C	JHP-	FLUX ACROSS THE H-ATPASE
C	JAE1-	FLUX ACROSS AE1
C	JTSC-	FLUX ACROSS TSC
C	JNHE3-	FLUX ACROSS NHE3
C	JNKCC-	FLUX ACROSS NKCC
C	JKCC-	FLUX ACROSS KCC
C
C       FEVM-   JUCTION VOLUME FLOW (ML/SEC)
C       FEKM-   JUNCTION SOLUTE FLUX (MMOL/SEC)
C       FIVM-   CELL APICAL VOLUME FLOW
C       FIKM-   CELL APICAL SOLUTE FLUX
C       FEVS-   CHANNEL BASEMENT MEMBRANE VOLUME FLOW
C       FEKS-   CHANNEL BASEMENT MEMBRANE SOLUTE FLUX
C       FIVS-   BASAL CELL VOLUME FLOW
C       FIKS-   BASAL CELL SOLUTE FLUX
C       CURE-   CHANNEL CURRENT (MAMP)
C       CURI-   CELL CURRENT
C       JV-     LATERAL CELL MEMBRANE VOLUME FLOW
C       JS-     LATERAL CELL MEMBRANE SOLUTE FLUX
C
C	TL-	TUBULE LENGTH
C	DX-	SPATIAL STEP OF INTEGRATION
C	RM0-	TUBULE RADIUS AT ZERO TRANSTUBULAR PRESSURE
C	MUM-	NUMBER OF TUBULES AT THE START OF THE SEGMENT
C	ETA-	TUBULE FLUID VISCOSITY
C	SM-	LUMINAL CIRCUMFERENCE
C	AM-	LUMINAL AREA
C	FVM-	LUMINAL VOLUME FLOW RATE
C	FKM-	LUMINAL SOLUTE FLUX RATE
C
	DOUBLE PRECISION RJ,TPS,TNS,QPS,QNS,ZIMPS
	DOUBLE PRECISION PHI,GAMMA,DGAM,PPHI,MAXPHI,
     1    DELGAM,NDERIV,RPHI
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       DGAM-   VARIABLE INCREMENT TO COMPUTE NUMERICAL DERIVATIVE
C       PPHI-   ERROR VECTOR AT GAMMA + DGAM
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C
	DIMENSION PHI(80),GAMMA(80),PPHI(80),NDERIV(80,80),
     1    RPHI(80),DELGAM(80)
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(80),PFL(80)
C
C
C
        COMMON SOLS,T,X,CHOP,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,BCO2,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   TL,DX,RM0,MUM,ETA,ZIMPS,
     1   SM,AM,FVM,FKM,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AMI,AIS,AI0,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/ ISOFM,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM,
     1   NNKCC,NKCL,JNKCC,JKCC
C
C IN THE GAMMA ARRAY); IF OFF (0) READS VM.
C
C
	OPEN (80,FILE='ahlmparam.dat')
	OPEN (81,FILE='ahlmresult.dat')
	OPEN (82,FILE='ahlmbound.dat')
	OPEN (83,FILE='ahlmptpick.dat')
	OPEN (84,FILE='ahlmguess.dat')
	OPEN (86,FILE='ahlmptim.dat')
	OPEN (87,FILE='ahlmfluxes.dat')
C
	OPEN (10,FILE='lumen.dat')
	OPEN (11,FILE='blood.dat')
	OPEN (12,FILE='xy.co2')
	OPEN (52,FILE='ahlcbnd.dat')
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
C
C  SOLUTE INDEX:
C	1-	NA+
C	2-	K+
C	3-	CL-
C	4-	HCO3-
C	5-	H2CO3
C	6-	CO2
C	7-	HPO4--
C	8-	H2PO4-
C	9-	UREA
C	10-	NH3
C	11-	NH4+
C	12-	H+
C
C
	SOLS=12
	NN=2*(2+SOLS)+1
C
	PRINT 27
   27   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
	READ *, EXP
	DO 600 COUNT=1,EXP
C
C BOUNDARY VALUES,ACTIVE TRANSPORT,AND THE TIME INCREMENT ARE
C READ. THE TIME STEP WILL CYCLE FROM THIS POINT.
C
	T=0
	TIME=0.D0
   40   T=T+1
C
	CALL AHLMNEWT(COUNT)
C
	IF (TAU.EQ.1) GO TO 570
	GO TO 600
C
C  FOR THE TRANSIENT EXPERIMENT,THE VARIABLES AT T=1 ARE SET;
C  IF NEEDED, THE GUESS FOR T=2 IS GENERATED.
C
  570   CALL AHLMRESET(2)
C
C  IF TIME IS UP,STOP.
  585   TLIM=TLIM-1
	IF (TLIM.GT.0) GO TO 40
	PRINT 590
  590   FORMAT (' TIME IS UP')
C
  600   CONTINUE
C
	WRITE (52,50) DT,FVM(X,T),
     1  VM(X,T),VS(X,T),VS(CHOP+1,T),
     1  PM(X,T),PS(X,T),PS(CHOP+1,T),
     1  IMPM(X,T),IMPS(X,T),IMPS(CHOP+1,T),
     1  (CM(I,X,T),CS(I,X,T),CS(I,CHOP+1,T),I=1,SOLS-1)
   50   FORMAT (2D12.4,/,3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
	STOP
	END
