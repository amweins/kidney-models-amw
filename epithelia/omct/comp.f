	PROGRAM COMP
C COMP:COMPARTMENT MODEL FOR EPITHELIAL FLUXES
C REVISED 12/12/94 (STARTING WITH PROXIMAL TUBULE MODEL IN NAH)
C TO REPRESENT THE CORTICAL COLLECTING DUCT EPITHELIUM
C WILL SUPERCEDE THE STRIETER MODEL BY ADDING CO2, H2CO3, NH3, NH4
C MODEL DERIVED FROM CCT-COMP FOR OMCD (PURE ALPHA CELL)
C
	INTEGER SOLS,TAU,T,TLIM,COUNT,EXP,SW1,NPAR
C
C       SOLS-   NUMBER OF SOLUTES
C       TAU-    INDICATES STEADY STATE (=0) OR TIMED (=1) EXPT.
C       T-      INDEX FOR OLD (=1) OR NEW (=2) TIME STEP
C       COUNT-  NUMBER OF EXPT. BEING SIMULATED
C       EXP-    TOTAL NUMBER OF EXPTS. IN THE RUN
C	NPAR-	NUMBER OF PARAMETERS ICREMENTED IN A TRANSIENT PROBLEM.
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY(4),KDHY(4),L0,L(2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(12),IMPM,LCHM,XM(12),
     1   VS,PS,CS(12),IMPS,LCHS,XS(12)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(2),PE(2),CE(12,2),LCHE,XE(12),
     1   FEVM(2),FEKM(12,2),FEVS(2),FEKS(12,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C
C   SUBSCRIPTS (CELL TYPE,SOLUTE,TIME)
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
C       DUMVM-  DUMMY VARIABLE HOLDING VM (OPEN CIRCUIT)
C
C       NP-     NA+ EXTRUSION OUT THE LATERAL CELL MEMBRANE
C	KNPN-	NA AFFINITY FOR THE PUMP
C	KNPK-	K AFFINITY FOR THE PUMP
C       KNH4-   RELATIVE AFFINITY FOR NH4+ FOR THE PUMP
C	NPHK-	ABUNDANCE OF H-K-ATPASE
C	NAE1-	ABUNDANCE OF AE1
C	ATMI-	ACTIVE TRANSPORT ACROSS THE LUMINAL MEMBRANE
C	ATIS-	ACTIVE TRANSPORT ACROSS THE PERITUBULAR MEMBRANE
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
C
	DOUBLE PRECISION PHI,GAMMA,DGAM,PPHI,MAXPHI,PREMAX,
     1    DELGAM,NDERIV,RPHI
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       DGAM-   VARIABLE INCREMENT TO COMPUTE NUMERICAL DERIVATIVE
C       PPHI-   ERROR VECTOR AT GAMMA + DGAM
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C       PREMAX- MAXPHI FOR THE PREVIOUS NEWTON ITERATE
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C
	DIMENSION PHI(60),GAMMA(60),PPHI(60),NDERIV(60,60),
     1    RPHI(60),DELGAM(60)
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(60),PFL(60)
C
C
        COMMON SOLS,T,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AI0,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
C
C SW1--IF ON (1) GIVES OPEN CIRCUIT CONDITIONS (VM IS A VARIABLE
C IN THE GAMMA ARRAY); IF OFF (0) READS VM.
C
C
	OPEN (20,FILE='param.dat')
	OPEN (21,FILE='result.dat')
	OPEN (22,FILE='bound.dat')
	OPEN (23,FILE='xy.dat')
	OPEN (24,FILE='guess.dat')
	OPEN (26,FILE='ptim.dat')
	OPEN (27,FILE='fluxes.dat')
	OPEN (28,FILE='table.dat')
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
	PRINT 20
   20   FORMAT (' TYPE 1/0 FOR OPEN/CLOSED CIRCUIT:  ',$)
	READ *,SW1
	NN=2*(2+SOLS)+SW1
C
	DGAM=1.D-6
	PTOL=1.D-12
	M1=60
	IERR=0
C
	READ (24,24) (GAMMA(J),J=1,NN)
   24   FORMAT (D30.20)
C
	PRINT 27
   27   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
	READ *, EXP
	DO 600 COUNT=1,EXP
C
	CALL PSET(TAU,TLIM)
	T=1
	CALL GAMSET(1,GAMMA)
C
C
C BOUNDARY VALUES,ACTIVE TRANSPORT,AND THE TIME INCREMENT ARE
C READ. THE TIME STEP WILL CYCLE FROM THIS POINT.
C
	T=0
	TIME=0.D0
   40   T=T+1
	READ(22,50) DT,DUMVM,PM,IMPM,VS,PS,IMPS,
     1  (CM(I),CS(I),I=1,SOLS-1)
   50   FORMAT (D12.4,/,6F8.4,/,(2F14.9))
	LCHM=PKC + DLOG10(CM(4)/CM(5))
	LCHS=PKC + DLOG10(CS(4)/CS(5))
	CM(SOLS)=10.**(-LCHM)
	CS(SOLS)=10.**(-LCHS)
	IF (SW1.EQ.0) VM=DUMVM
	TIME=TIME+DT*FLOAT(T-1)
	RTAU=DBLE(FLOAT(TAU))/DT
	IF (TAU.EQ.1) CALL PTIM(NPAR)
C
C   START THE SPATIAL ITERATION.
C
	PRINT 72, COUNT,TIME
   72   FORMAT(5X,'EXP=',I3,5X,'TIME=',F8.3)
   70   ITER=0
   80   ITER=ITER+1
	PRINT 82, ITER
   82   FORMAT(I5,$)
	IF (ITER.GT.1) PREMAX=MAXPHI
	CALL ERRVEC(PHI)
	MAXPHI=DABS(PHI(1))
	DO 230 J=2,NN
  230   IF (MAXPHI.LT.DABS(PHI(J))) MAXPHI=DABS(PHI(J))
	PRINT 232,MAXPHI
  232   FORMAT ('+',D12.4)
	IF(MAXPHI.LT.EPSI) GO TO 558
C
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
C
	CALL GAMSET(2,GAMMA)
C
	DO  900 K=1,NN
	GAMMA(K)=GAMMA(K)+DGAM
C
	CALL GAMSET(1,GAMMA)
C
	CALL ERRVEC(PPHI)
C
	DO 825 J=1,NN
  825   NDERIV(J,K)=(PPHI(J)-PHI(J))/DGAM
C
	GAMMA(K)=GAMMA(K)-DGAM
  900   CONTINUE
C
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C
  923   DO 924 J=1,NN
  924   RPHI(J)=PHI(J)
C
	CALL LES8(NDERIV,RPHI,PR,PFL,NN,SW,PTOL,DELGAM,M1)
	IF (SW) 930,600,930
C
  930   DO 940 J=1,NN
  940   GAMMA(J)=GAMMA(J)-DELGAM(J)
C
C
C RESET THE GUESSES IN THE NAMES BY WHICH WE KNOW THEM.
C
	CALL GAMSET(1,GAMMA)
C
	IF (ITER.LT.50) GO TO 80
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS')
	GO TO 600
C
C
C THE SPATIAL ITERATION CONVERGED,RESULTS WERE OUTPUT, AND
C THE PROGRAM PICKS UP AT THIS POINT:
C
  558   IF(T.EQ.1) CALL ARESUL
	CALL BRESUL
	CALL CRESUL
	CALL FRESUL
C	IF(EXP.EQ.1) CALL TABLE
	IF (TAU.EQ.1) GO TO 570
C
C  THUS IN THE STEADY STATE CASE WE ARE DONE.
	PRINT 560
  560   FORMAT (' PROBLEM SOLVED FOR THE STEADY STATE')
	GO TO 600
C
C  FOR THE TRANSIENT EXPERIMENT,THE VARIABLES AT T=1 ARE SET;
C  IF NEEDED, THE GUESS FOR T=2 IS GENERATED.
C
  570   CALL RESET
C
C  IF TIME IS UP,STOP.
  585   TLIM=TLIM-1
	IF (TLIM.GT.0) GO TO 40
	PRINT 590
  590   FORMAT (' TIME IS UP')
C
C
  600   CONTINUE
C
	CLOSE (24,STATUS='KEEP')
	OPEN (24,FILE='guess.new')
	WRITE (24,624) (GAMMA(J),J=1,NN)
  624   FORMAT(D30.20)
	STOP
	END
