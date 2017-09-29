	PROGRAM DN
C
C AMALGAM OF CCT, OMCT, AND IMCT - 2/18/95
C PLUS DCT AND CNT - 2/24/04
C PLUS AHLM AND AHLC - 11/30/08
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
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY(5),KDHY(5),L0,L(1601,2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(1601,2),PM(1601,2),CM(12,1601,2),
     1   IMPM(1601,2),LCHM(1601),XM(12,1601),
     1   VS(1601,2),PS(1601,2),CS(12,1601,2),
     1   IMPS(1601,2),LCHS(1601),XS(12,1601),ZIMPS,
     1   TL,DX,RM0,MUM,ETA,
     1   SM(1601,2),AM(1601,2),FVM(1601,2),FKM(13,1601,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(1601,2),MUA,CHVL0,CHVL(1601,2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(1601,2),PE(1601,2),CE(12,1601,2),LCHE(1601),XE(12,1601),
     1   FEVM(1601,2),FEKM(12,1601,2),FEVS(1601,2),FEKS(12,1601,2),
     1   CURE(1601)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,1601,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,1601,2),HCBUF(3,1601,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),
     1   ATMI(3,12,1601),ATIS(3,12,1601),ATIE(3,12,1601),
     1   VI(3,1601,2),PI(3,1601,2),CI(3,12,1601,2),IMP(3,1601),
     1   LCHI(3,1601),XI(3,12,1601),
     1   FIVM(3,1601,2),FIKM(3,12,1601,2),FIVS(3,1601,2),
     1   FIKS(3,12,1601,2),CURI(3,1601),JV(3,1601,2),JK(3,12,1601,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,1601,2),JHK(3,1601,2),JHP(3,1601,2),
     1   JAE1(3,1601,2),JTSC(3,1601,2),JNHE3(3,3,1601,2),
     1   NNKCC(3),NKCL(3),JNKCC(3,4,1601,2),JKCC(3,3,1601,2)
C
C   SUBSCRIPTS (CELL TYPE,SOLUTE,DISTANCE,TIME)
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
        COMMON SOLS,T,X,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,ZIMPS,
     1   TL,DX,RM0,MUM,ETA,
     1   SM,AM,FVM,FKM,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AI0,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,
     1   NNKCC,NKCL,JNKCC,JKCC
C
C
	OPEN (10,FILE='lumen.dat')
	OPEN (11,FILE='blood.dat')
	OPEN (12,FILE='xy.co2')
	OPEN (13,FILE='eqlib.dat')
	OPEN (15,FILE='eqlibo.dat')
	OPEN (16,FILE='xy.osm')
	OPEN (19,FILE='errlog')
C
	OPEN (80,FILE='ahlm/ahlmparam.dat')
	OPEN (81,FILE='ahlm/ahlmresult.dat')
	OPEN (82,FILE='ahlm/ahlmbound.dat')
	OPEN (83,FILE='ahlm/ahlmptpick.dat')
	OPEN (84,FILE='ahlm/ahlmguess.dat')
	OPEN (86,FILE='ahlm/ahlmptim.dat')
	OPEN (87,FILE='ahlm/ahlmfluxes.dat')
C
	OPEN (70,FILE='ahlc/ahlcparam.dat')
	OPEN (71,FILE='ahlc/ahlcresult.dat')
	OPEN (72,FILE='ahlc/ahlcbound.dat')
	OPEN (73,FILE='ahlc/ahlcptpick.dat')
	OPEN (74,FILE='ahlc/ahlcguess.dat')
	OPEN (76,FILE='ahlc/ahlcptim.dat')
	OPEN (77,FILE='ahlc/ahlcfluxes.dat')
C
	OPEN (60,FILE='dct/dcparam.dat')
	OPEN (61,FILE='dct/dcresult.dat')
	OPEN (62,FILE='dct/dcbound.dat')
	OPEN (63,FILE='dct/dcptpick.dat')
	OPEN (64,FILE='dct/dcguess.dat')
	OPEN (66,FILE='dct/dcptim.dat')
C
	OPEN (50,FILE='cnt/cnparam.dat')
	OPEN (51,FILE='cnt/cnresult.dat')
	OPEN (52,FILE='cnt/cnbound.dat')
	OPEN (53,FILE='cnt/cnptpick.dat')
	OPEN (54,FILE='cnt/cnguess.dat')
	OPEN (56,FILE='cnt/cnptim.dat')
C
	OPEN (40,FILE='cct/ccparam.dat')
	OPEN (41,FILE='cct/ccresult.dat')
	OPEN (42,FILE='cct/ccbound.dat')
	OPEN (43,FILE='cct/ccptpick.dat')
	OPEN (44,FILE='cct/ccguess.dat')
	OPEN (46,FILE='cct/ccptim.dat')
C
	OPEN (30,FILE='ompct/omparam.dat')
	OPEN (31,FILE='ompct/omresult.dat')
	OPEN (32,FILE='ompct/ombound.dat')
	OPEN (33,FILE='ompct/omptpick.dat')
	OPEN (34,FILE='ompct/omguess.dat')
	OPEN (36,FILE='ompct/omptim.dat')
C
	OPEN (20,FILE='imct/imparam.dat')
	OPEN (21,FILE='imct/imresult.dat')
	OPEN (22,FILE='imct/imbound.dat')
	OPEN (23,FILE='imct/imptpick.dat')
	OPEN (24,FILE='imct/imguess.dat')
	OPEN (26,FILE='imct/imptim.dat')
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
	ZIMPS=0.D0
	PRINT 27
   27   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
	READ *, EXP
	DO 9600 COUNT=1,EXP
C
C BOUNDARY VALUES,ACTIVE TRANSPORT,AND THE TIME INCREMENT ARE
C READ. THE TIME STEP WILL CYCLE FROM THIS POINT.
C
	T=0
	TIME=0.D0
   40   T=T+1
C
	CALL AHLMNEWT(COUNT)
	CALL AHLCNEWT(COUNT)
	CALL DCTNEWT(COUNT)
	CALL CNTNEWT(COUNT)
	CALL CCTNEWT(COUNT)
	CALL OMPCTNEWT(COUNT)
	CALL IMCTNEWT(COUNT)
C
	IF (TAU.EQ.0) THEN
C  THUS IN THE STEADY STATE CASE WE ARE DONE.
	PRINT 9560
 9560   FORMAT (' PROBLEM SOLVED FOR THE STEADY STATE')
	ELSE 
C  OTHERWISE TAKE ANOTHER TIME STEP OR, IF TIME IS UP,STOP.
	TLIM=TLIM-1
	IF (TLIM.GT.0) THEN
	DO 71 LL=1,54
   71   BACKSPACE 20
	DO 72 LL=1,39
   72   BACKSPACE 30
	DO 73 LL=1,99
   73   BACKSPACE 40
	DO 74 LL=1,99
   74   BACKSPACE 50
	DO 75 LL=1,49
   75   BACKSPACE 60
	DO 76 LL=1,47
   76   BACKSPACE 70
	DO 77 LL=1,47
   77   BACKSPACE 80
	GO TO 40
	ELSE
	PRINT 9570
 9570   FORMAT (' TIME IS UP')
	ENDIF
	ENDIF
C
 9600   CONTINUE
	STOP
	END
