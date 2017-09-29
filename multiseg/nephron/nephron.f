	PROGRAM NEPHRON
C
C 12/29/13  PROGRAM TO MERGE PN WITH DN TO CREATE A FULL NEPHRON
C Novel additions are 
C  Solving the 6 nephrons through dct,
C    then rapid mixing to enter a common CNT (no CO2 hydration/dehydration at the mixing point)
C  Entering CNT hydrostatic pressure will depend upon end-IMCD outlet pressure (nominally 0 mmHg)
C  Since CNT pressure depends upon CNT flow, it depends upon flow delivered by each nephron, 
C    so that nephrons are weakly coupled.
C  Since proximal nephron flux is pressure-dependent, initial PCT flow and pressure are 
C    subject to TGF and matching CNT inlet pressure, in all 13 global variables.
C  
C  Calculation plan: 
C  SFPCT -> SFPST -> SDHL -> AHLm -> AHLc -> DCT
C  JMPCT -> JMPST -> LDHLu -> LDHLl(2-6) -> tAHL(2-6) -> AHLm(2-6) -> AHLc(2-6) -> DCT(2-6)
C  Merged DCT -> CNT -> CCD -> OMCD -> IMCD
C 
C Carbonic Anhydrase
C   Clearly less in PST; there is some in sDHL; largely absent in IM.
C   Model parameters: 1% of full catalysis in PST, sDHL, and lDHLu; absent (0.01%) in lDHLl and tAHL
C
C 12/22/12 - ARCHIVAL VARIABLES ADDED IN ANTICIPATION OF MULTIPLE NEPHRONS
C THERE ARE TWO MEASURES OF DISTANCE ALONG THE NEPHRON:
C    X IS AN INTEGER VARIABLE, AND WILL NOW BE ZEROED WITH EACH SEGMENT
C    DIST IS A REAL VARIABLE THAT TRACKS DISTANCE ALONG THE NEPHRON
C    ALTHOUGH THIS WILL BECOME AMBIGUOUS AS DIFFERENT LENGTHS COALESCE IN CNT
C
C 2/8/13 - TGF ADJUSTS SNGFR ACCORDING TO MACULA DENSA [CL-]i 
C    AND EARLY PT PRESSURE ACCORDING TO DCT PRESSURE.
C
	INTEGER SOLS,TAU,T,TLIM,NPAR,X,CHOP,COUNT,EXP
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
     1   Z(15),RT,RTE,F,DIST,
     1   EPSI,DT,RTAU,TIME,L0,L(1601,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5)
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
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   TQM(81,2),RM(81,2),MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0(3),NP0(3),NNHE30(3),LPMI0(3),LPIS0(3),
     1   HMI0(3,15),HIS0(3,15),LMI0(3,15,15),LIS0(3,15,15),
     1   QIAMM0,MVL,MVD,RMT0
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
C       AME-    AREA OF TIGHT JUNCTION
C       AE-     CHANNEL AREA AT BASE
C       AMI-    LUMINAL CELL AREA
C       AIE-    LATERAL CELL AREA
C       AIS-    BASAL CELL AREA
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
C
C	JNAK-	FLUX (MMOL/S) ACROSS THE PERITUBULAR NA,K-ATPASE
C	JHK-	FLUX ACROSS THE LUMINAL H,K-ATPASE
C	JHP-	FLUX ACROSS THE H-ATPASE
C	JAE1-	FLUX ACROSS AE1
C	JTSC-	FLUX ACROSS TSC
C	JNHE3-	FLUX ACROSS NHE3
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
C KEY TO INDICES:
C  6 nephrons: SF (INEPH = 1) plus 5 JM (INEPH = 2-6)
C    SFPCT and SFPST in series with SDHL
C    JMPCT and JMPST in series with LDHLu and LDHLl
C  LDHLu will be restricted to OM and LDHLl to IM
C
C Plan for 24000 SF nephrons (@30 nl/min -> 720 mul/min) 
C    12400 JM nephrons (@60 nl/min -> 720 mul/min) 
C     Of the 12000 JM, turns can parallel CD coning:
C    1 mm 6400
C    2 mm 3200
C    3 mm 1600
C    4 mm  800
C    5 mm  400
C
C Segment numbers (ISEG)
C      1- PCT
C      2- PST
C      3- sDHL, lDHLu
C      4- lDHLl
C      5- tAHL
C      6- AHLm
C      7- AHLc
C      8- DCT
C      9- CNT
C      10-CCD
C      11-OMCD
C      12-IMCD
C
	INTEGER IARCH,INEPH,ISEG
C  IARCH = 0 (ARCHIVE) 
C        = 1 (RETRIEVE EVERYTHING) 
C        = 2 (INITIALIZE THE NEXT SEGMENT)
C        = 3 (RETRIEVE EVERYTHING EXCEPT X=1 VARIABLES)
C
	INTEGER 
     1   RX(6,12), RCHOP(6,12)
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   RDIST(6,12), REPSI(6,12), RL0(6,12),
     1   RL(6,12,901,2), RKHY(6,12,5), RKDHY(6,12,5)
C
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   RVM(6,12,901,2), RPM(6,12,901,2), RCM(6,12,15,901,2),
     1   RIMPM(6,12,901,2), RLCHM(6,12,901), RXM(6,12,15,901),
     1   RVS(6,12,901,2), RPS(6,12,901,2), RCS(6,12,15,901,2),
     1   RIMPS(6,12,901,2), RLCHS(6,12,901), RXS(6,12,15,901),
     1   RTL(6,12), RDX(6,12), RRM0(6,12),
     1   RMUM(6,12), RETA(6,12), RSM(6,12,901,2),
     1   RAM(6,12,901,2), RFVM(6,12,901,2), RFKM(6,12,16,901,2)
C
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   RAME(6,12), RAE0(6,12), RAE(6,12,901,2),
     1   RMUA(6,12), RCHVL0(6,12), RCHVL(6,12,901,2),
     1   RMUV(6,12), RLPME(6,12), RLPES(6,12),
     1   RSME(6,12,15), RSES(6,12,15), RHME(6,12,15),
     1   RHES(6,12,15), RVE(6,12,901,2), RPE(6,12,901,2),
     1   RCE(6,12,15,901,2), RLCHE(6,12,901), RXE(6,12,15,901),
     1   RFEVM(6,12,901,2), RFEKM(6,12,15,901,2), RFEVS(6,12,901,2),
     1   RFEKS(6,12,15,901,2), RCURE(6,12,901)
C
C CELL PARAMETERS
	DOUBLE PRECISION
     1   RAIE(6,12,3), RAMI(6,12,3), RAIS(6,12,3),
     1   RAI0(6,12,3), RCLVL0(6,12,3), RIMP0(6,12,3),
     1   RCLVL(6,12,3,901,2), RZIMP(6,12,3), RTBUF(6,12,3),
     1   RPKB(6,12,3), RCBUF(6,12,3,901,2), RHCBUF(6,12,3,901,2),
     1   RLPMI(6,12,3), RLPIS(6,12,3), RSMI(6,12,3,15),
     1   RSIS(6,12,3,15), RHMI(6,12,3,15), RHIS(6,12,3,15),
     1   RLMI(6,12,3,15,15), RLIS(6,12,3,15,15), RATMI(6,12,3,15,901),
     1   RATIS(6,12,3,15,901), RATIE(6,12,3,15,901), RVI(6,12,3,901,2),
     1   RPI(6,12,3,901,2), RCI(6,12,3,15,901,2), RIMP(6,12,3,901),
     1   RLCHI(6,12,3,901), RXI(6,12,3,15,901), RFIVM(6,12,3,901,2),
     1   RFIKM(6,12,3,15,901,2), RFIVS(6,12,3,901,2), 
     1   RFIKS(6,12,3,15,901,2), RCURI(6,12,3,901), 
     1   RJV(6,12,3,901,2), RJK(6,12,3,15,901,2)
C
C SPECIAL TRANSPORTERS
	CHARACTER*1 RISOFM(6,12)
	DOUBLE PRECISION
     1   RNP(6,12,3), RKNPN(6,12,3), RKNPK(6,12,3),
     1   RKNH4(6,12,3), RNPHK(6,12,3), RLHP(6,12,3),
     1   RXIHP(6,12,3), RXHP(6,12,3), RNAE1(6,12,3),
     1   RNTSC(6,12,3), RNNHE3(6,12,3), RJNAK(6,12,3,3,901,2),
     1   RJHK(6,12,3,901,2), RJHP(6,12,3,901,2), RJAE1(6,12,3,901,2),
     1   RJTSC(6,12,3,901,2), RJNHE3(6,12,3,3,901,2), RQIAMM(6,12),
     1   RNNKCC(6,12,3), RNKCL(6,12,3), RJNKCC(6,12,3,4,901,2),
     1   RJKCC(6,12,3,3,901,2)
C
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   RTQM(6,12,501,2), RRM(6,12,501,2), RMUR(6,12), RSCALMI(6,12),
     1   RSCALIS(6,12), RVM0(6,12), RAM0(6,12), RTQM0(6,12),
     1   RLHP0(6,12,3), RNP0(6,12,3), RNNHE30(6,12,3), RLPMI0(6,12,3),
     1   RLPIS0(6,12,3), RHMI0(6,12,3,15), RHIS0(6,12,3,15),
     1   RLMI0(6,12,3,15,15), RLIS0(6,12,3,15,15), RQIAMM0(6,12),
     1   RMVL(6,12), RMVD(6,12), RRMT0(6,12)
C
C TGF VARIABLES: STARTING WITH THE BASELINE SNGFR, 
C   IT IS ADJUSTED ACCORDING TO MACULA DENSA CELL [CL-] 
C INITIAL PT PRESSURE IS ADJUSTED ACCORDING TO A DISTAL NEPHRON RESISTANCE
C
	DOUBLE PRECISION 
     1   FVM0(2),AFVM(6),DFVM(6),AFVM0(2),DFVM0(2),MDCL(2),DMDCL(2),
     1   TGGAM(15),TGPHI(15),TGPPHI(15),TGEPSI,NDERIV(15,15),
     1   TGDGAM(15),TGDELG(15),PREMAX,MAXPHI,RPHI(15),DNR,RDERIV(15,15)
C	
C	FVM0-	SNGFR READ FROM BOUND.DAT-INDEX FOR SF (1) OR JM (2)
C	AFVM-	FIXED COMPONENT OF FVM0
C	DFVM-	COMPONENT OF FVM0 THAT CAN BE MODULATED BY TGF
C	MDCL-	MACULA DENSA CL- REFERENCE CONCENTRATION
C	DMDCL-	MACULA DENSA CL- RANGE YIELDING DFVM
C	TGGAM-	TG FEEDBACK VARIABLE (CHANGE IN FV)
C               1-      SF FVM(1)/FVM0(1)
C               2-      SF PM(1)
C               2J-1-   JM FVM(1)/FVM0(2), J=2,6
C               2J-     JM PM(1)
C               13-     CNT PM(1)
C	TGPHI-	TG FEEDBACK RELATION
C               2J-1-   TGF RELATIONS, J=1,6
C         TGPHI(2J-1) = FVM(0) - AFVM - (DFVM / DMDCL) * [MDCL - CIAHL]
C               2J-     MATCHING PRESSURES AT CNT ENTRY 
C         TGPHI(2J) = PM-END-DCT = PM-CNT
C  xxx          13-     IMCD PM(L)=0 mmHg   - This did not work
C         TGPHI(13) = CNT PM(1) - DNR*FVMCNT
C	TGEPSI-	TG FEEDBACK ERROR
C	NDERIV-	TG FEEDBACK JACOBIAN 
C	TGDGAM-	TG FEEDBACK VARIABLE INCREMENT
C	TGDELG-	TG FEEDBACK VARIABLE NEWTON CORRECTION
C       PREMAX- PRIOR MAXPHI TO ASSESS REUSE OF THE JACOBIAN
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(15),PFL(15)
C
C	
        COMMON SOLS,T,X,CHOP,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,L0,L,
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
	COMMON/KINET/ISOFM,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM,
     1   NNKCC,NKCL,JNKCC,JKCC
	COMMON/TORQUE/
     1   TQM,RM,MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0,NP0,NNHE30,LPMI0,LPIS0,
     1   HMI0,HIS0,LMI0,LIS0,QIAMM0,MVL,MVD,RMT0
C
C
	COMMON /ARCH/
     1   RX, RCHOP,
     1   RDIST, REPSI, RL0, RL, RKHY, RKDHY,
     1   RVM, RPM, RCM, RIMPM, RLCHM, RXM,
     1   RVS, RPS, RCS, RIMPS, RLCHS, RXS,
     1   RTL, RDX, RRM0, RMUM, RETA, RSM, RAM,
     1   RFVM, RFKM,
     1   RAME, RAE0, RAE, RMUA, RCHVL0, RCHVL,
     1   RMUV, RLPME, RLPES, RSME, RSES,
     1   RHME, RHES, RVE, RPE,
     1   RCE, RLCHE, RXE, RFEVM, RFEKM,
     1   RFEVS, RFEKS, RCURE,
     1   RAIE, RAMI, RAIS, RAI0, RCLVL0, RIMP0,
     1   RCLVL, RZIMP, RTBUF, RPKB, RCBUF, RHCBUF,
     1   RLPMI, RLPIS, RSMI, RSIS, RHMI, RHIS,
     1   RLMI, RLIS, RATMI, RATIS, RATIE, RVI, RPI,
     1   RCI, RIMP, RLCHI, RXI, RFIVM, RFIKM,
     1   RFIVS, RFIKS, RCURI, RJV, RJK, RISOFM,
     1   RNP, RKNPN, RKNPK, RKNH4, RNPHK, RLHP,
     1   RXIHP, RXHP, RNAE1, RNTSC, RNNHE3,
     1   RJNAK, RJHK, RJHP, RJAE1, RJTSC, RJNHE3,
     1   RQIAMM, RNNKCC, RNKCL, RJNKCC, RJKCC,
     1   RTQM, RRM, RMUR, RSCALMI, RSCALIS, 
     1   RVM0, RAM0, RTQM0, RLHP0, RNP0, RNNHE30,
     1   RLPMI0, RLPIS0, RHMI0, RHIS0, RLMI0, 
     1   RLIS0, RQIAMM0, RMVL, RMVD, RRMT0
        COMMON /TGF/
     1   FVM0,AFVM,DFVM,AFVM0,DFVM0,MDCL,DMDCL,DNR
C
C
	OPEN (10,FILE='lumen.dat')
	OPEN (11,FILE='blood.dat')
	OPEN (12,FILE='xy.co2')
	OPEN (13,FILE='tgf.dat')
	OPEN (14,FILE='guess.dat')
	OPEN (15,FILE='guess.new')
	OPEN (17,FILE='torque.dat')
	OPEN (19,FILE='errlog')
	OPEN (18,FILE='xy.dat')
C
C
	OPEN (160,FILE='sfpct/sfpctparam.dat')
	OPEN (161,FILE='sfpct/sfpctresult.dat')
	OPEN (162,FILE='sfpct/sfpctbound.dat')
	OPEN (163,FILE='sfpct/sfpctptpick.dat')
	OPEN (164,FILE='sfpct/sfpctguess.dat')
	OPEN (165,FILE='sfpct/sfpctguess.new')
	OPEN (166,FILE='sfpct/sfpctptim.dat')
	OPEN (167,FILE='sfpct/sfpctfluxes.dat')
	OPEN (168,FILE='sfpct/sfpctxy.qo2')
C
	OPEN (150,FILE='sfpst/sfpstparam.dat')
	OPEN (151,FILE='sfpst/sfpstresult.dat')
	OPEN (152,FILE='sfpst/sfpstbound.dat')
	OPEN (153,FILE='sfpst/sfpstptpick.dat')
	OPEN (154,FILE='sfpst/sfpstguess.dat')
	OPEN (155,FILE='sfpst/sfpstguess.new')
	OPEN (156,FILE='sfpst/sfpstptim.dat')
	OPEN (157,FILE='sfpst/sfpstfluxes.dat')
	OPEN (158,FILE='sfpst/sfpstxy.qo2')
C
	OPEN (140,FILE='jmpct/jmpctparam.dat')
	OPEN (141,FILE='jmpct/jmpctresult.dat')
	OPEN (142,FILE='jmpct/jmpctbound.dat')
	OPEN (143,FILE='jmpct/jmpctptpick.dat')
	OPEN (144,FILE='jmpct/jmpctguess.dat')
	OPEN (145,FILE='jmpct/jmpctguess.new')
	OPEN (146,FILE='jmpct/jmpctptim.dat')
	OPEN (147,FILE='jmpct/jmpctfluxes.dat')
	OPEN (148,FILE='jmpct/jmpctxy.qo2')
C
	OPEN (130,FILE='jmpst/jmpstparam.dat')
	OPEN (131,FILE='jmpst/jmpstresult.dat')
	OPEN (132,FILE='jmpst/jmpstbound.dat')
	OPEN (133,FILE='jmpst/jmpstptpick.dat')
	OPEN (134,FILE='jmpst/jmpstguess.dat')
	OPEN (135,FILE='jmpst/jmpstguess.new')
	OPEN (136,FILE='jmpst/jmpstptim.dat')
	OPEN (137,FILE='jmpst/jmpstfluxes.dat')
	OPEN (138,FILE='jmpst/jmpstxy.qo2')
C
	OPEN (120,FILE='sdhl/sdhlparam.dat')
	OPEN (121,FILE='sdhl/sdhlresult.dat')
	OPEN (122,FILE='sdhl/sdhlbound.dat')
	OPEN (123,FILE='sdhl/sdhlptpick.dat')
	OPEN (124,FILE='sdhl/sdhlguess.dat')
	OPEN (125,FILE='sdhl/sdhlguess.new')
	OPEN (126,FILE='sdhl/sdhlptim.dat')
	OPEN (127,FILE='sdhl/sdhlfluxes.dat')
C
	OPEN (110,FILE='ldhlu/ldhluparam.dat')
	OPEN (111,FILE='ldhlu/ldhluresult.dat')
	OPEN (112,FILE='ldhlu/ldhlubound.dat')
	OPEN (113,FILE='ldhlu/ldhluptpick.dat')
	OPEN (114,FILE='ldhlu/ldhluguess.dat')
	OPEN (115,FILE='ldhlu/ldhluguess.new')
	OPEN (116,FILE='ldhlu/ldhluptim.dat')
	OPEN (117,FILE='ldhlu/ldhlufluxes.dat')
C
	OPEN (100,FILE='ldhll/ldhllparam.dat')
	OPEN (101,FILE='ldhll/ldhllresult.dat')
	OPEN (102,FILE='ldhll/ldhllbound.dat')
	OPEN (103,FILE='ldhll/ldhllptpick.dat')
	OPEN (104,FILE='ldhll/ldhllguess.dat')
	OPEN (105,FILE='ldhll/ldhllguess.new')
	OPEN (106,FILE='ldhll/ldhllptim.dat')
	OPEN (107,FILE='ldhll/ldhllfluxes.dat')
C
	OPEN (90,FILE='tahl/tahlparam.dat')
	OPEN (91,FILE='tahl/tahlresult.dat')
	OPEN (92,FILE='tahl/tahlbound.dat')
	OPEN (93,FILE='tahl/tahlptpick.dat')
	OPEN (94,FILE='tahl/tahlguess.dat')
	OPEN (95,FILE='tahl/tahlguess.new')
	OPEN (96,FILE='tahl/tahlptim.dat')
	OPEN (97,FILE='tahl/tahlfluxes.dat')
C
	OPEN (80,FILE='ahlm/ahlmparam.dat')
	OPEN (81,FILE='ahlm/ahlmresult.dat')
	OPEN (82,FILE='ahlm/ahlmbound.dat')
	OPEN (83,FILE='ahlm/ahlmptpick.dat')
	OPEN (84,FILE='ahlm/ahlmguess.dat')
	OPEN (85,FILE='ahlm/ahlmguess.new')
	OPEN (86,FILE='ahlm/ahlmptim.dat')
	OPEN (87,FILE='ahlm/ahlmfluxes.dat')
C
	OPEN (70,FILE='ahlc/ahlcparam.dat')
	OPEN (71,FILE='ahlc/ahlcresult.dat')
	OPEN (72,FILE='ahlc/ahlcbound.dat')
	OPEN (73,FILE='ahlc/ahlcptpick.dat')
	OPEN (74,FILE='ahlc/ahlcguess.dat')
	OPEN (75,FILE='ahlc/ahlcguess.new')
	OPEN (76,FILE='ahlc/ahlcptim.dat')
	OPEN (77,FILE='ahlc/ahlcfluxes.dat')
C
	OPEN (60,FILE='dct/dcparam.dat')
	OPEN (61,FILE='dct/dcresult.dat')
	OPEN (62,FILE='dct/dcbound.dat')
	OPEN (63,FILE='dct/dcptpick.dat')
	OPEN (66,FILE='dct/dcptim.dat')
	OPEN (64,FILE='dct/dcguess.dat')
	OPEN (65,FILE='dct/dcguess.new')
	OPEN (67,FILE='dct/dcfluxes.dat')
C
	OPEN (50,FILE='cnt/cnparam.dat')
	OPEN (51,FILE='cnt/cnresult.dat')
	OPEN (52,FILE='cnt/cnbound.dat')
	OPEN (53,FILE='cnt/cnptpick.dat')
	OPEN (54,FILE='cnt/cnguess.dat')
	OPEN (55,FILE='cnt/cnguess.new')
	OPEN (56,FILE='cnt/cnptim.dat')
	OPEN (57,FILE='cnt/cnfluxes.dat')
C
	OPEN (40,FILE='cct/ccparam.dat')
	OPEN (41,FILE='cct/ccresult.dat')
	OPEN (42,FILE='cct/ccbound.dat')
	OPEN (43,FILE='cct/ccptpick.dat')
	OPEN (44,FILE='cct/ccguess.dat')
	OPEN (45,FILE='cct/ccguess.new')
	OPEN (46,FILE='cct/ccptim.dat')
	OPEN (47,FILE='cct/ccfluxes.dat')
C
	OPEN (30,FILE='ompct/omparam.dat')
	OPEN (31,FILE='ompct/omresult.dat')
	OPEN (32,FILE='ompct/ombound.dat')
	OPEN (33,FILE='ompct/omptpick.dat')
	OPEN (34,FILE='ompct/omguess.dat')
	OPEN (35,FILE='ompct/omguess.new')
	OPEN (36,FILE='ompct/omptim.dat')
	OPEN (37,FILE='ompct/omfluxes.dat')
C
	OPEN (20,FILE='imct/imparam.dat')
	OPEN (21,FILE='imct/imresult.dat')
	OPEN (22,FILE='imct/imbound.dat')
	OPEN (23,FILE='imct/imptpick.dat')
	OPEN (24,FILE='imct/imguess.dat')
	OPEN (25,FILE='imct/imguess.new')
	OPEN (26,FILE='imct/imptim.dat')
	OPEN (27,FILE='imct/imfluxes.dat')
C
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
C	13-	HCO2-
C	14-	H2CO2
C	15-	GLUC
C
C  IN DHL AND BEYOND, THE SOLUTES ARE TRUNCATED TO THE FIRST 12
C    WITH HCO2- ADDED TO CL- AND GLUC BECOMING IMPM
C
	ZIMPS=0.D0
	ZIMPM=0.D0
	TAU=0
	T=1
	TIME=0.D0
C
        NN=13
	PTOL=1.D-10
	M1=15
	IERR=0
C
        READ(14,20) (TGGAM(J),J=1,NN)
   20   FORMAT(D30.20)
        DO 25 J=1,NN
   25   TGDGAM(J)=0.01*TGGAM(J)
        CALL GAMSET(1,TGGAM)
C
	PRINT 30
   30   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
	READ *, EXP
	DO 600 COUNT=1,EXP
C
C TG FEEDBACK SETTINGS: (Set outside proximal tubule since they involve AHL)
        READ(13,35) TGEPSI,DNR,
     1   (AFVM0(J),DFVM0(J),MDCL(J),DMDCL(J),J=1,2)
   35   FORMAT(2D12.4,/,4D12.4,/,4D12.4)
C
C INSERT THE TGF ITERATIONS HERE
C  CALLS TO THE *NEWT PROGRAMS SHOULD NOT REREAD PARAMETERS OR BOUNDARY DATA
C	
  170   ITER=0
  180   ITER=ITER+1
	PRINT 182, ITER
  182   FORMAT(16X,'TGF ITER=',I5,$)
        IF (ITER.EQ.1) THEN
	CALL FLOWS(COUNT,0)
        ELSE
	CALL FLOWS(COUNT,2)
        ENDIF
        CALL ERRVEC(TGPHI)
	PRINT 184, (TGGAM(K),K=1,NN)
  184   FORMAT (/,8X,'TGGAM =',/,6D18.10,/,6D18.10,/,D18.10)
	PRINT 185, (TGPHI(K),K=1,NN)
  185   FORMAT (/,8X,'TGPHI =',/,6D18.10,/,6D18.10,/,D18.10)
C
        PREMAX=MAXPHI
	MAXPHI=DABS(TGPHI(1))
	MAXJ=1
	DO 230 J=2,NN
	IF (MAXPHI.LT.DABS(TGPHI(J))) MAXJ=J
  230   IF (MAXPHI.LT.DABS(TGPHI(J))) MAXPHI=DABS(TGPHI(J))
	PRINT 232, MAXJ, MAXPHI
  232   FORMAT ('+',I5,D12.4)
	IF(MAXPHI.LT.TGEPSI) GO TO 401
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
	IF (ITER.EQ.1) GO TO 250
        IF (MAXPHI.LE.0.1D0*PREMAX) GO TO 923
  250   CONTINUE
C
	DO  900 K=1,NN
	TGGAM(K)=TGGAM(K)+TGDGAM(K)
C
	CALL GAMSET(1,TGGAM)
C
	CALL FLOWS(COUNT,2)
	CALL ERRVEC(TGPPHI)
C        PRINT 186, K
C	PRINT 184, (TGGAM(KK),KK=1,NN)
C        PRINT 185, (TGPPHI(KK),KK=1,NN)
  186   FORMAT(2X,'K=',I2,$)
C
	DO 825 J=1,NN
  825   NDERIV(J,K)=(TGPPHI(J)-TGPHI(J))/TGDGAM(K)
C
	TGGAM(K)=TGGAM(K)-TGDGAM(K)
  900   CONTINUE
C
C        do 904 j=1,nn
C  904   write (18,905) (nderiv(j,k),k=1,nn)
C  905   format (13d10.2)
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C
  923   DO 924 J=1,NN
  924   RPHI(J)=TGPHI(J)
C
        DO 925 K=1,NN
        DO 925 J=1,NN
  925   RDERIV(J,K)=NDERIV(J,K)
C
	CALL LES8(RDERIV,RPHI,PR,PFL,NN,SW,PTOL,TGDELG,M1)
	IF (SW) 930,600,930
C
  930   DO 940 J=1,NN
  940   TGGAM(J)=TGGAM(J)-TGDELG(J)
C	PRINT 184, (TGGAM(KK),KK=1,NN)
C
C RESET THE GUESSES IN THE NAMES BY WHICH WE KNOW THEM.
C
	CALL GAMSET(1,TGGAM)
C
	IF (ITER.LT.50) GO TO 180
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS IN NEPHRON')
	GO TO 600
C
C THE SPATIAL STEP CONVERGED,AND THE PROGRAM PICKS UP AT THIS POINT:
C
  188   CONTINUE
C
  401   INEPH=1
        SOLS=15
	CALL ARCHIVE(1,INEPH,1)
	DIST= -0.5*TL
	X= 1 - CHOP/2
	DO 411 LL=1,3
	DIST=DIST + TL*0.5
	X= X + CHOP/2
	IF(LL.EQ.1) CALL SFPCTRESA
	CALL SFPCTRESB
	CALL SFPCTRESC
  411   CALL SFPCTRESF
	CALL SFPCTRESD
        CALL SFPCTPICK
	CALL SFPCTORQUE(3)
C
	CALL ARCHIVE(1,INEPH,2)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 412 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL SFPSTRESA
	CALL SFPSTRESB
	CALL SFPSTRESC
  412   CALL SFPSTRESF
	CALL SFPSTRESD
        CALL SFPSTPICK
C
	SOLS = 12
	CALL ARCHIVE(1,INEPH,3)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 413 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL SDHLRESA
	CALL SDHLRESB
  413   CALL SDHLRESC
	CALL SDHLRESF
	CALL SDHLRESD
        CALL SDHLPICK
C
	CALL ARCHIVE(1,INEPH,6)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 414 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL AHLMRESA
	CALL AHLMRESB
  414   CALL AHLMRESC
	CALL AHLMRESF
	CALL AHLMRESD
        CALL AHLMPICK
C
	CALL ARCHIVE(1,INEPH,7)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 415 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL AHLCRESA
	CALL AHLCRESB
  415   CALL AHLCRESC
	CALL AHLCRESF
	CALL AHLCRESD
        CALL AHLCPICK
C
	CALL ARCHIVE(1,INEPH,8)
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 416 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL DCTRESA
	CALL DCTRESB
  416   CALL DCTRESC
        CALL DCTRESD
        CALL DCTPICK
C
C JUXTAMEDULLARY NEPHRON
C
        DO 450 J=2,6
  430   INEPH=J
        SOLS=15
	CALL ARCHIVE(1,INEPH,1)
	DIST= -0.5*TL
	X= 1 - CHOP/2
	DO 431 LL=1,3
	DIST=DIST + TL*0.5
	X= X + CHOP/2
	IF(LL.EQ.1) CALL JMPCTRESA
	CALL JMPCTRESB
	CALL JMPCTRESC
  431   CALL JMPCTRESF
	CALL JMPCTRESD
        CALL JMPCTPICK
	CALL JMPCTORQUE(3)
C
	CALL ARCHIVE(1,INEPH,2)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 432 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL JMPSTRESA
	CALL JMPSTRESB
	CALL JMPSTRESC
  432   CALL JMPSTRESF
	CALL JMPSTRESD
        CALL JMPSTPICK
C
	SOLS = 12
	CALL ARCHIVE(1,INEPH,3)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 433 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL LDHLURESA
	CALL LDHLURESB
  433   CALL LDHLURESC
	CALL LDHLURESF
	CALL LDHLURESD
	CALL LDHLUPICK
C
	SOLS = 12
	CALL ARCHIVE(1,INEPH,4)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 434 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL LDHLLRESA
	CALL LDHLLRESB
  434   CALL LDHLLRESC
	CALL LDHLLRESF
	CALL LDHLLRESD
	CALL LDHLLPICK
C
	SOLS = 12
	CALL ARCHIVE(1,INEPH,5)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 435 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL TAHLRESA
	CALL TAHLRESB
  435   CALL TAHLRESC
	CALL TAHLRESF
	CALL TAHLRESD
	CALL TAHLPICK
C
	CALL ARCHIVE(1,INEPH,6)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 436 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL AHLMRESA
	CALL AHLMRESB
  436   CALL AHLMRESC
	CALL AHLMRESF
	CALL AHLMRESD
	CALL AHLMPICK
C
	CALL ARCHIVE(1,INEPH,7)
	DIST= DIST - 0.5*TL
	X= 1 - CHOP/2
	DO 437 LL=1,3
	DIST=DIST + TL*0.5
	X=X + CHOP/2
	IF(LL.EQ.1) CALL AHLCRESA
	CALL AHLCRESB
  437   CALL AHLCRESC
	CALL AHLCRESF
	CALL AHLCRESD
	CALL AHLCPICK
C
	CALL ARCHIVE(1,INEPH,8)
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 438 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL DCTRESA
	CALL DCTRESB
  438   CALL DCTRESC
        CALL DCTRESD
        CALL DCTPICK
C
  450   CONTINUE
C
        INEPH=1
C
	CALL ARCHIVE(1,INEPH,9)
        X=1
        CALL CNTFLUX
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 461 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL CNTRESA
	CALL CNTRESB
  461   CALL CNTRESC
        CALL CNTRESD
        CALL CNTPICK
C
	CALL ARCHIVE(1,INEPH,10)
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 463 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL CCTRESA
	CALL CCTRESB
  463   CALL CCTRESC
        CALL CCTRESD
        CALL CCTPICK
C
	CALL ARCHIVE(1,INEPH,11)
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 465 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL OMPCTRESA
	CALL OMPCTRESB
  465   CALL OMPCTRESC
        CALL OMPCTRESD
        CALL OMPCTPICK
C
	CALL ARCHIVE(1,INEPH,12)
	DIST=DIST - TL*0.5
	X=1-(CHOP/2)
	DO 467 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL IMCTRESA
	CALL IMCTRESB
  467   CALL IMCTRESC
        CALL IMCTRESD
        CALL IMCTPICK
C
C
C  THUS IN THE STEADY STATE CASE WE ARE DONE.
	PRINT 260, COUNT
  260   FORMAT (' PROBLEM ',I2,' SOLVED FOR THE STEADY STATE')
C
        WRITE (15,20) (TGGAM(J),J=1,NN)
  600   CONTINUE
	STOP
	END
C
