	SUBROUTINE VASC(COUNT,IFL)
C
C PROGRAM TO CONSOLIDATE THE VR INTO FULL MEDULLARY VASCULATURE 
C Plan for 36000 OMDVR = 24000 sDVR + 12000 lDVR ~ about 3x JM glomeruli
C           7200 MRDVR (DVR from 24000 SF glom = 72000, so this is 10%)  
C Plan for 72000 OMAVR, 48000 IMAVR, 14400 MRAVR; AVR are 2x the DVR number
C  OM and IM DVR derive from 12000 JM nephrons; MR DVR from 2400 SF nephrons
C  Of the IM vessels, turns parallel those of the nephron
C         IMDVR       IMAVR
C    1 mm  6000       12000
C    2 mm  3000        6000
C    3 mm  1500        3000
C    4 mm   750        1500
C    5 mm   375         750
C
C PROGRAM TO BUILD A VAS RECTUM:
C   OUTER MEDULLARY 1,2, descending 1 or 2 mm; 
C   INNER MEDULLARY 3,7, descending 1,...,5 mm into IM;
C   AND MEDULLARY RAY, 8, 2 mm in length.
C  NOTE THAT OMDVR FEEDING IMDVR IS DIFFERENT FROM OMDVR WHICH TURNS IN OM
C
C ARCHIVAL VARIABLES IN ANTICIPATION OF MULTIPLE VESSELS
C THERE ARE TWO MEASURES OF DISTANCE ALONG THE VESSEL:
C    X IS AN INTEGER VARIABLE, AND WILL NOW BE ZEROED WITH EACH SEGMENT
C    DIST IS A REAL VARIABLE THAT TRACKS DISTANCE ALONG THE NEPHRON
C
C SPECIFIED BOUNDARY DATA
C
	INTEGER OMCHOP,IMCHOP,MRCHOP,SLICE
        DOUBLE PRECISION
     1   CSMR(15,201,2),CSOM(15,101,2),CSIM(15,501,2),
     1   PSMR(201,2),PSOM(101,2),PSIM(501,2),
     1   CSGAM(15,51),PSGAM(51),KEPSI,
     1   PS0(3),IMPM0(2),IMPS0(3),CM0(15,2),CS0(16,3), 
     1   FVC0(3),PC0(3),IMPC0(3),HCT0(3),CC0(16,3)
C
C       CSMR,CSOM,CSIM- CS CONCENTRATIONS AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C       PSMR,PSOM,PSIM- PS VALUES AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C       CSGAM-  GUESS CS CONCENTRATIONS AT SLICE BOUNDARIES
C       PSGAM-  GUESS PS VALUES AT SLICE BOUNDARIES
C       KEPSI-   ERROR TOLERANCE FOR THE KIDNEY ITERATIONS
C
C	FVMO(1) and FVM0(2)- SFPCT and JMPCT(2) flows for all SF and JM(2) nephrons
C	PSO- 	SFPCT, JMPCT, and DCT/CNT vascular pressures (uniform)
C	IMPM0- filtered impermeant osmolality
C	IMPS0- vascular impermeant osmolality for SFPCT, JMPCT, and DCT/CNT (mM)
C	CM0-	glomerular filtrate composition
C	CS0-	vascular plasma composition for SFPCT, JMPCT, and DCT/CNT
C	CSMR-	cortical composition within medullary ray
C
C	FVC0-	capillary plasma flow for individual OMDVR, OIDVR, and MRDVR
C	PC0-	capillary pressures
C	IMPC0- 	capillary protein concentration (gm/dl)
C	HCT0-	hematocrits
C	CC0-	capillary solute composition
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
C
	CHARACTER*2 SWVR
C	SWVR-	DIRECTING THE CALCULATION TO OM, IM, or MR vessels 
        DOUBLE PRECISION RADVR(8),RDNUM(8),RANUM(8)
C       RADVR- RATIO OF AVR NUMBER TO DVR NUMBER FOR EACH VESSEL
C       RDNUM- NUMBER OF EACH CLASS OF DVR
C       RANUM- NUMBER OF EACH CLASS OF AVR
C
C KEY TO INDICES:
C  8 vessels: OM (IVR = 1,2) plus 5 IM (IVR = 3-7), plus 1 MR (IVR = 8)
C
C Segment numbers (ISEG)
C      1- OMDVR, MRDVR, OIDVR
C      2- IMDVR
C      3- IMAVR
C      4- OMAVR, MRAVR, OIAVR
C
	INTEGER IARCH,IVR,ISEG
C  IARCH = 0 (ARCHIVE)   = 1 (RETRIEVE)   = 2 (INITIALIZE THE NEXT SEGMENT)
C
	INTEGER 
     1   RVX(8,4), RVCHOP(8,4)
C PARAMETERS
	DOUBLE PRECISION
     1   RVEPSI(8,4),RVDIST(8,4),RVKHY(8,4),RVKDHY(8,4),
     1   RVDX(8,4),RVETA(8,4),RCL(8,4),RRC0(8,4),RMUC(8,4),
     1   RLPCS(8,4),RSCS(8,4,15),RLCS(8,4,15),RHCS(8,4,15),
     1   RLPCSE(8,4),RSCSE(8,4,15),RLCSE(8,4,15),RHCSE(8,4,15),
     1   RLPCSI(8,4),RSCSI(8,4,15),RLCSI(8,4,15),RHCSI(8,4,15)
C VARIABLES
	DOUBLE PRECISION
     1   RVVS(8,4,901,2),RVPS(8,4,901,2),RVCS(8,4,16,901,2),
     1   RVIMPS(8,4,901,2),RVLCHS(8,4,901),RVXS(8,4,15,901),
     1   RVC(8,4,901,2),RPC(8,4,901,2),RCC(8,4,31,901,2),
     1   RIMPC(8,4,901,2),RLCHC(8,4,901),RXC(8,4,15,901),
     1   RSAT(8,4,901,2),RCCN(8,4,901,2),
     1   RCCX(8,4,901,2),RCCY(8,4,901,2),RPOX(8,4,901,2),
     1   RSC(8,4,901,2),RAC(8,4,901,2),RFVC(8,4,901,2),
     1   RFKC(8,4,31,901,2),RFVCS(8,4,901,2),RFKCS(8,4,16,901,2),
     1   RFVCSE(8,4,901,2),RFKCSE(8,4,16,901,2),RFVCSI(8,4,901,2),
     1   RFKCSI(8,4,16,901,2),RHCTC(8,4,901,2),RFBC(8,4,901,2)
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
	COMMON /ARCHVR/
     1   RVX,RVCHOP,RADVR,RDNUM,RANUM,
     1   RVEPSI,RVDIST,RVKHY,RVKDHY,
     1   RVDX,RVETA,RCL,RRC0,RMUC,
     1   RLPCS,RSCS,RLCS,RHCS,
     1   RLPCSE,RSCSE,RLCSE,RHCSE,
     1   RLPCSI,RSCSI,RLCSI,RHCSI,
     1   RVVS,RVPS,RVCS,RVIMPS,RVLCHS,RVXS,
     1   RVC,RPC,RCC,RIMPC,RLCHC,RXC,
     1   RSAT,RCCN,RCCX,RCCY,RSC,RAC,
     1   RFVC,RFKC,RFVCS,RFKCS,
     1   RFVCSE,RFKCSE,RFVCSI,
     1   RFKCSI,RHCTC,RFBC,RPOX
        COMMON /BOUND/
     1   OMCHOP,IMCHOP,MRCHOP,SLICE,
     1   CSMR,CSOM,CSIM,PSMR,PSOM,PSIM,
     1   TGEPSI,TGGAM,CSGAM,PSGAM,KEPSI,
     1   PS0,IMPM0,IMPS0,CM0,CS0, 
     1   FVC0,PC0,IMPC0,HCT0,CC0
C
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
C
        T=1
	SOLS=15
        RVCHOP(1,1)=OMCHOP/2
        RVCHOP(1,4)=OMCHOP/2
        DO 10 IVR=2,7
        RVCHOP(IVR,1)=OMCHOP        
   10   RVCHOP(IVR,4)=OMCHOP        
        RVCHOP(8,1)=MRCHOP
        RVCHOP(8,4)=MRCHOP
        DO 15 IVR=3,7
        RVCHOP(IVR,2)=(IVR-2)*(IMCHOP/5)
   15   RVCHOP(IVR,3)=(IVR-2)*(IMCHOP/5)
C
        CALL VFLOWS(COUNT,IFL)
        RETURN
C
	END
