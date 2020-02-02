	PROGRAM KPCR
C  PROGRAM TO GENERATE KPARAM.DAT FROM KPARAM.TEM
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),PKC,PKF,PKN,PKP,KCO2,
     1   LCHS,TPS,TNS,TFS,QPS,QNS,QFS
C
C SPECIFIED BOUNDARY DATA
C
	INTEGER OMCHOP,IMCHOP,MRCHOP,SLICE,SOLS
C
        DOUBLE PRECISION KEPSI,FVM0(2),
     1   PS0(3),IMPM0(2),IMPS0(3),CM0(15,2),CS0(16,3), 
     1   FVC0(3),PC0(3),IMPC0(3),HCT0(3),CC0(16,3),BCO2(3)
C
C	FVMO(1) and FVM0(2)- SFPCT and JMPCT(2) flows for all SF and JM(2) nephrons
C	PSO- 	SFPCT, JMPCT, and DCT/CNT vascular pressures (uniform)
C	IMPM0- filtered impermeant osmolality
C	IMPS0- vascular impermeant osmolality for SFPCT, JMPCT, and DCT/CNT (mM)
C	CM0-	glomerular filtrate composition
C	CS0-	vascular plasma composition for SFPCT, JMPCT, and DCT/CNT
C
C	FVC0-	capillary plasma flow for individual OMDVR, OIDVR, and MRDVR
C	PC0-	capillary pressures
C	IMPC0- 	capillary protein concentration (gm/dl)
C	HCT0-	hematocrits
C	CC0-	capillary solute composition
C	BCO2-	coefficients for CO2 production
C 
C
C LOAD ARCHIVED VARIABLES FOR TUBULES AND VESSELS WITH BOUNDARY DATA:
C
C Boundary dat, specified and variable, will be set at the top level (kidney) and placed
C   in archive variables.
C                        Lumen (Init)            IS (Init)            IS (End)
C       SFPCT,JMPCT         KBSET            KBSET (Cortex)         KBSET (Cortex)
C       SFPST,JMPST           -              KGAM (OM - init)       KGAM (OM - mid)
C       sDHL,lDHLu            -              KGAM (OM - mid)        KGAM (OM - end)
C       lDHLl,tAHL            -              KGAM (OM - end)        KGAM (IM)
C       AHLm                  -              KGAM (OM - end)        KGAM (MR - end)
C       AHLc                  -              KGAM (MR - end)        KGAM (MR - init)
C       DCT,CNT               -              KBSET (Cortex)         KBSET (Cortex)
C       CCT                   -              KGAM (MR - init)       KGAM (MR - end)
C       OMCT                  -              KGAM (MR - end)        KGAM (OM - end)
C       IMCT                  -              KGAM (OM - end)        KGAM (IM - end)
C
C       OMDVR,OMAVR (IS)    KBSET            KGAM (OM - init)       KGAM (OM - end)
C       OIDVR,OIAVR         KBSET            KGAM (MR - end)        KGAM (OM - end)
C       IMDVR,IMAVR           -              KGAM (OM - end)        KGAM (IM)
C       MRDVR,MRAVR (IS)    KBSET            KGAM (MR - init)       KGAM (MR - end)
C
C
C FOR THE BASELINE CONFIRGURATION:
C  SLICE DEMARCATIONS AT 0 - 7 MM FROM OM TO PAPILLA
C  THE INITIAL POINT AT 0 MM CORRESPONDS TO A SINGLE MR SLICE, 2MM IN THICKNESS
C Interstitial unknows are pressure and concentration at 2 points in OM, 5 points
C  in IM, and one point in MR.  The x=0 point in OM abuts PST; there is a second
c  x=0 point in OM that is identified with MR=2mm, and this abuts AHL and OMCD.
C The points MR(1) and OM(1) will be assigned cortical concentrations, so there
C  are a total of 8 sets of IS unknowns.
C
C                                       MR(1)
C
C                                       MR(2)
C                _____________________________________ (continuity MR->OM)
C                        OM(1)          MR(2)
C                             \        /
C                             OM(2)=MR(3)
C                             OM(3)=MR(4)
C                _____________________________________ (continuity OM->IM)
C                                OM(3)
C                                IM(1)
C                                 ...
C                                IM(5)
C
C
        OPEN (302,FILE='kparam.tem')
        OPEN (303,FILE='kparam.dat')
C
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
	Z(1)=1.D0
	Z(2)=1.D0
	Z(3)=-1.D0
	Z(4)=-1.D0
	Z(5)=0.D0
	Z(6)=0.D0
	Z(7)=-2.D0
	Z(8)=-1.D0
	Z(9)=0.D0
	Z(10)=0.D0
	Z(11)=1.D0
	Z(12)=1.D0
	Z(13)=-1.D0
	Z(14)=0.D0
	Z(15)=0.D0
C
	PKC=3.57
	PKF=3.76
	PKP=6.8
	PKN=9.15
        KCO2=340.D0
        SOLS=15
C
C
C PARAMETERS SUBJECT TO VARIATION ARE READ
C
	READ (302,20) 
     1  SLICE,OMCHOP,IMCHOP,MRCHOP,KEPSI,
     1  FVM0(1),FVM0(2),IMPM0(1),IMPM0(2),
     1  PS0(1),PS0(2),PS0(3),
     1  IMPS0(1),IMPS0(2),IMPS0(3),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=1,11),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=13,15)
	READ (302,22) 
     1  FVC0(1),FVC0(2),FVC0(3),
     1  PC0(1),PC0(2),PC0(3),
     1  CC0(16,1),CC0(16,2),CC0(16,3),
     1  HCT0(1),HCT0(2),HCT0(3),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=1,11),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=13,15)
	READ (302,24) (BCO2(I),I=1,3)
   20   FORMAT (4I5,D12.4,/,2D12.4,2F8.4,/,3F8.4,/,3F8.4,/,(5F14.9))
   22   FORMAT (3D12.4,/,3F8.4,/,3F8.4,/,4F8.4,/,(3F14.9))
   24   FORMAT (3F8.4)
C
C
        DO 34 K=1,2
	CM0(5,K)=CM0(6,K)/KCO2
	LCHS = PKC + DLOG10(CM0(4,K)/CM0(5,K))
	CM0(12,K)=10.**(-LCHS)
C
	TPS=CM0(7,K)+CM0(8,K)
	TNS=CM0(10,K)+CM0(11,K)
	TFS=CM0(13,K)+CM0(14,K)
C
	QPS=10.**(LCHS-PKP)
	QNS=10.**(LCHS-PKN)
	QFS=10.**(LCHS-PKF)
C
	CM0(8,K)=TPS/(1.+QPS)
	CM0(7,K)=TPS-CM0(8,K)
	CM0(11,K)=TNS/(1.+QNS)
	CM0(10,K)=TNS-CM0(11,K)
	CM0(14,K)=TFS/(1.+QFS)
	CM0(13,K)=TFS-CM0(14,K)
C
	CM0(3,K)= CM0(1,K)*Z(1) + CM0(2,K)*Z(2)
	DO 33 I=4,SOLS
   33   CM0(3,K)=CM0(3,K)+CM0(I,K)*Z(I)
   34   CONTINUE
C
C
        DO 44 K=1,3
	CS0(5,K)=CS0(6,K)/KCO2
	LCHS = PKC + DLOG10(CS0(4,K)/CS0(5,K))
	CS0(12,K)=10.**(-LCHS)
C
	TPS=CS0(7,K)+CS0(8,K)
	TNS=CS0(10,K)+CS0(11,K)
	TFS=CS0(13,K)+CS0(14,K)
C
	QPS=10.**(LCHS-PKP)
	QNS=10.**(LCHS-PKN)
	QFS=10.**(LCHS-PKF)
C
	CS0(8,K)=TPS/(1.+QPS)
	CS0(7,K)=TPS-CS0(8,K)
	CS0(11,K)=TNS/(1.+QNS)
	CS0(10,K)=TNS-CS0(11,K)
	CS0(14,K)=TFS/(1.+QFS)
	CS0(13,K)=TFS-CS0(14,K)
C
	CS0(3,K)= CS0(1,K)*Z(1) + CS0(2,K)*Z(2)
	DO 43 I=4,SOLS
   43   CS0(3,K)=CS0(3,K)+CS0(I,K)*Z(I)
   44   CONTINUE
C
C
        DO 54 K=1,3
	CC0(5,K)=CC0(6,K)/KCO2
	LCHS = PKC + DLOG10(CC0(4,K)/CC0(5,K))
	CC0(12,K)=10.**(-LCHS)
C
	TPS=CC0(7,K)+CC0(8,K)
	TNS=CC0(10,K)+CC0(11,K)
	TFS=CC0(13,K)+CC0(14,K)
C
	QPS=10.**(LCHS-PKP)
	QNS=10.**(LCHS-PKN)
	QFS=10.**(LCHS-PKF)
C
	CC0(8,K)=TPS/(1.+QPS)
	CC0(7,K)=TPS-CC0(8,K)
	CC0(11,K)=TNS/(1.+QNS)
	CC0(10,K)=TNS-CC0(11,K)
	CC0(14,K)=TFS/(1.+QFS)
	CC0(13,K)=TFS-CC0(14,K)
C
	CC0(3,K)= CC0(1,K)*Z(1) + CC0(2,K)*Z(2)
	DO 53 I=4,SOLS
   53   CC0(3,K)=CC0(3,K)+CC0(I,K)*Z(I)
   54   CONTINUE
C
C
	WRITE (303,20) 
     1  SLICE,OMCHOP,IMCHOP,MRCHOP,KEPSI,
     1  FVM0(1),FVM0(2),IMPM0(1),IMPM0(2),
     1  PS0(1),PS0(2),PS0(3),
     1  IMPS0(1),IMPS0(2),IMPS0(3),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=1,11),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=13,15)
	WRITE (303,22) 
     1  FVC0(1),FVC0(2),FVC0(3),
     1  PC0(1),PC0(2),PC0(3),
     1  CC0(16,1),CC0(16,2),CC0(16,3),
     1  HCT0(1),HCT0(2),HCT0(3),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=1,11),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=13,15)
	WRITE (303,24) (BCO2(I),I=1,3)
C
        STOP
        END
