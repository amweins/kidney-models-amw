	SUBROUTINE KIS
C
C  Program to start with interstitial variables at the nodes,
C   loaded into archived variables, supplied by kgam.f or updated in the search,
C   and to interpolate SLICE values to OMCHOP+IMCHOP+MRCHOP values of each variable.
C
C Linear interpolation for F(J), with J on the interval (IJ,LJ) looks like  
C         RJ = (J-IJ)/(LJ-IJ)
C    F(J) = (1-RJ)*F(IJ) + RJ*F(LJ)
C A "backward" scheme might be RJ = 1.0
C
C LOCAL VARIABLES
C
	INTEGER SWG,SOLS,X,
     1    CHOP,OMCHOP,IMCHOP,MRCHOP,SLICE
C
C	X-	INDICATES SPATIAL STEP
C	CHOP-	SPATIAL CHOP OF TUBULE
C       OMCHOP,IMCHOP,MRCHOP - REGIONAL CHOP FOR OM, IM AND MR
C       SLICE-  OM+IM+MR DISCRETIZATION OF GUESSES 
C               (REFINED TO OMCHOP+IMCHOP+MRCHOP BY LINEAR INTERPOLATION)
C               GAMMA VECTOR IS BASED ON SLICE
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),PKC,PKF,PKN,PKP,
     1   RJ,LCHS(501),TPS,TNS,TFS,QPS,QNS,QFS
C
C SPECIFIED BOUNDARY DATA
C
        DOUBLE PRECISION
     1   CSMR(15,201,2),CSOM(15,101,2),CSIM(15,501,2),
     1   PSMR(201,2),PSOM(101,2),PSIM(501,2),
     1   CSGAM(15,51),PSGAM(51),KEPSI,
     1   PS0(3),IMPM0(2),IMPS0(3),CM0(15,2),CS0(16,3), 
     1   FVC0(3),PC0(3),IMPC0(3),HCT0(3),CC0(16,3)
C
C       CSMR,CSOM,CSIM- CS CONCENTRATIONS AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C       PSMR,PSOM,PSIM- PS VALUES AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C         2 pathways through OM: OM1 -> OM2 -> OM3  (apply to PST and DHL; OMDVR and OMAVR))
C         and MR2 -> OM2 -> OM3 (apply to AHLm and OMCD; OIDVR and OIAVR) Held in MR variables
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
C ARCHIVED NEPHRON VARIABLES
C
	INTEGER 
     1   RX(6,12), RCHOP(6,12)
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   RDIST(6,12), REPSI(6,12), RL0(6,12),
     1   RL(6,12,901,2), RKHY(6,12,5), RKDHY(6,12,5), RBCO2(6,12,3)
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
     1   TGDGAM(15),TGDELG(15),MAXPHI,RPHI(15),DNR
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
C
C VASC:
C
        DOUBLE PRECISION RADVR(8),RDNUM(8),RANUM(8)
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
C	
	COMMON /ARCH/
     1   RX, RCHOP,
     1   RDIST, REPSI, RL0, RL, RKHY, RKDHY, RBCO2,
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
C	
        COMMON /BOUND/
     1   OMCHOP,IMCHOP,MRCHOP,SLICE,
     1   CSMR,CSOM,CSIM,PSMR,PSOM,PSIM,
     1   TGEPSI,TGGAM,CSGAM,PSGAM,KEPSI,
     1   PS0,IMPM0,IMPS0,CM0,CS0, 
     1   FVC0,PC0,IMPC0,HCT0,CC0
C
C FOR THE BASELINE CONFIRGURATION:
C  SLICE DEMARCATIONS AT 0 - 7 MM FROM OM TO PAPILLA
C  THE INITIAL POINT AT 0 MM CORRESPONDS TO A SINGLE MR SLICE, 2MM IN THICKNESS
C Interstitial unknows are pressure and concentration at 3 points in OM, 5 points
C  in IM, and one point in MR.  The x=0 point in OM abuts PST; there is a second
c  x=0 point in OM that is identified with MR=2mm, and this abuts AHL and OMCD.
C
	NTUB=6
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
        SOLS=15
C
C UNLOAD SLICE DATA FROM THE ARCHIVED VARIABLES.
C  After the linear interpolation of the interior points, there is pH correction
C  Note that this is NOT done for the nodal points, as this would erase GAMMA choices
C
C OM-MR
        INEPH=1
        ISEG=11
        KOM=OMCHOP/2
        PSMR(MRCHOP+1,1)=RPS(INEPH,ISEG,1,1)
        PSMR(MRCHOP+KOM+1,1)=RPS(INEPH,ISEG,KOM+1,1)
        PSMR(MRCHOP+OMCHOP+1,1)=RPS(INEPH,ISEG,OMCHOP+1,1)
        DO 240 ISOL=1,SOLS
        CSMR(ISOL,MRCHOP+1,1)=RCS(INEPH,ISEG,ISOL,1,1)
        CSMR(ISOL,MRCHOP+KOM+1,1)=RCS(INEPH,ISEG,ISOL,KOM+1,1)
  240   CSMR(ISOL,MRCHOP+OMCHOP+1,1)=RCS(INEPH,ISEG,ISOL,OMCHOP+1,1)
C
        IJ=MRCHOP+1-KOM
        LJ=MRCHOP+1
        DO 245 K=1,2
        IJ=IJ+KOM
        LJ=LJ+KOM
	DO 244 J=IJ+1,LJ-1
C	RJ=DBLE(FLOAT(J-IJ)/FLOAT(KOM))
C        RJ=0.5 + 0.5*RJ
        RJ=1.D0
	PSMR(J,1)=(1.D0-RJ)*PSMR(IJ,1) + RJ*PSMR(LJ,1)
	DO 242 I=1,SOLS
  242   CSMR(I,J,1)=(1.-RJ)*CSMR(I,IJ,1) + RJ*CSMR(I,LJ,1)
	LCHS(J)=PKC + DLOG10(CSMR(4,J,1)/CSMR(5,J,1))
  244   CONTINUE
  245   CONTINUE
C
C OM-OM
        INEPH=1
        KOM=OMCHOP/2
        ISEG=2
        PSOM(1,1)=RPS(INEPH,ISEG,1,1)
        PSOM(KOM+1,1)=RPS(INEPH,ISEG,KOM+1,1)
        ISEG=3
        PSOM(OMCHOP+1,1)=RPS(INEPH,ISEG,KOM+1,1)
        DO 250 ISOL=1,SOLS
        ISEG=2
        CSOM(ISOL,1,1)=RCS(INEPH,ISEG,ISOL,1,1)
        CSOM(ISOL,KOM+1,1)=RCS(INEPH,ISEG,ISOL,KOM+1,1)
        ISEG=3
  250   CSOM(ISOL,OMCHOP+1,1)=RCS(INEPH,ISEG,ISOL,KOM+1,1)
C
        IJ=1-KOM
        LJ=1
        DO 255 K=1,2
        IJ=IJ+KOM
        LJ=LJ+KOM
	DO 254 J=IJ+1,LJ-1
C	RJ=DBLE(FLOAT(J-IJ)/FLOAT(KOM))
C        RJ=0.5 + 0.5*RJ
        RJ=1.D0
	PSOM(J,1)=(1.D0-RJ)*PSOM(IJ,1) + RJ*PSOM(LJ,1)
	DO 252 I=1,SOLS
  252   CSOM(I,J,1)=(1.-RJ)*CSOM(I,IJ,1) + RJ*CSOM(I,LJ,1)
	LCHS(J)=PKC + DLOG10(CSOM(4,J,1)/CSOM(5,J,1))
  254   CONTINUE
  255   CONTINUE
C
C IM
        ISEG=12
        KIM=IMCHOP/5
        PSIM(1,1)=RPS(INEPH,ISEG,1,1)
        DO 269 K=1,5
  269   PSIM(K*KIM+1,1)=RPS(INEPH,ISEG,K*KIM+1,1)
        DO 271 ISOL=1,SOLS
        CSIM(ISOL,1,1)=RCS(INEPH,ISEG,ISOL,1,1)
        DO 270 K=1,5
  270   CSIM(ISOL,K*KIM+1,1)=RCS(INEPH,ISEG,ISOL,K*KIM+1,1)
  271   CONTINUE
C
        IJ=1-KIM
        LJ=1
        DO 278 K=1,5
        IJ=IJ+KIM
        LJ=LJ+KIM
	DO 277 J=IJ+1,LJ-1
C	RJ=DBLE(FLOAT(J-IJ)/FLOAT(KIM))
C        RJ=0.5 + 0.5*RJ
        RJ=1.D0
	PSIM(J,1)=(1.D0-RJ)*PSIM(IJ,1) + RJ*PSIM(LJ,1)
	DO 275 I=1,SOLS
  275   CSIM(I,J,1)=(1.-RJ)*CSIM(I,IJ,1) + RJ*CSIM(I,LJ,1)
	LCHS(J)=PKC + DLOG10(CSIM(4,J,1)/CSIM(5,J,1))
  277   CONTINUE
  278   CONTINUE
C
C MR
        ISEG=10
        PSMR(1,1)=RPS(INEPH,ISEG,1,1)
        PSMR(MRCHOP+1,1)=RPS(INEPH,ISEG,MRCHOP+1,1)
        DO 280 ISOL=1,SOLS
        CSMR(ISOL,1,1)=RCS(INEPH,ISEG,ISOL,1,1)
  280   CSMR(ISOL,MRCHOP+1,1)=RCS(INEPH,ISEG,ISOL,MRCHOP+1,1)
C
        IJ=1
        LJ=1+MRCHOP
	DO 287 J=IJ+1,LJ-1
C	RJ=DBLE(FLOAT(J-IJ)/FLOAT(MRCHOP))
C        RJ=0.5 + 0.5*RJ
        RJ=1.D0
	PSMR(J,1)=(1.D0-RJ)*PSMR(IJ,1) + RJ*PSMR(LJ,1)
	DO 285 I=1,SOLS
  285   CSMR(I,J,1)=(1.-RJ)*CSMR(I,IJ,1) + RJ*CSMR(I,LJ,1)
	LCHS(J)=PKC + DLOG10(CSMR(4,J,1)/CSMR(5,J,1))
  287   CONTINUE
C
C LOAD THESE INTO THE ARCHIVED VARIABLES.
C    NOTE: PST and sDHL are assumed to have chop of KOM for each
C
C SFPST,JMPST
        ISEG=2
        DO 155 INEPH=1,6
        DO 152 J=1,KOM+1
        RPS(INEPH,ISEG,J,1)=PSOM(J,1)
        DO 150 ISOL=1,SOLS
  150   RCS(INEPH,ISEG,ISOL,J,1)=CSOM(ISOL,J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  152   CONTINUE
  155   CONTINUE
C
C sDHL,lDHLu
        ISEG=3
        DO 165 INEPH=1,6
        DO 164 J=1,KOM+1
        RPS(INEPH,ISEG,J,1)=PSOM(KOM+J,1)
        DO 160 ISOL=1,SOLS
  160   RCS(INEPH,ISEG,ISOL,J,1)=CSOM(ISOL,KOM+J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  164   CONTINUE
  165   CONTINUE
C 
C lDHLl
        ISEG=4
        DO 175 INEPH=2,6
        LJ=(INEPH-1)*KIM+1
        DO 174 J=1,LJ
        RPS(INEPH,ISEG,J,1)=PSIM(J,1)
        DO 170 ISOL=1,SOLS
  170   RCS(INEPH,ISEG,ISOL,J,1)=CSIM(ISOL,J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  174   CONTINUE
  175   CONTINUE
C`
C tAHL
        ISEG=5
        DO 145 INEPH=2,6
        LJ=(INEPH-1)*KIM+1
        DO 144 J=1,LJ
        RPS(INEPH,ISEG,J,1)=PSIM(LJ+1-J,1)
        DO 140 ISOL=1,SOLS
  140   RCS(INEPH,ISEG,ISOL,J,1)=CSIM(ISOL,LJ+1-J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  144   CONTINUE
  145   CONTINUE
C
C AHLm (Experiences MR conditions within OM)
        ISEG=6
        DO 185 INEPH=1,6
        DO 184 J=1,OMCHOP+1
        RPS(INEPH,ISEG,J,1)=PSMR(MRCHOP+OMCHOP+2-J,1)
        DO 180 ISOL=1,SOLS
  180   RCS(INEPH,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+OMCHOP+2-J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  184   CONTINUE
  185   CONTINUE
C
C May need to avoid a jump at the last step AHLM(CHOP)
        DO 186 INEPH=1,6
        J=OMCHOP+1
        RPS(INEPH,ISEG,J,1)=PSMR(MRCHOP+OMCHOP+3-J,1)
        DO 187 ISOL=1,SOLS
  187   RCS(INEPH,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+OMCHOP+3-J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  186   CONTINUE
C
C AHLc
        ISEG=7
        DO 195 INEPH=1,6
        DO 194 J=1,MRCHOP+1
        RPS(INEPH,ISEG,J,1)=PSMR(MRCHOP+2-J,1)
        DO 190 ISOL=1,SOLS
  190   RCS(INEPH,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+2-J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  194   CONTINUE
  195   CONTINUE
C
C CCT
        ISEG=10
        INEPH=1
        DO 304 J=1,MRCHOP+1
        RPS(INEPH,ISEG,J,1)=PSMR(J,1)
        DO 300 ISOL=1,SOLS
  300   RCS(INEPH,ISEG,ISOL,J,1)=CSMR(ISOL,J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  304   CONTINUE
C
C OMCT
        ISEG=11
        INEPH=1
        DO 314 J=1,OMCHOP+1
        RPS(INEPH,ISEG,J,1)=PSMR(MRCHOP+J,1)
        DO 310 ISOL=1,SOLS
  310   RCS(INEPH,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  314   CONTINUE
C
C IMCT
        ISEG=12
        INEPH=1
        DO 324 J=1,IMCHOP+1
        RPS(INEPH,ISEG,J,1)=PSIM(J,1)
        DO 320 ISOL=1,SOLS
  320   RCS(INEPH,ISEG,ISOL,J,1)=CSIM(ISOL,J,1)
        RLCHS(INEPH,ISEG,J)=PKC+
     1    DLOG10(RCS(INEPH,ISEG,4,J,1)/RCS(INEPH,ISEG,5,J,1))
	RCS(INEPH,ISEG,12,J,1)=10.**(-RLCHS(INEPH,ISEG,J))
  324   CONTINUE
C
C
C OMDVR
        ISEG=1
        IVR=1
        DO 334 J=1,KOM+1
        RVPS(IVR,ISEG,J,1)=PSOM(J,1)
        DO 330 ISOL=1,SOLS
  330   RVCS(IVR,ISEG,ISOL,J,1)=CSOM(ISOL,J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  334   CONTINUE
C
        IVR=2
        DO 338 J=1,OMCHOP+1
        RVPS(IVR,ISEG,J,1)=PSOM(J,1)
        DO 337 ISOL=1,SOLS
  337   RVCS(IVR,ISEG,ISOL,J,1)=CSOM(ISOL,J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  338   CONTINUE
C
C OIDVR
        ISEG=1
        DO 345 IVR=3,7
        DO 344 J=1,OMCHOP+1
        RVPS(IVR,ISEG,J,1)=PSMR(MRCHOP+J,1)
        DO 340 ISOL=1,SOLS
  340   RVCS(IVR,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  344   CONTINUE
  345   CONTINUE
C
C MRDVR
        ISEG=1
        IVR=8
        DO 355 J=1,MRCHOP+1
        RVPS(IVR,ISEG,J,1)=PSMR(J,1)
        DO 350 ISOL=1,SOLS
  350   RVCS(IVR,ISEG,ISOL,J,1)=CSMR(ISOL,J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  355   CONTINUE
C
C IMDVR
        ISEG=2
        DO 365 IVR=3,7
        LJ=(IVR-2)*KIM+1
        DO 364 J=1,LJ
        RVPS(IVR,ISEG,J,1)=PSIM(J,1)
        DO 360 ISOL=1,SOLS
  360   RVCS(IVR,ISEG,ISOL,J,1)=CSIM(ISOL,J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  364   CONTINUE
  365   CONTINUE
C
C OMAVR
        ISEG=4
        IVR=1
        DO 434 J=1,KOM+1
        RVPS(IVR,ISEG,J,1)=PSOM(KOM+2-J,1)
        DO 430 ISOL=1,SOLS
  430   RVCS(IVR,ISEG,ISOL,J,1)=CSOM(ISOL,KOM+2-J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  434   CONTINUE
C
        IVR=2
        DO 438 J=1,OMCHOP+1
        RVPS(IVR,ISEG,J,1)=PSOM(OMCHOP+2-J,1)
        DO 437 ISOL=1,SOLS
  437   RVCS(IVR,ISEG,ISOL,J,1)=CSOM(ISOL,OMCHOP+2-J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  438   CONTINUE
C
C OIAVR
        ISEG=4
        DO 445 IVR=3,7
        DO 444 J=1,OMCHOP+1
        RVPS(IVR,ISEG,J,1)=PSMR(MRCHOP+OMCHOP+2-J,1)
        DO 440 ISOL=1,SOLS
  440   RVCS(IVR,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+OMCHOP+2-J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1   DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  444   CONTINUE
  445   CONTINUE
C
C MRAVR
        ISEG=4
        IVR=8
        DO 455 J=1,MRCHOP+1
        RVPS(IVR,ISEG,J,1)=PSMR(MRCHOP+2-J,1)
        DO 450 ISOL=1,SOLS
  450   RVCS(IVR,ISEG,ISOL,J,1)=CSMR(ISOL,MRCHOP+2-J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  455   CONTINUE
C
C IMAVR
        ISEG=3
        DO 465 IVR=3,7
        LJ=(IVR-2)*KIM+1
        DO 464 J=1,LJ
        RVPS(IVR,ISEG,J,1)=PSIM(LJ+1-J,1)
        DO 460 ISOL=1,SOLS
  460   RVCS(IVR,ISEG,ISOL,J,1)=CSIM(ISOL,LJ+1-J,1)
        RVLCHS(IVR,ISEG,J)=PKC+
     1    DLOG10(RVCS(IVR,ISEG,4,J,1)/RVCS(IVR,ISEG,5,J,1))
	RVCS(IVR,ISEG,12,J,1)=10.**(-RVLCHS(IVR,ISEG,J))
  464   CONTINUE
  465   CONTINUE
C
        RETURN
        END
