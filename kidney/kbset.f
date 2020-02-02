	SUBROUTINE KBSET
C
C  PROGRAM TO LOAD PARAMETERS TO KIDNEY.F, NOTABLLY COMPOSITION OF CORTICAL PLASMA.
C
C LOCAL VARIABLES
C
	INTEGER SOLS,X,CHOP,
     1    OMCHOP,IMCHOP,MRCHOP,SLICE
C
C	X-	INDICATES SPATIAL STEP
C	CHOP-	SPATIAL CHOP OF TUBULE
C       OMCHOP,IMCHOP,MRCHOP - REGIONAL CHOP FOR OM, IM AND MR
C       SLICE-  OM+IM+MR DISCRETIZATION OF GUESSES 
C               (REFINED TO OMCHOP+IMCHOP+MRCHOP BY LINEAR INTERPOLATION)
C
C SPECIFIED BOUNDARY DATA
C
        DOUBLE PRECISION
     1   CSMR(15,201,2),CSOM(15,101,2),CSIM(15,501,2),
     1   PSMR(201,2),PSOM(101,2),PSIM(501,2),
     1   CSGAM(15,51),PSGAM(51),KEPSI,
     1   PS0(3),IMPM0(2),IMPS0(3),CM0(15,2),CS0(16,3), 
     1   FVC0(3),PC0(3),IMPC0(3),HCT0(3),CC0(16,3),BCO2(3)
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
C	BCO2-	coefficients for CO2 production
C 
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
C
C PARAMETERS SUBJECT TO VARIATION ARE READ
C
	READ (302,30) 
     1  SLICE,OMCHOP,IMCHOP,MRCHOP,KEPSI,
     1  FVM0(1),FVM0(2),IMPM0(1),IMPM0(2),
     1  PS0(1),PS0(2),PS0(3),
     1  IMPS0(1),IMPS0(2),IMPS0(3),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=1,11),
     1  (CM0(I,1),CM0(I,2), CS0(I,1),CS0(I,2),CS0(I,3), I=13,15)
	READ (302,32) 
     1  FVC0(1),FVC0(2),FVC0(3),
     1  PC0(1),PC0(2),PC0(3),
     1  CC0(16,1),CC0(16,2),CC0(16,3),
     1  HCT0(1),HCT0(2),HCT0(3),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=1,11),
     1  (CC0(I,1),CC0(I,2),CC0(I,3),I=13,15)
	READ (302,34) (BCO2(I),I=1,3)
   30   FORMAT (4I5,D12.4,/,2D12.4,2F8.4,/,3F8.4,/,3F8.4,/,(5F14.9))
   32   FORMAT (3D12.4,/,3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
   34   FORMAT (3F8.4)
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
C SFPCT
        ISEG=1
        INEPH=1
        RFVM(INEPH,ISEG,1,1)=FVM0(1)
        RIMPM(INEPH,ISEG,1,1)=IMPM0(1)
        RPS(INEPH,ISEG,1,1)=PS0(1)
        RIMPS(INEPH,ISEG,1,1)=IMPS0(1)
        DO 40 ISOL=1,11
        RCM(INEPH,ISEG,ISOL,1,1)=CM0(ISOL,1)
   40   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,1)
        DO 41 ISOL=13,15
        RCM(INEPH,ISEG,ISOL,1,1)=CM0(ISOL,1)
   41   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,1)
C
C JMPCT
        DO 52 INEPH=2,6
        RFVM(INEPH,ISEG,1,1)=FVM0(2)
        RIMPM(INEPH,ISEG,1,1)=IMPM0(2)
        RPS(INEPH,ISEG,1,1)=PS0(2)
        RIMPS(INEPH,ISEG,1,1)=IMPS0(2)
        DO 50 ISOL=1,11
        RCM(INEPH,ISEG,ISOL,1,1)=CM0(ISOL,2)
   50   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,2)
        DO 51 ISOL=13,15
        RCM(INEPH,ISEG,ISOL,1,1)=CM0(ISOL,2)
   51   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,2)
   52   CONTINUE
C
C DCT
        ISEG=8
        DO 62 INEPH=1,6
        RPS(INEPH,ISEG,1,1)=PS0(3)
        RIMPS(INEPH,ISEG,1,1)=IMPS0(3)
        DO 60 ISOL=1,11
   60   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,3)
   62   CONTINUE
C
C CNT
        ISEG=9
        INEPH=1
        RPS(INEPH,ISEG,1,1)=PS0(3)
        RIMPS(INEPH,ISEG,1,1)=IMPS0(3)
        DO 70 ISOL=1,11
   70   RCS(INEPH,ISEG,ISOL,1,1)=CS0(ISOL,3)
C
C OMDVR
        ISEG=1
        DO 102 IVR=1,2
        RFVC(IVR,ISEG,1,1)=FVC0(1)
        RPC(IVR,ISEG,1,1)=PC0(1)
        RIMPC(IVR,ISEG,1,1)=IMPC0(1)
        RHCTC(IVR,ISEG,1,1)=HCT0(1)
        DO 100 ISOL=1,11
  100   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,1)
        DO 101 ISOL=13,16
  101   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,1)
  102   CONTINUE
C
C OIDVR
        DO 112 IVR=3,7
        RFVC(IVR,ISEG,1,1)=FVC0(2)
        RPC(IVR,ISEG,1,1)=PC0(2)
        RIMPC(IVR,ISEG,1,1)=IMPC0(2)
        RHCTC(IVR,ISEG,1,1)=HCT0(2)
        DO 110 ISOL=1,11
  110   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,2)
        DO 111 ISOL=13,16
  111   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,2)
  112   CONTINUE
C
C MRDVR
        IVR=8
        RFVC(IVR,ISEG,1,1)=FVC0(3)
        RPC(IVR,ISEG,1,1)=PC0(3)
        RIMPC(IVR,ISEG,1,1)=IMPC0(3)
        RHCTC(IVR,ISEG,1,1)=HCT0(3)
        DO 120 ISOL=1,11
  120   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,3)
        DO 121 ISOL=13,16
  121   RCC(IVR,ISEG,ISOL,1,1)=CC0(ISOL,3)
C
C Uniform CO2 Production
        DO 140 INEPH=1,6
        DO 140 ISEG=1,12
        DO 140 IPUMP=1,3
  140   RBCO2(INEPH,ISEG,IPUMP)=BCO2(IPUMP)
C
        RETURN
        END
