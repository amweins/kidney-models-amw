        SUBROUTINE GAMSET(SWG,TGGAM)
C PROGRAM TO SET AND RESET THE GAMMA ARRAY OF INDEPENDENT VARIABLES
C FOR SWG=1 THE GAMMA ARE DEFINE; FOR SWG=2, THE NAMED VARIABLES ARE RESET
C GAMSET WILL OPERATE AT THE LEVEL OF THE MAIN PROGRAM, AND DEAL WITH ARCHIVED VARIABLES
C
	INTEGER SWG
C
C KEY TO INDICES:
C  6 nephrons: SF (INEPH = 1) plus 5 JM (INEPH = 2-6)
C    SFPCT and SFPST in series with SDHL
C    JMPCT and JMPST in series with LDHLu and LDHLl
C  LDHLu will be restricted to OM and LDHLl to IM
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
     1   TGGAM(15),DNR
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
	NTUB=6
	GO TO (100,200) SWG
C
  100   INEPH=1
        RFVM(INEPH,1,1,1)=TGGAM(1)*FVM0(1)
        RPM(INEPH,1,1,1)=TGGAM(2)
        DO 120 INEPH=2,NTUB
        RFVM(INEPH,1,1,1)=TGGAM(2*INEPH-1)*FVM0(2)
        RPM(INEPH,1,1,1)=TGGAM(2*INEPH)
  120   CONTINUE
        RPM(1,9,1,1)=TGGAM(13)
C
	RETURN
C
C
  200   INEPH=1
        TGGAM(1)=RFVM(INEPH,1,1,1)/FVM0(1)
        TGGAM(2)=RPM(INEPH,1,1,1)
        DO 220 INEPH=2,NTUB
        TGGAM(2*INEPH-1)=RFVM(INEPH,1,1,1)/FVM0(2)
        TGGAM(2*INEPH)=RPM(INEPH,1,1,1)
  220   CONTINUE
        TGGAM(13)=RPM(1,9,1,1)
C
	RETURN
	END
