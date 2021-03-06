	SUBROUTINE ARCHIVE(IARCH,INEPH,ISEG)
C
C PROGRAM TO ARCHIVE OR RETRIEVE A SEGMENTAL SOLUTION FOR EACH NEPHRON
C
C KEY TO INDICES:
C  6 nephrons: SF (INEPH = 1) plus 5 JM (INEPH = 2-6)
C    SFPCT and SFPST in series with SDHL
C    JMPCT and JMPST in series with LDHLu and LDHLl
C  LDHLu will be restricted to OM and LDHLl to IM
C
C Plan for 24000 SF nephrons 12000 JM nephrons
C  Of the 12000 JM, turns can parallel CD coning:
C    1 mm 6000
C    2 mm 3000
C    3 mm 1500
C    4 mm  750
C    5 mm  375
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
C  IARCH = 0 (ARCHIVE)   = 1 (RETRIEVE)   = 2 (INITIALIZE THE NEXT SEGMENT)
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
     1   RTQM(6,12,901,2), RRM(6,12,901,2), RMUR(6,12), RSCALMI(6,12),
     1   RSCALIS(6,12), RVM0(6,12), RAM0(6,12), RTQM0(6,12),
     1   RLHP0(6,12,3), RNP0(6,12,3), RNNHE30(6,12,3), RLPMI0(6,12,3),
     1   RLPIS0(6,12,3), RHMI0(6,12,3,15), RHIS0(6,12,3,15),
     1   RLMI0(6,12,3,15,15), RLIS0(6,12,3,15,15), RQIAMM0(6,12),
     1   RMVL(6,12), RMVD(6,12), RRMT0(6,12)
C
C MODEL VARIABLES
C
	INTEGER SOLS,TAU,T,TLIM,COUNT,EXP,NPAR,X,CHOP
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
	COMMON/KINET/ ISOFM,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM,
     1   NNKCC,NKCL,JNKCC,JKCC
	COMMON/TORQUE/
     1   TQM,RM,MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0,NP0,NNHE30,LPMI0,LPIS0,
     1   HMI0,HIS0,LMI0,LIS0,QIAMM0,MVL,MVD,RMT0
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
C
	GOTO (100,200,300) IARCH+1
C
  100   CONTINUE
C
C GENERAL PARAMETERS
	RCHOP(INEPH,ISEG) = CHOP
	REPSI(INEPH,ISEG) = EPSI
	RL0(INEPH,ISEG) = L0
	RTL(INEPH,ISEG) = TL
	RDX(INEPH,ISEG) = DX
	RRM0(INEPH,ISEG) = RM0
	RMUM(INEPH,ISEG) = MUM
	RETA(INEPH,ISEG) = ETA
	RAME(INEPH,ISEG) = AME
	RAE0(INEPH,ISEG) = AE0
	RMUA(INEPH,ISEG) = MUA
	RCHVL0(INEPH,ISEG) = CHVL0
	RMUV(INEPH,ISEG) = MUV
	RLPME(INEPH,ISEG) = LPME
	RLPES(INEPH,ISEG) = LPES
C
	DO 110 ISITE = 1,5
	RKHY(INEPH,ISEG,ISITE) = KHY(ISITE)
	RKDHY(INEPH,ISEG,ISITE) = KDHY(ISITE)
  110   CONTINUE
C
C LUMINAL AND PERITUBULAR PARAMETERS
C
	DO 120 IX = 1,CHOP+1
	RL(INEPH,ISEG,IX,T) = L(IX,T)
	RVM(INEPH,ISEG,IX,T) = VM(IX,T)
	RPM(INEPH,ISEG,IX,T) = PM(IX,T)
	RIMPM(INEPH,ISEG,IX,T) = IMPM(IX,T)
	RLCHM(INEPH,ISEG,IX) = LCHM(IX)
	RVS(INEPH,ISEG,IX,T) = VS(IX,T)
	RPS(INEPH,ISEG,IX,T) = PS(IX,T)
	RIMPS(INEPH,ISEG,IX,T) = IMPS(IX,T)
	RLCHS(INEPH,ISEG,IX) = LCHS(IX)
	RSM(INEPH,ISEG,IX,T) = SM(IX,T)
	RAM(INEPH,ISEG,IX,T) = AM(IX,T)
	RFVM(INEPH,ISEG,IX,T) = FVM(IX,T)
	DO 119 ISOL = 1,SOLS
	RCM(INEPH,ISEG,ISOL,IX,T) = CM(ISOL,IX,T)
	RXM(INEPH,ISEG,ISOL,IX) = XM(ISOL,IX)
	RCS(INEPH,ISEG,ISOL,IX,T) = CS(ISOL,IX,T)
	RXS(INEPH,ISEG,ISOL,IX) = XS(ISOL,IX)
	RFKM(INEPH,ISEG,ISOL,IX,T) = FKM(ISOL,IX,T)
  119   CONTINUE
	RFKM(INEPH,ISEG,SOLS+1,IX,T) = FKM(SOLS+1,IX,T)
  120   CONTINUE
C
C INTERSPACE PARAMETERS
	DO 125 ISOL = 1,SOLS
	RSME(INEPH,ISEG,ISOL) = SME(ISOL)
	RSES(INEPH,ISEG,ISOL) = SES(ISOL)
	RHME(INEPH,ISEG,ISOL) = HME(ISOL)
	RHES(INEPH,ISEG,ISOL) = HES(ISOL)
  125   CONTINUE
	DO 130 IX = 1,CHOP+1
	RAE(INEPH,ISEG,IX,T) = AE(IX,T)
	RCHVL(INEPH,ISEG,IX,T) = CHVL(IX,T)
	RVE(INEPH,ISEG,IX,T) = VE(IX,T)
	RPE(INEPH,ISEG,IX,T) = PE(IX,T)
	RLCHE(INEPH,ISEG,IX) = LCHE(IX)
	RFEVM(INEPH,ISEG,IX,T) = FEVM(IX,T)
	RFEVS(INEPH,ISEG,IX,T) = FEVS(IX,T)
	RCURE(INEPH,ISEG,IX) = CURE(IX)
	DO 129 ISOL = 1,SOLS
	RXE(INEPH,ISEG,ISOL,IX) = XE(ISOL,IX)
	RCE(INEPH,ISEG,ISOL,IX,T) = CE(ISOL,IX,T)
	RFEKM(INEPH,ISEG,ISOL,IX,T) = FEKM(ISOL,IX,T)
	RFEKS(INEPH,ISEG,ISOL,IX,T) = FEKS(ISOL,IX,T)
  129   CONTINUE
  130   CONTINUE
C
C CELL PARAMETERS
	DO 141 ICELL = 1,3
	RAIE(INEPH,ISEG,ICELL) = AIE(ICELL)
	RAMI(INEPH,ISEG,ICELL) = AMI(ICELL)
	RAIS(INEPH,ISEG,ICELL) = AIS(ICELL)
	RAI0(INEPH,ISEG,ICELL) = AI0(ICELL)
	RCLVL0(INEPH,ISEG,ICELL) = CLVL0(ICELL)
	RIMP0(INEPH,ISEG,ICELL) = IMP0(ICELL)
	RZIMP(INEPH,ISEG,ICELL) = ZIMP(ICELL)
	RTBUF(INEPH,ISEG,ICELL) = TBUF(ICELL)
	RPKB(INEPH,ISEG,ICELL) = PKB(ICELL)
	RLPMI(INEPH,ISEG,ICELL) = LPMI(ICELL)
	RLPIS(INEPH,ISEG,ICELL) = LPIS(ICELL)
	DO 135 ISOL = 1,SOLS
	RSMI(INEPH,ISEG,ICELL,ISOL) = SMI(ICELL,ISOL)
	RSIS(INEPH,ISEG,ICELL,ISOL) = SIS(ICELL,ISOL)
	RHMI(INEPH,ISEG,ICELL,ISOL) = HMI(ICELL,ISOL)
	RHIS(INEPH,ISEG,ICELL,ISOL) = HIS(ICELL,ISOL)
	DO 134 JSOL = 1,SOLS
	RLMI(INEPH,ISEG,ICELL,ISOL,JSOL) = LMI(ICELL,ISOL,JSOL)
	RLIS(INEPH,ISEG,ICELL,ISOL,JSOL) = LIS(ICELL,ISOL,JSOL)
  134   CONTINUE
  135   CONTINUE
	DO 140 IX = 1,CHOP+1
	RCLVL(INEPH,ISEG,ICELL,IX,T) = CLVL(ICELL,IX,T)
	RCBUF(INEPH,ISEG,ICELL,IX,T) = CBUF(ICELL,IX,T)
	RHCBUF(INEPH,ISEG,ICELL,IX,T) = HCBUF(ICELL,IX,T)
	RVI(INEPH,ISEG,ICELL,IX,T) = VI(ICELL,IX,T)
	RPI(INEPH,ISEG,ICELL,IX,T) = PI(ICELL,IX,T)
	RIMP(INEPH,ISEG,ICELL,IX) = IMP(ICELL,IX)
	RLCHI(INEPH,ISEG,ICELL,IX) = LCHI(ICELL,IX)
	RFIVM(INEPH,ISEG,ICELL,IX,T) = FIVM(ICELL,IX,T)
	RFIVS(INEPH,ISEG,ICELL,IX,T) = FIVS(ICELL,IX,T)
	RJV(INEPH,ISEG,ICELL,IX,T) = JV(ICELL,IX,T)
	RCURI(INEPH,ISEG,ICELL,IX) = CURI(ICELL,IX)
	DO 139 ISOL = 1,SOLS
	RCI(INEPH,ISEG,ICELL,ISOL,IX,T) = CI(ICELL,ISOL,IX,T)
	RXI(INEPH,ISEG,ICELL,ISOL,IX) = XI(ICELL,ISOL,IX)
	RFIKM(INEPH,ISEG,ICELL,ISOL,IX,T) = FIKM(ICELL,ISOL,IX,T)
	RFIKS(INEPH,ISEG,ICELL,ISOL,IX,T) = FIKS(ICELL,ISOL,IX,T)
	RJK(INEPH,ISEG,ICELL,ISOL,IX,T) = JK(ICELL,ISOL,IX,T)
	RATMI(INEPH,ISEG,ICELL,ISOL,IX) = ATMI(ICELL,ISOL,IX)
	RATIS(INEPH,ISEG,ICELL,ISOL,IX) = ATIS(ICELL,ISOL,IX)
	RATIE(INEPH,ISEG,ICELL,ISOL,IX) = ATIE(ICELL,ISOL,IX)
  139   CONTINUE
  140   CONTINUE
  141   CONTINUE
C
C SPECIAL TRANSPORTERS
	RISOFM(INEPH,ISEG) = ISOFM
	RQIAMM(INEPH,ISEG) = QIAMM
	DO 151 ICELL = 1,3
	RNP(INEPH,ISEG,ICELL) = NP(ICELL)
	RKNPN(INEPH,ISEG,ICELL) = KNPN(ICELL)
	RKNPK(INEPH,ISEG,ICELL) = KNPK(ICELL)
	RKNH4(INEPH,ISEG,ICELL) = KNH4(ICELL)
	RNPHK(INEPH,ISEG,ICELL) = NPHK(ICELL)
	RLHP(INEPH,ISEG,ICELL) = LHP(ICELL)
	RXIHP(INEPH,ISEG,ICELL) = XIHP(ICELL)
	RXHP(INEPH,ISEG,ICELL) = XHP(ICELL)
	RNAE1(INEPH,ISEG,ICELL) = NAE1(ICELL)
	RNTSC(INEPH,ISEG,ICELL) = NTSC(ICELL)
	RNNHE3(INEPH,ISEG,ICELL) = NNHE3(ICELL)
	RNNKCC(INEPH,ISEG,ICELL) = NNKCC(ICELL)
	RNKCL(INEPH,ISEG,ICELL) = NKCL(ICELL)
	DO 150 IX = 1,CHOP+1
	RJHK(INEPH,ISEG,ICELL,IX,T) = JHK(ICELL,IX,T)
	RJHP(INEPH,ISEG,ICELL,IX,T) = JHP(ICELL,IX,T)
	RJAE1(INEPH,ISEG,ICELL,IX,T) = JAE1(ICELL,IX,T)
	RJTSC(INEPH,ISEG,ICELL,IX,T) = JTSC(ICELL,IX,T)
	DO 148 ISPEC = 1,3
	RJNAK(INEPH,ISEG,ICELL,ISPEC,IX,T) = JNAK(ICELL,ISPEC,IX,T)
	RJNHE3(INEPH,ISEG,ICELL,ISPEC,IX,T) = JNHE3(ICELL,ISPEC,IX,T)
	RJKCC(INEPH,ISEG,ICELL,ISPEC,IX,T) = JKCC(ICELL,ISPEC,IX,T)
  148   CONTINUE
	DO 149 ISPEC = 1,4
	RJNKCC(INEPH,ISEG,ICELL,ISPEC,IX,T) = JNKCC(ICELL,ISPEC,IX,T)
  149   CONTINUE
  150   CONTINUE
  151   CONTINUE
C
C TORQUE AND COMPLIANCE VARIABLES
	RMUR(INEPH,ISEG) = MUR
	RSCALMI(INEPH,ISEG) = SCALMI
	RSCALIS(INEPH,ISEG) = SCALIS
	RVM0(INEPH,ISEG) = VM0
	RAM0(INEPH,ISEG) = AM0
	RTQM0(INEPH,ISEG) = TQM0
	RQIAMM0(INEPH,ISEG) = QIAMM0
	RMVL(INEPH,ISEG) = MVL
	RMVD(INEPH,ISEG) = MVD
	RRMT0(INEPH,ISEG) = RMT0
	DO 160 IX = 1,CHOP+1
	RTQM(INEPH,ISEG,IX,T) = TQM(IX,T)
	RRM(INEPH,ISEG,IX,T) = RM(IX,T)
  160   CONTINUE
	DO 161 ICELL = 1,3
	RLHP0(INEPH,ISEG,ICELL) = LHP0(ICELL)
	RNP0(INEPH,ISEG,ICELL) = NP0(ICELL)
	RNNHE30(INEPH,ISEG,ICELL) = NNHE30(ICELL)
	RLPMI0(INEPH,ISEG,ICELL) = LPMI0(ICELL)
	RLPIS0(INEPH,ISEG,ICELL) = LPIS0(ICELL)
	DO 165 ISOL = 1,SOLS
	RHMI0(INEPH,ISEG,ICELL,ISOL) = HMI0(ICELL,ISOL)
	RHIS0(INEPH,ISEG,ICELL,ISOL) = HIS0(ICELL,ISOL)
	DO 164 JSOL = 1,SOLS
	RLMI0(INEPH,ISEG,ICELL,ISOL,JSOL) = LMI0(ICELL,ISOL,JSOL)
	RLIS0(INEPH,ISEG,ICELL,ISOL,JSOL) = LIS0(ICELL,ISOL,JSOL)
  164   CONTINUE
  165   CONTINUE
  161   CONTINUE
C
	RETURN
C
  200   CONTINUE
C
C GENERAL PARAMETERS
	CHOP = RCHOP(INEPH,ISEG)
	EPSI = REPSI(INEPH,ISEG)
	L0 = RL0(INEPH,ISEG)
	TL = RTL(INEPH,ISEG)
	DX = RDX(INEPH,ISEG)
	RM0 = RRM0(INEPH,ISEG)
	MUM = RMUM(INEPH,ISEG)
	ETA = RETA(INEPH,ISEG)
	AME = RAME(INEPH,ISEG)
	AE0 = RAE0(INEPH,ISEG)
	MUA = RMUA(INEPH,ISEG)
	CHVL0 = RCHVL0(INEPH,ISEG)
	MUV = RMUV(INEPH,ISEG)
	LPME = RLPME(INEPH,ISEG)
	LPES = RLPES(INEPH,ISEG)
C
	DO 210 ISITE = 1,5
	KHY(ISITE) = RKHY(INEPH,ISEG,ISITE)
	KDHY(ISITE) = RKDHY(INEPH,ISEG,ISITE)
  210   CONTINUE
C
C LUMINAL AND PERITUBULAR PARAMETERS
C
	DO 220 IX = 1,CHOP+1
	L(IX,T) = RL(INEPH,ISEG,IX,T)
	VM(IX,T) = RVM(INEPH,ISEG,IX,T)
	PM(IX,T) = RPM(INEPH,ISEG,IX,T)
	IMPM(IX,T) = RIMPM(INEPH,ISEG,IX,T)
	LCHM(IX) = RLCHM(INEPH,ISEG,IX)
	VS(IX,T) = RVS(INEPH,ISEG,IX,T)
	PS(IX,T) = RPS(INEPH,ISEG,IX,T)
	IMPS(IX,T) = RIMPS(INEPH,ISEG,IX,T)
	LCHS(IX) = RLCHS(INEPH,ISEG,IX)
	SM(IX,T) = RSM(INEPH,ISEG,IX,T)
	AM(IX,T) = RAM(INEPH,ISEG,IX,T)
	FVM(IX,T) = RFVM(INEPH,ISEG,IX,T)
	DO 219 ISOL = 1,SOLS
	CM(ISOL,IX,T) = RCM(INEPH,ISEG,ISOL,IX,T)
	XM(ISOL,IX) = RXM(INEPH,ISEG,ISOL,IX)
	CS(ISOL,IX,T) = RCS(INEPH,ISEG,ISOL,IX,T)
	XS(ISOL,IX) = RXS(INEPH,ISEG,ISOL,IX)
	FKM(ISOL,IX,T) = RFKM(INEPH,ISEG,ISOL,IX,T)
  219   CONTINUE
	FKM(SOLS+1,IX,T) = RFKM(INEPH,ISEG,SOLS+1,IX,T)
  220   CONTINUE
C
C INTERSPACE PARAMETERS
	DO 225 ISOL = 1,SOLS
	SME(ISOL) = RSME(INEPH,ISEG,ISOL)
	SES(ISOL) = RSES(INEPH,ISEG,ISOL)
	HME(ISOL) = RHME(INEPH,ISEG,ISOL)
	HES(ISOL) = RHES(INEPH,ISEG,ISOL)
  225   CONTINUE
	DO 230 IX = 1,CHOP+1
	AE(IX,T) = RAE(INEPH,ISEG,IX,T)
	CHVL(IX,T) = RCHVL(INEPH,ISEG,IX,T)
	VE(IX,T) = RVE(INEPH,ISEG,IX,T)
	PE(IX,T) = RPE(INEPH,ISEG,IX,T)
	LCHE(IX) = RLCHE(INEPH,ISEG,IX)
	FEVM(IX,T) = RFEVM(INEPH,ISEG,IX,T)
	FEVS(IX,T) = RFEVS(INEPH,ISEG,IX,T)
	CURE(IX) = RCURE(INEPH,ISEG,IX)
	DO 229 ISOL = 1,SOLS
	CE(ISOL,IX,T) = RCE(INEPH,ISEG,ISOL,IX,T)
	XE(ISOL,IX) = RXE(INEPH,ISEG,ISOL,IX)
	FEKM(ISOL,IX,T) = RFEKM(INEPH,ISEG,ISOL,IX,T)
	FEKS(ISOL,IX,T) = RFEKS(INEPH,ISEG,ISOL,IX,T)
  229   CONTINUE
  230   CONTINUE
C
C CELL PARAMETERS
	DO 241 ICELL = 1,3
	AIE(ICELL) = RAIE(INEPH,ISEG,ICELL)
	AMI(ICELL) = RAMI(INEPH,ISEG,ICELL)
	AIS(ICELL) = RAIS(INEPH,ISEG,ICELL)
	AI0(ICELL) = RAI0(INEPH,ISEG,ICELL)
	CLVL0(ICELL) = RCLVL0(INEPH,ISEG,ICELL)
	IMP0(ICELL) = RIMP0(INEPH,ISEG,ICELL)
	ZIMP(ICELL) = RZIMP(INEPH,ISEG,ICELL)
	TBUF(ICELL) = RTBUF(INEPH,ISEG,ICELL)
	PKB(ICELL) = RPKB(INEPH,ISEG,ICELL)
	LPMI(ICELL) = RLPMI(INEPH,ISEG,ICELL)
	LPIS(ICELL) = RLPIS(INEPH,ISEG,ICELL)
	DO 235 ISOL = 1,SOLS
	SMI(ICELL,ISOL) = RSMI(INEPH,ISEG,ICELL,ISOL)
	SIS(ICELL,ISOL) = RSIS(INEPH,ISEG,ICELL,ISOL)
	HMI(ICELL,ISOL) = RHMI(INEPH,ISEG,ICELL,ISOL)
	HIS(ICELL,ISOL) = RHIS(INEPH,ISEG,ICELL,ISOL)
	DO 234 JSOL = 1,SOLS
	LMI(ICELL,ISOL,JSOL) = RLMI(INEPH,ISEG,ICELL,ISOL,JSOL)
	LIS(ICELL,ISOL,JSOL) = RLIS(INEPH,ISEG,ICELL,ISOL,JSOL)
  234   CONTINUE
  235   CONTINUE
	DO 240 IX = 1,CHOP+1
	CLVL(ICELL,IX,T) = RCLVL(INEPH,ISEG,ICELL,IX,T)
	CBUF(ICELL,IX,T) = RCBUF(INEPH,ISEG,ICELL,IX,T)
	HCBUF(ICELL,IX,T) = RHCBUF(INEPH,ISEG,ICELL,IX,T)
	VI(ICELL,IX,T) = RVI(INEPH,ISEG,ICELL,IX,T)
	PI(ICELL,IX,T) = RPI(INEPH,ISEG,ICELL,IX,T)
	IMP(ICELL,IX) = RIMP(INEPH,ISEG,ICELL,IX)
	LCHI(ICELL,IX) = RLCHI(INEPH,ISEG,ICELL,IX)
	FIVM(ICELL,IX,T) = RFIVM(INEPH,ISEG,ICELL,IX,T)
	FIVS(ICELL,IX,T) = RFIVS(INEPH,ISEG,ICELL,IX,T)
	JV(ICELL,IX,T) = RJV(INEPH,ISEG,ICELL,IX,T)
	CURI(ICELL,IX) = RCURI(INEPH,ISEG,ICELL,IX)
	DO 239 ISOL = 1,SOLS
	CI(ICELL,ISOL,IX,T) = RCI(INEPH,ISEG,ICELL,ISOL,IX,T)
	XI(ICELL,ISOL,IX) = RXI(INEPH,ISEG,ICELL,ISOL,IX)
	FIKM(ICELL,ISOL,IX,T) = RFIKM(INEPH,ISEG,ICELL,ISOL,IX,T)
	FIKS(ICELL,ISOL,IX,T) = RFIKS(INEPH,ISEG,ICELL,ISOL,IX,T)
	JK(ICELL,ISOL,IX,T) = RJK(INEPH,ISEG,ICELL,ISOL,IX,T)
	ATMI(ICELL,ISOL,IX) = RATMI(INEPH,ISEG,ICELL,ISOL,IX)
	ATIS(ICELL,ISOL,IX) = RATIS(INEPH,ISEG,ICELL,ISOL,IX)
	ATIE(ICELL,ISOL,IX) = RATIE(INEPH,ISEG,ICELL,ISOL,IX)
  239   CONTINUE
  240   CONTINUE
  241   CONTINUE
C
C SPECIAL TRANSPORTERS
	ISOFM = RISOFM(INEPH,ISEG)
	QIAMM = RQIAMM(INEPH,ISEG)
	DO 251 ICELL = 1,3
	NP(ICELL) = RNP(INEPH,ISEG,ICELL)
	KNPN(ICELL) = RKNPN(INEPH,ISEG,ICELL)
	KNPK(ICELL) = RKNPK(INEPH,ISEG,ICELL)
	KNH4(ICELL) = RKNH4(INEPH,ISEG,ICELL)
	NPHK(ICELL) = RNPHK(INEPH,ISEG,ICELL)
	LHP(ICELL) = RLHP(INEPH,ISEG,ICELL)
	XIHP(ICELL) = RXIHP(INEPH,ISEG,ICELL)
	XHP(ICELL) = RXHP(INEPH,ISEG,ICELL)
	NAE1(ICELL) = RNAE1(INEPH,ISEG,ICELL)
	NTSC(ICELL) = RNTSC(INEPH,ISEG,ICELL)
	NNHE3(ICELL) = RNNHE3(INEPH,ISEG,ICELL)
	NNKCC(ICELL) = RNNKCC(INEPH,ISEG,ICELL)
	NKCL(ICELL) = RNKCL(INEPH,ISEG,ICELL)
	DO 250 IX = 1,CHOP+1
	JHK(ICELL,IX,T) = RJHK(INEPH,ISEG,ICELL,IX,T)
	JHP(ICELL,IX,T) = RJHP(INEPH,ISEG,ICELL,IX,T)
	JAE1(ICELL,IX,T) = RJAE1(INEPH,ISEG,ICELL,IX,T)
	JTSC(ICELL,IX,T) = RJTSC(INEPH,ISEG,ICELL,IX,T)
	DO 248 ISPEC = 1,3
	JNAK(ICELL,ISPEC,IX,T) = RJNAK(INEPH,ISEG,ICELL,ISPEC,IX,T)
	JNHE3(ICELL,ISPEC,IX,T) = RJNHE3(INEPH,ISEG,ICELL,ISPEC,IX,T)
	JKCC(ICELL,ISPEC,IX,T) = RJKCC(INEPH,ISEG,ICELL,ISPEC,IX,T)
  248   CONTINUE
	DO 249 ISPEC = 1,4
	JNKCC(ICELL,ISPEC,IX,T) = RJNKCC(INEPH,ISEG,ICELL,ISPEC,IX,T)
  249   CONTINUE
  250   CONTINUE
  251   CONTINUE
C
C TORQUE AND COMPLIANCE VARIABLES
	MUR = RMUR(INEPH,ISEG)
	SCALMI = RSCALMI(INEPH,ISEG)
	SCALIS = RSCALIS(INEPH,ISEG)
	VM0 = RVM0(INEPH,ISEG)
	AM0 = RAM0(INEPH,ISEG)
	TQM0 = RTQM0(INEPH,ISEG)
	QIAMM0 = RQIAMM0(INEPH,ISEG)
	MVL = RMVL(INEPH,ISEG)
	MVD = RMVD(INEPH,ISEG)
	RMT0 = RRMT0(INEPH,ISEG)
	DO 260 IX = 1,CHOP+1
	TQM(IX,T) = RTQM(INEPH,ISEG,IX,T)
	RM(IX,T) = RRM(INEPH,ISEG,IX,T)
  260   CONTINUE
	DO 261 ICELL = 1,3
	LHP0(ICELL) = RLHP0(INEPH,ISEG,ICELL)
	NP0(ICELL) = RNP0(INEPH,ISEG,ICELL)
	NNHE30(ICELL) = RNNHE30(INEPH,ISEG,ICELL)
	LPMI0(ICELL) = RLPMI0(INEPH,ISEG,ICELL)
	LPIS0(ICELL) = RLPIS0(INEPH,ISEG,ICELL)
	DO 265 ISOL = 1,SOLS
	HMI0(ICELL,ISOL) = RHMI0(INEPH,ISEG,ICELL,ISOL)
	HIS0(ICELL,ISOL) = RHIS0(INEPH,ISEG,ICELL,ISOL)
	DO 264 JSOL = 1,SOLS
	LMI0(ICELL,ISOL,JSOL) = RLMI0(INEPH,ISEG,ICELL,ISOL,JSOL)
	LIS0(ICELL,ISOL,JSOL) = RLIS0(INEPH,ISEG,ICELL,ISOL,JSOL)
  264   CONTINUE
  265   CONTINUE
  261   CONTINUE
C
	RETURN
C
C CONTROL PASSES HERE TO INITIALIZE THE NEXT SEGMENT FOR INTEGRATION
C
  300   CONTINUE
	IX = RCHOP(INEPH,ISEG)+1
C
C LUMINAL VARIABLES
	VM(1,T) = RVM(INEPH,ISEG,IX,T)
	PM(1,T) = RPM(INEPH,ISEG,IX,T)
	IMPM(1,T) = RIMPM(INEPH,ISEG,IX,T)
	LCHM(1) = RLCHM(INEPH,ISEG,IX)
	FVM(1,T) = RFVM(INEPH,ISEG,IX,T)
	DO 319 ISOL = 1,SOLS
	CM(ISOL,1,T) = RCM(INEPH,ISEG,ISOL,IX,T)
	FKM(ISOL,1,T) = RFKM(INEPH,ISEG,ISOL,IX,T)
  319   CONTINUE
	FKM(SOLS+1,1,T) = RFKM(INEPH,ISEG,SOLS+1,IX,T)
C
C INTERSPACE PARAMETERS
	AE(1,T) = RAE(INEPH,ISEG,IX,T)
	CHVL(1,T) = RCHVL(INEPH,ISEG,IX,T)
	VE(1,T) = RVE(INEPH,ISEG,IX,T)
	PE(1,T) = RPE(INEPH,ISEG,IX,T)
	LCHE(1) = RLCHE(INEPH,ISEG,IX)
	FEVM(1,T) = RFEVM(INEPH,ISEG,IX,T)
	FEVS(1,T) = RFEVS(INEPH,ISEG,IX,T)
	DO 329 ISOL = 1,SOLS
	CE(ISOL,1,T) = RCE(INEPH,ISEG,ISOL,IX,T)
	FEKM(ISOL,1,T) = RFEKM(INEPH,ISEG,ISOL,IX,T)
	FEKS(ISOL,1,T) = RFEKS(INEPH,ISEG,ISOL,IX,T)
  329   CONTINUE
C
C CELL PARAMETERS
	DO 341 ICELL = 1,3
	CLVL(ICELL,1,T) = RCLVL(INEPH,ISEG,ICELL,IX,T)
	CBUF(ICELL,1,T) = RCBUF(INEPH,ISEG,ICELL,IX,T)
	HCBUF(ICELL,1,T) = RHCBUF(INEPH,ISEG,ICELL,IX,T)
	VI(ICELL,1,T) = RVI(INEPH,ISEG,ICELL,IX,T)
	PI(ICELL,1,T) = RPI(INEPH,ISEG,ICELL,IX,T)
	IMP(ICELL,1) = RIMP(INEPH,ISEG,ICELL,IX)
	LCHI(ICELL,1) = RLCHI(INEPH,ISEG,ICELL,IX)
	FIVM(ICELL,1,T) = RFIVM(INEPH,ISEG,ICELL,IX,T)
	FIVS(ICELL,1,T) = RFIVS(INEPH,ISEG,ICELL,IX,T)
	JV(ICELL,1,T) = RJV(INEPH,ISEG,ICELL,IX,T)
	DO 339 ISOL = 1,SOLS
	CI(ICELL,ISOL,1,T) = RCI(INEPH,ISEG,ICELL,ISOL,IX,T)
	XI(ICELL,ISOL,1) = RXI(INEPH,ISEG,ICELL,ISOL,IX)
	FIKM(ICELL,ISOL,1,T) = RFIKM(INEPH,ISEG,ICELL,ISOL,IX,T)
	FIKS(ICELL,ISOL,1,T) = RFIKS(INEPH,ISEG,ICELL,ISOL,IX,T)
	JK(ICELL,ISOL,1,T) = RJK(INEPH,ISEG,ICELL,ISOL,IX,T)
  339   CONTINUE
  341   CONTINUE
C
	RETURN
	END
