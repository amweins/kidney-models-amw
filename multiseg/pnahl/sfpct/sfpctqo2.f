	SUBROUTINE SFPCTQO2
C
	INTEGER SOLS,X,T,CHOP
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
	CHARACTER*5 SOL(15)
C
	DOUBLE PRECISION 
     1    PIKM(3,15),PJK(3,15),PIKS(3,15),
     1    CIKM(3,15),CJK(3,15),CIKS(3,15),
     1    GIKM(3,15),GJK(3,15),GIKS(3,15)
C
C	PIKM-	ACTIVE TRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	PJK-	ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	PIKS-	ACTIVE TRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	CIKM-	COTRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	CJK-	COTRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	CIKS-	COTRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	GIKM-	CHANNEL FLUX OF SOLUTE ACROSS APICAL MEMBRANE
C	GJK-	CHANNEL FLUX OF SOLUTE ACROSS LATERAL MEMBRANE	
C	GIKS-	CHANNEL FLUX OF SOLUTE ACROSS BASAL MEMBRANE
C
	DOUBLE PRECISION FLUXV,FLUXS(15),SCALE
C
C	FLUXV-	TOTAL EPITHELIAL VOLUME FLUX
C	FLUXS-	TOTAL EPITHELIAL SOLUTE FLUX
C	SCALE-  FACTOR TO TRANSLATE FROM MMOL/S*CM2 TO PMOL/MIN*MM
C		=PI*DIAM*(1.0E+09)*60/10
C			DIAM = 0.0030 cm  => SCALE = 5.65E+07
C			DIAM = 0.0016 cm  => SCALE = 3.02E+07
C
	DOUBLE PRECISION JNAGLU,JNAPI,JCLF,JCLB,JKCL,JNAB,JNACLB
C
	DO 700 II=1,2
	IF (II.EQ.1) SCALE=1.0
C SCALE FOR A 30 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=5.65D+07
C SCALE FOR A 25 MICRON TUBULE DIAMETER
	IF (II.EQ.2) SCALE=4.71D+07
C SCALE FOR A 16 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=3.02D+07
C SCALE FOR A 15 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=2.83D+07
C
	SOL(1)='  NA '
	SOL(2)='  K  '
	SOL(3)='  CL '
	SOL(4)=' HCO3'
	SOL(5)='H2CO3'
	SOL(6)=' CO2 '
	SOL(7)=' HPO4'
	SOL(8)='H2PO4'
	SOL(9)=' UREA'
	SOL(10)=' NH3 '
	SOL(11)=' NH4 '
	SOL(12)='  H  '
	SOL(13)=' HCO2'
	SOL(14)='H2CO2'
	SOL(15)=' GLUC'
C
	NCELLS=1
C
	FEVM(X,T)=FEVM(X,T)*SCALE
	FEVS(X,T)=FEVS(X,T)*SCALE
	DO 520 I=1,SOLS
	FEKM(I,X,T)=FEKM(I,X,T)*SCALE
  520   FEKS(I,X,T)=FEKS(I,X,T)*SCALE
	DO 530 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)*SCALE
	JV(K,X,T)=JV(K,X,T)*SCALE
	FIVS(K,X,T)=FIVS(K,X,T)*SCALE
	DO 530 I=1,SOLS
	FIKM(K,I,X,T)=FIKM(K,I,X,T)*SCALE
	JK(K,I,X,T)=JK(K,I,X,T)*SCALE
  530   FIKS(K,I,X,T)=FIKS(K,I,X,T)*SCALE
C
	FLUXV=FEVM(X,T)
	DO 550 I=1,SOLS
  550   FLUXS(I)=FEKM(I,X,T)
	DO 560 K=1,1
	FLUXV=FLUXV + FIVM(K,X,T)
	DO 555 I=1,SOLS
  555   FLUXS(I)=FLUXS(I) + FIKM(K,I,X,T)
  560	CONTINUE
C
C NET COTRANSPORTERS
C
	DO 110 K=1,1
	DO 110 I=1,SOLS
	CIKM(K,I)=0.
	CIKS(K,I)=0.
  110   CJK(K,I)=0.
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AMI(K)*(XM(J,X)-XI(K,J,X))*SCALE
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AIS(K)*(XI(K,J,X)-XS(J,X))*SCALE
  120   CJK(K,I)=CJK(K,I)+LIS(K,I,J)*AIE(K)*(XI(K,J,X)-XE(J,X))*SCALE
C
	CIKM(1,1)=CIKM(1,1) + JNHE3(1,1,X,T)*SCALE
	CIKM(1,11)=CIKM(1,11) + JNHE3(1,3,X,T)*SCALE
	CIKM(1,12)=CIKM(1,12) + JNHE3(1,2,X,T)*SCALE
C
C THIS PORTION OF THE CODE BACK-CALCULATES THE INDIVIDUAL COTRANSPORTERS OF THE CELL:
C	
	JNAGLU= CIKM(1,15)
	JNAPI=  CIKM(1,8)
	JCLF=  -CIKM(1,13)
	JCLB=  -CIKM(1,4)
	JKCL=   CIKS(1,2)+CJK(1,2)
	JNACLB=-(CIKS(1,3)+CJK(1,3) - JKCL)
	JNAB=  (CIKS(1,4)+CJK(1,4) - 2.0*JNACLB)/3.0
C
C PUMPS
C
	DO 130 K=1,1
	DO 130 I=1,SOLS
	PIKM(K,I)=0.
	PIKS(K,I)=0.
  130   PJK(K,I)=0.
C
	DO 140 K=1,1
	DO 140 I=1,SOLS
	PIKM(K,I)=PIKM(K,I)+AMI(K)*ATMI(K,I,X)*SCALE
	PIKS(K,I)=PIKS(K,I)+AIS(K)*ATIS(K,I,X)*SCALE
  140   PJK(K,I)=PJK(K,I)+AIE(K)*ATIS(K,I,X)*SCALE
C
	DO 150 K=1,1
	DO 150 I=1,SOLS
	GIKM(K,I)=FIKM(K,I,X,T)-PIKM(K,I)-CIKM(K,I)
	GIKS(K,I)=FIKS(K,I,X,T)-PIKS(K,I)-CIKS(K,I)
  150   GJK(K,I)=JK(K,I,X,T)-PJK(K,I)-CJK(K,I)
C
	IF (II.EQ.2) WRITE(168,500) DX*FLOAT(X-1),
     1    FEKM(1,X,T),FIKM(1,1,X,T),FLUXS(1),FEKM(1,X,T)/FLUXS(1),
     1    PIKS(1,1)+PJK(1,1),(PIKS(1,1)+PJK(1,1))/FLUXS(1)
  500   FORMAT (8D16.8)
C
	FEVM(X,T)=FEVM(X,T)/SCALE
	FEVS(X,T)=FEVS(X,T)/SCALE
	DO 620 I=1,SOLS
	FEKM(I,X,T)=FEKM(I,X,T)/SCALE
  620   FEKS(I,X,T)=FEKS(I,X,T)/SCALE
	DO 630 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)/SCALE
	JV(K,X,T)=JV(K,X,T)/SCALE
	FIVS(K,X,T)=FIVS(K,X,T)/SCALE
	DO 630 I=1,SOLS
	FIKM(K,I,X,T)=FIKM(K,I,X,T)/SCALE
	JK(K,I,X,T)=JK(K,I,X,T)/SCALE
  630   FIKS(K,I,X,T)=FIKS(K,I,X,T)/SCALE
C
  700	CONTINUE
	RETURN
	END
