	SUBROUTINE AHLRESC
C
	INTEGER SOLS,X,T
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY(5),KDHY(5),L0,L(81,2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(81,2),PM(81,2),CM(12,81,2),
     1   IMPM(81,2),LCHM(81),XM(12,81),
     1   VS(81,2),PS(81,2),CS(12,81,2),
     1   IMPS(81,2),LCHS(81),XS(12,81),
     1   TL,DX,RM0,MUM,ETA,
     1   SM(81,2),AM(81,2),FVM(81,2),FKM(13,81,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(81,2),MUA,CHVL0,CHVL(81,2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(81,2),PE(81,2),CE(12,81,2),LCHE(81),XE(12,81),
     1   FEVM(81,2),FEKM(12,81,2),FEVS(81,2),FEKS(12,81,2),CURE(81)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,81,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,81,2),HCBUF(3,81,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),
     1   ATMI(3,12,81),ATIS(3,12,81),ATIE(3,12,81),
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2),
     1   NNKCC(3),NKCL(3),JNKCC(3,4,81,2),JKCC(3,3,81,2)
C
        COMMON SOLS,T,X,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
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
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION 
     1    PIKM(3,12),PJK(3,12),PIKS(3,12),
     1    CIKM(3,12),CJK(3,12),CIKS(3,12)
C
C	PIKM-	ACTIVE TRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	PJK-	ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	PIKS-	ACTIVE TRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	CIKM-	COTRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	CJK-	COTRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	CIKS-	COTRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C
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
C
	NCELLS=1
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
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AI0(K)*(XM(J,X)-XI(K,J,X))
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AI0(K)*(XI(K,J,X)-XS(J,X))
  120   CJK(K,I)=CJK(K,I)+LIS(K,I,J)*AIE(K)*(XI(K,J,X)-XE(J,X))
C
C ADD IN THE INDIVIDUAL COTRANSPORTERS OF THE CELL
C
C  NKCC COTRANSPORTER
	CIKM(1,1) = CIKM(1,1) + JNKCC(1,1,X,T)
	CIKM(1,2) = CIKM(1,2) + JNKCC(1,2,X,T)
	CIKM(1,3) = CIKM(1,3) + JNKCC(1,3,X,T)
	CIKM(1,11) = CIKM(1,11) + JNKCC(1,4,X,T)
C
C  KCC COTRANSPORTER
	CJK(1,2) = CJK(1,2) + JKCC(1,1,X,T)*(AIE(1)/(AIE(1)+AI0(1)))
	CJK(1,3) = CJK(1,3) + JKCC(1,2,X,T)*(AIE(1)/(AIE(1)+AI0(1)))
	CJK(1,11) = CJK(1,11) + JKCC(1,3,X,T)*(AIE(1)/(AIE(1)+AI0(1)))
C
	CIKS(1,2) = CIKS(1,2) + JKCC(1,1,X,T)*(AI0(1)/(AIE(1)+AI0(1)))
	CIKS(1,3) = CIKS(1,3) + JKCC(1,2,X,T)*(AI0(1)/(AIE(1)+AI0(1)))
	CIKS(1,11) = CIKS(1,11) + JKCC(1,3,X,T)*(AI0(1)/(AIE(1)+AI0(1)))
C
C  TSC COTRANSPORTER
	CIKM(1,1)=CIKM(1,1) + JTSC(1,X,T)
	CIKM(1,3)=CIKM(1,3) + JTSC(1,X,T)
C
C  NHE3 COTRANSPORTER
	CIKM(1,1)=CIKM(1,1) + JNHE3(1,1,X,T)
	CIKM(1,11)=CIKM(1,11) + JNHE3(1,3,X,T)
	CIKM(1,12)=CIKM(1,12) + JNHE3(1,2,X,T)
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
	PIKM(K,I)=PIKM(K,I)+AI0(K)*ATMI(K,I,X)
	PIKS(K,I)=PIKS(K,I)+AI0(K)*ATIS(K,I,X)
  140   PJK(K,I)=PJK(K,I)+AIE(K)*ATIE(K,I,X)
C
C
	WRITE (21,480) (SOL(I),
     1  (CIKM(K,I),CJK(K,I),CIKS(K,I),K=1,1),I=1,SOLS)
	WRITE (21,490) (SOL(I),
     1  (PIKM(K,I),PJK(K,I),PIKS(K,I),K=1,1),I=1,SOLS)
C
  480   FORMAT(1H1,//,14X,'CIKMP',6X,'CJKP',9X,'CIKSP',7X,
     1  /,(4X,A5,3D12.4))
  490   FORMAT(//,14X,'PIKMP',6X,'PJKP',9X,'PIKSP',7X,
     1  /,(4X,A5,3D12.4))
C  480   FORMAT(1H1,//,14X,'CIKMP',6X,'CJKP',9X,'CIKSP',7X,
C     1  'CIKMA',6X,'CJKA',9X,'CIKSA',7X,
C     1  'CIKMB',6X,'CJKB',9X,'CIKSB',7X,/,(4X,A5,3D12.4))
C  490   FORMAT(//,14X,'PIKMP',6X,'PJKP',9X,'PIKSP',7X,
C     1  'PIKMA',6X,'PJKA',9X,'PIKSA',7X,
C     1  'PIKMB',6X,'PJKB',9X,'PIKSB',7X,/,(4X,A5,3D12.4))
C
	RETURN
	END
