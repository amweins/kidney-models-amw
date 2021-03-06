	SUBROUTINE CNTRESC
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
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12,81),ATIS(3,12,81),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
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
     1   LMI,LIS,ATMI,ATIS,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION SS(4),VAE1,
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
C
C NET COTRANSPORTERS
C
	DO 110 K=1,3
	DO 110 I=1,SOLS
	CIKM(K,I)=0.
	CIKS(K,I)=0.
  110   CJK(K,I)=0.
C
	DO 120 K=1,3
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AI0(K)*(XM(J,X)-XI(K,J,X))
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AI0(K)*(XI(K,J,X)-XS(J,X))
  120   CJK(K,I)=CJK(K,I)+LIS(K,I,J)*AIE(K)*(XI(K,J,X)-XE(J,X))
C
C AE1 COTRANSPORTERS
C
 	SS(1)=CI(2,3,X,T)
 	SS(2)=CI(2,4,X,T)
 	SS(3)=CE(3,X,T)
 	SS(4)=CE(4,X,T)
 	CALL AE1(VAE1,SS)
 	CJK(2,3)=CJK(2,3)-NAE1(2)*AIE(2)*VAE1
 	CJK(2,4)=CJK(2,4)+NAE1(2)*AIE(2)*VAE1
 	SS(3)=CS(3,X,T)
 	SS(4)=CS(4,X,T)
 	CALL AE1(VAE1,SS)
 	CIKS(2,3)=CIKS(2,3)-NAE1(2)*AI0(2)*VAE1
 	CIKS(2,4)=CIKS(2,4)+NAE1(2)*AI0(2)*VAE1
C
	SS(1)=CM(3,X,T)
	SS(2)=CM(4,X,T)
	SS(3)=CI(3,3,X,T)
	SS(4)=CI(3,4,X,T)
	CALL AE1(VAE1,SS)
	CIKM(3,3)=CIKM(3,3)-NAE1(3)*AI0(3)*VAE1
	CIKM(3,4)=CIKM(3,4)+NAE1(3)*AI0(3)*VAE1
C
C
C PUMPS
C
	DO 130 K=1,3
	DO 130 I=1,SOLS
	PIKM(K,I)=0.
	PIKS(K,I)=0.
  130   PJK(K,I)=0.
C
	DO 140 K=1,3
	DO 140 I=1,SOLS
	PIKM(K,I)=PIKM(K,I)+AI0(K)*ATMI(K,I,X)
	PIKS(K,I)=PIKS(K,I)+AI0(K)*ATIS(K,I,X)
  140   PJK(K,I)=PJK(K,I)+AIE(K)*ATIS(K,I,X)
C
C
	WRITE (21,480) (SOL(I),
     1  (CIKM(K,I),CJK(K,I),CIKS(K,I),K=1,3),I=1,SOLS)
	WRITE (21,490) (SOL(I),
     1  (PIKM(K,I),PJK(K,I),PIKS(K,I),K=1,3),I=1,SOLS)
C
  480   FORMAT(1H1,//,14X,'CIKMP',6X,'CJKP',9X,'CIKSP',7X,
     1  'CIKMA',6X,'CJKA',9X,'CIKSA',7X,
     1  'CIKMB',6X,'CJKB',9X,'CIKSB',7X,/,(4X,A5,9D12.4))
  490   FORMAT(//,14X,'PIKMP',6X,'PJKP',9X,'PIKSP',7X,
     1  'PIKMA',6X,'PJKA',9X,'PIKSA',7X,
     1  'PIKMB',6X,'PJKB',9X,'PIKSB',7X,/,(4X,A5,9D12.4))
C
	RETURN
	END
