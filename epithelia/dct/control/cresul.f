	SUBROUTINE CRESUL
C
	INTEGER SOLS,T
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY(4),KDHY(4),L0,L(2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(12),IMPM,LCHM,XM(12),
     1   VS,PS,CS(12),IMPS,LCHS,XS(12)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(2),PE(2),CE(12,2),LCHE,XE(12),
     1   FEVM(2),FEKM(12,2),FEVS(2),FEKS(12,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3),JHK(3),JHP(3),JAE1(3),JTSC(3),JNHE3(3,3)
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
        COMMON SOLS,T,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
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
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3
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
	DO 110 K=1,1
	DO 110 I=1,SOLS
	CIKM(K,I)=0.
	CIKS(K,I)=0.
  110   CJK(K,I)=0.
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AI0(K)*(XM(J)-XI(K,J))
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AI0(K)*(XI(K,J)-XS(J))
  120   CJK(K,I)=CJK(K,I)+LIS(K,I,J)*AIE(K)*(XI(K,J)-XE(J))
C
C
	CIKM(1,1)=CIKM(1,1) + JTSC(1)
	CIKM(1,3)=CIKM(1,3) + JTSC(1)
	CIKM(1,1)=CIKM(1,1) + JNHE3(1,1)
	CIKM(1,11)=CIKM(1,11) + JNHE3(1,3)
	CIKM(1,12)=CIKM(1,12) + JNHE3(1,2)
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
	PIKM(K,I)=PIKM(K,I)+AI0(K)*ATMI(K,I)
	PIKS(K,I)=PIKS(K,I)+AI0(K)*ATIS(K,I)
  140   PJK(K,I)=PJK(K,I)+AIE(K)*ATIS(K,I)
C
C
	WRITE (21,480) (SOL(I),
     1  (CIKM(K,I),CJK(K,I),CIKS(K,I),K=1,1),I=1,SOLS)
	WRITE (21,490) (SOL(I),
     1  (PIKM(K,I),PJK(K,I),PIKS(K,I),K=1,1),I=1,SOLS)
C
  480   FORMAT(1H1,//,14X,
     1  'CIKM ',6X,'CJK ',9X,'CIKS ',7X,/,(4X,A5,3D12.4))
  490   FORMAT(//,14X,
     1  'PIKM ',6X,'PJK ',9X,'PIKS ',7X,/,(4X,A5,3D12.4))
C
	RETURN
	END