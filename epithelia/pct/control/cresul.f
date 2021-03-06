	SUBROUTINE CRESUL
C
	INTEGER SOLS,T,TAU,TLIM
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,L0,L(2),
     1   PKC,PKF,PKN,PKP,KHY(4),KDHY(4)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(15),IMPM,LCHM,XM(15),
     1   VS,PS,CS(15),IMPS,LCHS,XS(15)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(2),PE(2),CE(15,2),LCHE,XE(15),
     1   FEVM(2),FEKM(15,2),FEVS(2),FEKS(15,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),ATMI(3,15),ATIS(3,15),
     1   VI(3,2),PI(3,2),CI(3,15,2),IMP(3),LCHI(3),XI(3,15),
     1   FIVM(3,2),FIKM(3,15,2),FIVS(3,2),FIKS(3,15,2),CURI(3),
     1   JV(3,2),JK(3,15,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3),JHK(3),JHP(3),JAE1(3),JTSC(3),JNHE3(3,3),
     1   QIAMM
C COTRANSPORTER INFORMATION
	INTEGER  
     1   NMI,ORDMI(9),IMI(9,9),JMI(9,9),MMI(9,9),
     1   NIS,ORDIS(9),IIS(9,9),JIS(9,9),MIS(9,9)
	DOUBLE PRECISION RMI(9),RIS(9)
C CONTROL VARIABLES
	DOUBLE PRECISION 
     1   HCON(10),HCON0(10),DHCON(10,4),CON(4),CON0(10,4)
C
	CHARACTER*5 SOL(15)
C
	DOUBLE PRECISION 
     1    PIKM(3,15),PJK(3,15),PIKS(3,15),
     1    CIKM(3,15),CJK(3,15),CIKS(3,15)
C
C	PIKM-	ACTIVE TRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	PJK-	ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	PIKS-	ACTIVE TRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	CIKM-	COTRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	CJK-	COTRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	CIKS-	COTRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C
        COMMON SOLS,T,TAU,TLIM,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AMI,AIS,CLVL0,IMP0,CLVL,
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
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM
	COMMON/COTR/
     1   NMI,ORDMI,IMI,JMI,MMI,
     1   NIS,ORDIS,IIS,JIS,MIS,RMI,RIS,
     1   HCON,HCON0,DHCON,CON,CON0
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
	SOL(13)=' HCO2'
	SOL(14)='H2CO2'
	SOL(15)=' GLUC'
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
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AMI(K)*(XM(J)-XI(K,J))
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AIS(K)*(XI(K,J)-XS(J))
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
	PIKM(K,I)=PIKM(K,I)+AMI(K)*ATMI(K,I)
	PIKS(K,I)=PIKS(K,I)+AIS(K)*ATIS(K,I)
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
