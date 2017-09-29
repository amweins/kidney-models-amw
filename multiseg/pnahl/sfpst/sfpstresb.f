	SUBROUTINE SFPSTRESB
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
	CHARACTER*5 SOL(15)
C
	DOUBLE PRECISION QE(15),QI(3,15),QEV,QIV(3),
     1  OSME,OSMI(3),OSMM,OSMS
C
C       QE-     DERIVATIVE OF SPECIES CONTENT OF THE CHANNEL (MMOL/S)
C       QI-     DERIVATIVE OF SPECIES CONTENT OF THE CELL (MMOL/S)
C       QEV-    DERIVATIVE VOLUME CONTENT OF THE CHANNEL (ML/S)
C       QIV-    DERIVATIVE VOLUME CONTENT OF THE  CELL (ML/S)
C       OSM.-   OSMOLALITY OF INDICATED COMPARTMENT
C	       .= E,I,M,S
C
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
	NCELLS=1
C
C DERIVED QUANTITIES ARE COMPUTED:OSMOLALITY,CURRENT FLOW,AND THE
C NET QUANTITY "Q(I)" OF SOLUTE INTRODUCED INTO THE COMPARTMENT
C
  300   DO 310 I=1,SOLS
	QE(I)=FEKM(I,X,T)-FEKS(I,X,T)
	DO 310 K=1,1
	QE(I)=QE(I)+JK(K,I,X,T)
  310   QI(K,I)=FIKM(K,I,X,T)-FIKS(K,I,X,T)-JK(K,I,X,T)
	QEV=FEVM(X,T)-FEVS(X,T)
	DO 312 K=1,1
	QEV=QEV+JV(K,X,T)
  312   QIV(K)=FIVM(K,X,T)-FIVS(K,X,T)-JV(K,X,T)
C
	OSME=0.D0
	OSMM=IMPM(X,T)
	OSMS=IMPS(X,T)
	DO 318 K=1,1
  318   OSMI(K)=IMP(K,X)
	DO 320 I=1,SOLS
	OSME=OSME+CE(I,X,T)
	OSMM=OSMM+CM(I,X,T)
	OSMS=OSMS+CS(I,X,T)
	DO 320 K=1,1
  320   OSMI(K)=OSMI(K)+CI(K,I,X,T)
C
C
	WRITE (151,396) TIME,DIST,AM(X,T),SM(X,T)
  396   FORMAT (1H1,/,5X,'TIME=',F15.4,7X,'DIST=',F8.6,
     1  7X,'AM=',D10.4,7X,'SM=',D10.4)
C
C
	WRITE (151,398) L(X,T),CHVL(X,T),(CLVL(K,X,T),K=1,1),AE(X,T)
	WRITE (151,400) VM(X,T),VE(X,T),(VI(K,X,T),K=1,1),VS(X,T)
	WRITE (151,402) PM(X,T),PE(X,T),(PI(K,X,T),K=1,1),PS(X,T)
C
	WRITE (151,410) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=1,4)	
	WRITE (151,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=5,5)	
	WRITE (151,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=6,9)	
	WRITE (151,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=10,10)	
	WRITE (151,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=11,11)
	WRITE (151,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=12,12)
	WRITE (151,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=13,13)
	WRITE (151,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=14,14)
	WRITE (151,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=15,15)
	WRITE (151,420) IMPM(X,T),(IMP(K,X),K=1,1),IMPS(X,T)
	WRITE (151,430) OSMM,OSME,(OSMI(K),K=1,1),OSMS
	WRITE (151,440) LCHM(X),LCHE(X),(LCHI(K,X),K=1,1),LCHS(X)
	WRITE (151,442) (CBUF(K,X,T),K=1,1),(HCBUF(K,X,T),K=1,1)
C
	WRITE (151,460) (SOL(I),XM(I,X)-XS(I,X),XE(I,X)-XS(I,X),
     1      (XI(K,I,X)-XS(I,X),K=1,1),I=1,SOLS)
	WRITE (151,470) (FIVM(K,X,T),JV(K,X,T),FIVS(K,X,T),K=1,1)
	WRITE (151,480) (SOL(I),
     1  (FIKM(K,I,X,T),JK(K,I,X,T),FIKS(K,I,X,T),K=1,1),I=1,SOLS)
	WRITE (151,485) (CURI(K,X),K=1,1)
	WRITE (151,490) FEVM(X,T),FEVS(X,T),QEV,(QIV(K),K=1,1)
	WRITE (151,500) (SOL(I),FEKM(I,X,T),FEKS(I,X,T),QE(I),
     1  (QI(K,I),K=1,1),I=1,SOLS)
	WRITE (151,505) CURE(X)
C
C
  398   FORMAT(//,14X,'L',8X,'CHVL',7X,'CLVLP',
     1  9X,'AE',/,6X,6D12.4)
  400   FORMAT(//,16X,'VM',10X,'VE',9X,'VIP',
     1  10X,'VS',/,8X,6F12.3)
  402   FORMAT(//,16X,'PM',10X,'PE',9X,'PIP',
     1  10X,'PS',/,8X,6F12.3)
C
  410   FORMAT(//,16X,'CM',10X,'CE',9X,'CIP',
     1  10X,'CS',/,(4X,A5,4F12.6))
  412   FORMAT(4X,A5,1X,4D12.4)
  414   FORMAT(4X,A5,4F12.6)
  420   FORMAT(5X,'IMP',1X,F12.6,12X,4F12.6)
  430   FORMAT(5X,'OSM',1X,6F12.6)
  440   FORMAT(6X,'PH',1X,6F12.6)
  442   FORMAT(4X,'CBUF',1X,24X,F12.6,
     1    /,4X,'HCBUF',24X,F12.6)
C
  460   FORMAT(//,15X,'XM-XS',7X,'XE-XS',6X,'XIP-XS',
     1  /,(4X,A5,3F12.4))
  470   FORMAT(1H1,//,14X,'FIVMP',7X,'JVP',9X,'FIVSP',7X,
     1  /,9X,9D12.4)
  480   FORMAT(//,14X,'FIKMP',7X,'JKP',9X,'FIKSP',7X,
     1  /,(4X,A5,3D12.4))
  485   FORMAT(4X,'CUR',2X,D12.4,24X,D12.4,24X,D12.4)
  490   FORMAT(//,15X,'FEVM',8X,'FEVS',20X,'QEV',9X,
     1  'QIVP',/,9X,2D12.4,12X,4D12.4)
  500   FORMAT(//,15X,'FEKM',8X,'FEKS',20X,'QE',10X,
     1  'QIP',/,(4X,A5,2D12.4,12X,2D12.4))
  505   FORMAT(4X,'CUR',2X,D12.4)
C
C
C  398   FORMAT(//,14X,'L',8X,'CHVL',7X,'CLVLP',7X,'CLVLA',
C     1  7X,'CLVLB',9X,'AE',/,6X,6D12.4)
C  400   FORMAT(//,16X,'VM',10X,'VE',9X,'VIP',9X,'VIA',9X,
C     1  'VIB',10X,'VS',/,8X,6F12.3)
C  402   FORMAT(//,16X,'PM',10X,'PE',9X,'PIP',9X,'PIA',9X,
C     1  'PIB',10X,'PS',/,8X,6F12.3)
C
C  410   FORMAT(//,16X,'CM',10X,'CE',9X,'CIP',9X,'CIA',9X,
C     1  'CIB',10X,'CS',/,(4X,A5,4F12.6))
C  412   FORMAT(4X,A5,1X,4D12.4)
C  414   FORMAT(4X,A5,4F12.6)
C  420   FORMAT(5X,'IMP',1X,F12.6,12X,4F12.6)
C  430   FORMAT(5X,'OSM',1X,6F12.6)
C  440   FORMAT(6X,'PH',1X,6F12.6)
C  442   FORMAT(4X,'CBUF',1X,24X,F12.6,
C     1    /,4X,'HCBUF',24X,F12.6)
C
C  460   FORMAT(//,15X,'XM-XS',7X,'XE-XS',6X,'XIP-XS',
C     1  6X,'XIA-XS',6X,'XIB-XS',/,(4X,A5,3F12.4))
C  470   FORMAT(1H1,//,14X,'FIVMP',7X,'JVP',9X,'FIVSP',7X,
C     1  'FIVMA',7X,'JVA',9X,'FIVSA',7X,
C     1  'FIVMB',7X,'JVB',9X,'FIVSB',7X,/,9X,9D12.4)
C  480   FORMAT(//,14X,'FIKMP',7X,'JKP',9X,'FIKSP',7X,
C     1  'FIKMA',7X,'JKA',9X,'FIKSA',7X,
C     1  'FIKMB',7X,'JKB',9X,'FIKSB',7X,/,(4X,A5,3D12.4))
C  485   FORMAT(4X,'CUR',2X,D12.4,24X,D12.4,24X,D12.4)
C  490   FORMAT(//,15X,'FEVM',8X,'FEVS',20X,'QEV',9X,
C     1  'QIVP',8X,'QIVA',8X,'QIVB',/,9X,2D12.4,12X,4D12.4)
C  500   FORMAT(//,15X,'FEKM',8X,'FEKS',20X,'QE',10X,
C     1  'QIP',9X,'QIA',9X,'QIB',/,(4X,A5,2D12.4,12X,2D12.4))
C  505   FORMAT(4X,'CUR',2X,D12.4)
C
	RETURN
	END
