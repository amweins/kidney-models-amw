	SUBROUTINE DCTRESB
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
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2)
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
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3
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
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION QE(12),QI(3,12),QEV,QIV(3),
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
	WRITE (21,396) TIME,DIST,AM(X,T),SM(X,T)
  396   FORMAT (1H1,/,5X,'TIME=',F15.4,7X,'DIST=',F8.6,
     1  7X,'AM=',D10.4,7X,'SM=',D10.4)
C
C
	WRITE (21,398) L(X,T),CHVL(X,T),(CLVL(K,X,T),K=1,1),AE(X,T)
	WRITE (21,400) VM(X,T),VE(X,T),(VI(K,X,T),K=1,1),VS(X,T)
	WRITE (21,402) PM(X,T),PE(X,T),(PI(K,X,T),K=1,1),PS(X,T)
C
	WRITE (21,410) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=1,4)	
	WRITE (21,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=5,5)	
	WRITE (21,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),CS(I,X,T),I=6,SOLS-3)	
	WRITE (21,412) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),
     1   CS(I,X,T),I=SOLS-2,SOLS-2)	
	WRITE (21,414) (SOL(I),
     1  CM(I,X,T),CE(I,X,T),(CI(K,I,X,T),K=1,1),
     1   CS(I,X,T),I=SOLS-1,SOLS-1)	
	WRITE (21,412) SOL(SOLS),
     1  CM(SOLS,X,T),CE(SOLS,X,T),(CI(K,SOLS,X,T),K=1,1),CS(SOLS,X,T)
	WRITE (21,420) IMPM(X,T),(IMP(K,X),K=1,1),IMPS(X,T)
	WRITE (21,430) OSMM,OSME,(OSMI(K),K=1,1),OSMS
	WRITE (21,440) LCHM(X),LCHE(X),(LCHI(K,X),K=1,1),LCHS(X)
	WRITE (21,442) (CBUF(K,X,T),K=1,1),(HCBUF(K,X,T),K=1,1)
C
	WRITE (21,460) (SOL(I),XM(I,X)-XS(I,X),XE(I,X)-XS(I,X),
     1      (XI(K,I,X)-XS(I,X),K=1,1),I=1,SOLS)
	WRITE (21,470) (FIVM(K,X,T),JV(K,X,T),FIVS(K,X,T),K=1,1)
	WRITE (21,480) (SOL(I),
     1  (FIKM(K,I,X,T),JK(K,I,X,T),FIKS(K,I,X,T),K=1,1),I=1,SOLS)
	WRITE (21,485) (CURI(K,X),K=1,1)
	WRITE (21,490) FEVM(X,T),FEVS(X,T),QEV,(QIV(K),K=1,1)
	WRITE (21,500) (SOL(I),FEKM(I,X,T),FEKS(I,X,T),QE(I),
     1  (QI(K,I),K=1,1),I=1,SOLS)
	WRITE (21,505) CURE(X)
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
C     1  6X,'XIA-XS',6X,'XIB-XS',/,(4X,A5,3F12.4,))
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
