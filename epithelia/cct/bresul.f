	SUBROUTINE BRESUL
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
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION QE(12),QI(3,12),QEV,QIV(3),
     1  OSME,OSMI(3),OSMM,OSMS,FLUXV,FLUXS(12),HMI11
C
C
C       QE-     DERIVATIVE OF SPECIES CONTENT OF THE CHANNEL (MMOL/S)
C       QI-     DERIVATIVE OF SPECIES CONTENT OF THE CELL (MMOL/S)
C       QEV-    DERIVATIVE VOLUME CONTENT OF THE CHANNEL (ML/S)
C       QIV-    DERIVATIVE VOLUME CONTENT OF THE  CELL (ML/S)
C       OSM.-   OSMOLALITY OF INDICATED COMPARTMENT
C	       .= E,I,M,S
C	FLUXV-	TOTAL EPITHELIAL VOLUME FLUX
C	FLUXS-	TOTAL EPITHELIAL SOLUTE FLUX
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
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
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
C DERIVED QUANTITIES ARE COMPUTED:OSMOLALITY,CURRENT FLOW,AND THE
C NET QUANTITY "Q(I)" OF SOLUTE INTRODUCED INTO THE COMPARTMENT
C
  300   DO 310 I=1,SOLS
	QE(I)=FEKM(I,T)-FEKS(I,T)
	DO 310 K=1,3
	QE(I)=QE(I)+JK(K,I,T)
  310   QI(K,I)=FIKM(K,I,T)-FIKS(K,I,T)-JK(K,I,T)
	QEV=FEVM(T)-FEVS(T)
	DO 312 K=1,3
	QEV=QEV+JV(K,T)
  312   QIV(K)=FIVM(K,T)-FIVS(K,T)-JV(K,T)
C
	LCHM=-1.*DLOG10(CM(SOLS))
	LCHS=-1.*DLOG10(CS(SOLS))
	LCHE=-1.*DLOG10(CE(SOLS,T))
	DO 314 K=1,3
  314   LCHI(K)=-1.*DLOG10(CI(K,SOLS,T))
C
	OSME=0.D0
	OSMM=IMPM
	OSMS=IMPS
	DO 318 K=1,3
  318   OSMI(K)=IMP(K)
	DO 320 I=1,SOLS
	OSME=OSME+CE(I,T)
	OSMM=OSMM+CM(I)
	OSMS=OSMS+CS(I)
	DO 320 K=1,3
  320   OSMI(K)=OSMI(K)+CI(K,I,T)
C
	CURE=0.D0
	DO 340 K=1,3
  340   CURI(K)=0.D0
	DO 345 I=1,SOLS
	CURE=CURE+F*Z(I)*FEKM(I,T)
	DO 345 K=1,3
  345   CURI(K)=CURI(K)+F*Z(I)*FIKM(K,I,T)
C
C
	WRITE (21,396) TIME
  396   FORMAT (1H1,/,5X,'TIME=',F9.4)
C
	WRITE (21,398) L(T),CHVL(T),(CLVL(K,T),K=1,3),AE(T)
	WRITE (21,400) VM,VE(T),(VI(K,T),K=1,3),VS
	WRITE (21,402) PM,PE(T),(PI(K,T),K=1,3),PS
C
	WRITE (21,410) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,3),CS(I),I=1,4)	
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,3),CS(I),I=5,5)	
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,3),CS(I),I=6,SOLS-3)	
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,3),CS(I),I=SOLS-2,SOLS-2)	
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,3),CS(I),I=SOLS-1,SOLS-1)	
	WRITE (21,412) SOL(SOLS),
     1  CM(SOLS),CE(SOLS,T),(CI(K,SOLS,T),K=1,3),CS(SOLS)
	WRITE (21,420) IMPM,(IMP(K),K=1,3),IMPS
	WRITE (21,430) OSMM,OSME,(OSMI(K),K=1,3),OSMS
	WRITE (21,440) LCHM,LCHE,(LCHI(K),K=1,3),LCHS
	WRITE (21,442) (CBUF(K,T),K=1,3),(HCBUF(K,T),K=1,3)
C
	WRITE (21,460) (SOL(I),XM(I)-XS(I),XE(I)-XS(I),
     1      (XI(K,I)-XS(I),K=1,3),I=1,SOLS)
	WRITE (21,470) (FIVM(K,T),JV(K,T),FIVS(K,T),K=1,3)
	WRITE (21,480) (SOL(I),
     1  (FIKM(K,I,T),JK(K,I,T),FIKS(K,I,T),K=1,3),I=1,SOLS)
	WRITE (21,485) (CURI(K),K=1,3)
	WRITE (21,490) FEVM(T),FEVS(T),QEV,(QIV(K),K=1,3)
	WRITE (21,500) (SOL(I),FEKM(I,T),FEKS(I,T),QE(I),
     1  (QI(K,I),K=1,3),I=1,SOLS)
	WRITE (21,505) CURE
C
C
  398   FORMAT(//,14X,'L',8X,'CHVL',7X,'CLVLP',7X,'CLVLA',
     1  7X,'CLVLB',9X,'AE',/,6X,6D12.4)
  400   FORMAT(//,16X,'VM',10X,'VE',9X,'VIP',9X,'VIA',9X,
     1  'VIB',10X,'VS',/,8X,6F12.3)
  402   FORMAT(//,16X,'PM',10X,'PE',9X,'PIP',9X,'PIA',9X,
     1  'PIB',10X,'PS',/,8X,6F12.3)
C
  410   FORMAT(//,16X,'CM',10X,'CE',9X,'CIP',9X,'CIA',9X,
     1  'CIB',10X,'CS',/,(4X,A5,6F12.6))
  412   FORMAT(4X,A5,1X,6D12.4)
  414   FORMAT(4X,A5,6F12.6)
  420   FORMAT(5X,'IMP',1X,F12.6,12X,4F12.6)
  430   FORMAT(5X,'OSM',1X,6F12.6)
  440   FORMAT(6X,'PH',1X,6F12.6)
  442   FORMAT(4X,'CBUF',1X,24X,3F12.6,/,4X,'HCBUF',24X,3F12.6)
C
  460   FORMAT(//,15X,'XM-XS',7X,'XE-XS',6X,'XIP-XS',
     1  6X,'XIA-XS',6X,'XIB-XS',/,(4X,A5,5F12.4))
  470   FORMAT(1H1,//,14X,'FIVMP',7X,'JVP',9X,'FIVSP',7X,
     1  'FIVMA',7X,'JVA',9X,'FIVSA',7X,
     1  'FIVMB',7X,'JVB',9X,'FIVSB',7X,/,9X,9D12.4)
  480   FORMAT(//,14X,'FIKMP',7X,'JKP',9X,'FIKSP',7X,
     1  'FIKMA',7X,'JKA',9X,'FIKSA',7X,
     1  'FIKMB',7X,'JKB',9X,'FIKSB',7X,/,(4X,A5,9D12.4))
  485   FORMAT(4X,'CUR',2X,D12.4,24X,D12.4,24X,D12.4)
  490   FORMAT(//,15X,'FEVM',8X,'FEVS',20X,'QEV',9X,
     1  'QIVP',8X,'QIVA',8X,'QIVB',/,9X,2D12.4,12X,4D12.4)
  500   FORMAT(//,15X,'FEKM',8X,'FEKS',20X,'QE',10X,
     1  'QIP',9X,'QIA',9X,'QIB',/,(4X,A5,2D12.4,12X,4D12.4))
  505   FORMAT(4X,'CUR',2X,D12.4)
C
C OUTPUT FOR AN AUXILLIARY FILE--E.G. PERMEABILITY CALCULATIONS
C
	FLUXV=FEVM(T)
	DO 550 I=1,SOLS
  550   FLUXS(I)=FEKM(I,T)
	DO 560 K=1,3
	FLUXV=FLUXV + FIVM(K,T)
	DO 555 I=1,SOLS
  555   FLUXS(I)=FLUXS(I) + FIKM(K,I,T)
  560	CONTINUE
C
C These are for the permeability calculations.
C	WRITE (23,570) FLUXV,(FLUXS(I),I=1,4),
C     1   FLUXS(6),FLUXS(9),FLUXS(11),FLUXS(10)
C  570   FORMAT (8E16.8,/,8E16.8)
C
C These are for the luminal membrane display
C	HMI11=HMI(1,1)*((1.+CM(1)/.03)**-1)*
C     1      (1.-CI(1,1,T)/.05)*1.666667
C	WRITE (23,570) CM(1),VM,VM-VI(1,T),HMI11,
C     1   (XM(I)-XI(1,I),I=1,2),
C     1   (1.D6*FIKM(1,I,T),I=1,2),
C     1   (1.D6*FLUXS(I),I=1,2)
C  570   FORMAT (8E16.8)
C
C These are as needed.
C        WRITE (23,570) VM
C     1  (1.D6*FEKM(I,T),I=1,3),1.D6*(FEKM(4,T)-FEKM(12,T)),
C     1  (1.D6*FLUXS(I),I=1,3),1.D6*(FLUXS(4)-FLUXS(12))
  570   FORMAT (8E16.8)
C
	RETURN
	END
