	SUBROUTINE BRESUL
C
	INTEGER SOLS,T
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
C
	CHARACTER*5 SOL(15)
C
	DOUBLE PRECISION QE(15),QI(3,15),QEV,QIV(3),
     1  OSME,OSMI(3),OSMM,OSMS,FLUXV,FLUXS(15),
     1  CURME,CURES,CURMI,CURIE,CURIS
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
C DERIVED QUANTITIES ARE COMPUTED:OSMOLALITY,CURRENT FLOW,AND THE
C NET QUANTITY "Q(I)" OF SOLUTE INTRODUCED INTO THE COMPARTMENT
C
  300   DO 310 I=1,SOLS
	QE(I)=FEKM(I,T)-FEKS(I,T)
	DO 310 K=1,1
	QE(I)=QE(I)+JK(K,I,T)
  310   QI(K,I)=FIKM(K,I,T)-FIKS(K,I,T)-JK(K,I,T)
	QEV=FEVM(T)-FEVS(T)
	DO 312 K=1,1
	QEV=QEV+JV(K,T)
  312   QIV(K)=FIVM(K,T)-FIVS(K,T)-JV(K,T)
C
	OSME=0.D0
	OSMM=IMPM
	OSMS=IMPS
	DO 318 K=1,1
  318   OSMI(K)=IMP(K)
	DO 320 I=1,SOLS
	OSME=OSME+CE(I,T)
	OSMM=OSMM+CM(I)
	OSMS=OSMS+CS(I)
	DO 320 K=1,1
  320   OSMI(K)=OSMI(K)+CI(K,I,T)
C
        CURME=0.D0
        CURES=0.D0
        CURMI=0.D0
        CURIE=0.D0
        CURIS=0.D0
	CURE=0.D0
	DO 340 K=1,1
  340   CURI(K)=0.D0
	DO 345 I=1,SOLS
	CURME=CURME+F*Z(I)*FEKM(I,T)
	CURES=CURES+F*Z(I)*FEKS(I,T)
	CURE=CURE+F*Z(I)*FEKM(I,T)
	DO 345 K=1,1
	CURMI=CURMI+F*Z(I)*FIKM(K,I,T)
	CURIE=CURIE+F*Z(I)*JK(K,I,T)
	CURIS=CURIS+F*Z(I)*FIKS(K,I,T)
  345   CURI(K)=CURI(K)+F*Z(I)*FIKM(K,I,T)
C
C
	WRITE (21,396) TIME
  396   FORMAT (1H1,/,5X,'TIME=',F9.4)
C
	WRITE (21,398) L(T),CHVL(T),(CLVL(K,T),K=1,1),AE(T)
	WRITE (21,400) VM,VE(T),(VI(K,T),K=1,1),VS
	WRITE (21,402) PM,PE(T),(PI(K,T),K=1,1),PS
C
	WRITE (21,410) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=1,4)	
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=5,5)	
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=6,9)	
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=10,10)	
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=11,11)
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=12,12)
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=13,13)
	WRITE (21,412) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=14,14)
	WRITE (21,414) (SOL(I),
     1  CM(I),CE(I,T),(CI(K,I,T),K=1,1),CS(I),I=15,15)
	WRITE (21,420) IMPM,(IMP(K),K=1,1),IMPS
	WRITE (21,430) OSMM,OSME,(OSMI(K),K=1,1),OSMS
	WRITE (21,440) LCHM,LCHE,(LCHI(K),K=1,1),LCHS
	WRITE (21,442) (CBUF(K,T),K=1,1),(HCBUF(K,T),K=1,1)
C
	WRITE (21,460) (SOL(I),XM(I)-XS(I),XE(I)-XS(I),
     1      (XI(K,I)-XS(I),K=1,1),I=1,SOLS)
	WRITE (21,470) (FIVM(K,T),JV(K,T),FIVS(K,T),K=1,1)
	WRITE (21,480) (SOL(I),
     1  (FIKM(K,I,T),JK(K,I,T),FIKS(K,I,T),K=1,1),I=1,SOLS)
C	WRITE (21,485) (CURI(K),K=1,1)
	WRITE (21,485) CURMI,CURIE,CURIS
	WRITE (21,490) FEVM(T),FEVS(T),QEV,(QIV(K),K=1,1)
	WRITE (21,500) (SOL(I),FEKM(I,T),FEKS(I,T),QE(I),
     1  (QI(K,I),K=1,1),I=1,SOLS)
C	WRITE (21,505) CURE
	WRITE (21,505) CURME,CURES
C
C
  398   FORMAT(//,14X,'L',8X,'CHVL',7X,'CLVLA',
     1  9X,'AE',/,6X,4D12.4)
  400   FORMAT(//,16X,'VM',10X,'VE',9X,'VIA',
     1  10X,'VS',/,8X,4F12.3)
  402   FORMAT(//,16X,'PM',10X,'PE',9X,'PIA',
     1  10X,'PS',/,8X,4F12.3)
C
  410   FORMAT(//,16X,'CM',10X,'CE',9X,'CIA',
     1  10X,'CS',/,(4X,A5,4F12.6))
  412   FORMAT(4X,A5,1X,4D12.4)
  414   FORMAT(4X,A5,4F12.6)
  420   FORMAT(5X,'IMP',1X,F12.6,12X,2F12.6)
  430   FORMAT(5X,'OSM',1X,4F12.6)
  440   FORMAT(6X,'PH',1X,4F12.6)
  442   FORMAT(4X,'CBUF',1X,24X,1F12.6,/,4X,'HCBUF',24X,1F12.6)
C
  460   FORMAT(//,15X,'XM-XS',7X,'XE-XS',
     1  6X,' XI-XS',/,(4X,A5,3F12.4))
  470   FORMAT(1H1,//,14X,
     1  'FIVM ',7X,'JV ',9X,'FIVS ',7X,/,9X,3D12.4)
  480   FORMAT(//,14X,
     1  'FIKM ',7X,'JK ',9X,'FIKS ',7X,/,(4X,A5,3D12.4))
  485   FORMAT(4X,'CUR',2X,3D12.4)
  490   FORMAT(//,15X,'FEVM',8X,'FEVS',20X,'QEV',9X,
     1  'QIV ',/,9X,2D12.4,12X,2D12.4)
  500   FORMAT(//,15X,'FEKM',8X,'FEKS',20X,'QE',10X,
     1  'QI ',/,(4X,A5,2D12.4,12X,2D12.4))
  505   FORMAT(4X,'CUR',2X,3D12.4)
C
C
C OUTPUT FOR AN AUXILLIARY FILE--E.G. PERMEABILITY CALCULATIONS
C
	FLUXV=FEVM(T)
	DO 550 I=1,SOLS
  550   FLUXS(I)=FEKM(I,T)
	DO 560 K=1,1
	FLUXV=FLUXV + FIVM(K,T)
	DO 555 I=1,SOLS
  555   FLUXS(I)=FLUXS(I) + FIKM(K,I,T)
  560	CONTINUE
C
C	WRITE (23,570) FLUXV,(FLUXS(I),I=1,SOLS)
  570   FORMAT (8E16.8)
C
	RETURN
	END
