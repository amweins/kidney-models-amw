	SUBROUTINE IMDRESD
C
C
	INTEGER SOLS,T,X,CHOP
C
	DOUBLE PRECISION
     1   Z(15),ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15),
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
C 
	DOUBLE PRECISION
     1   VC(901,2),PC(901,2),CC(31,901,2),IMPC(901,2),
     1   LCHC(901),XC(15,901),
     1   VS(901,2),PS(901,2),CS(16,901,2),
     1   IMPS(901,2),LCHS(901),XS(15,901),
     1   CCS(15),SAT(901,2),CCN(901,2),CCX(901,2),CCY(901,2),
     1   SC(901,2),AC(901,2),FVC(901,2),FKC(31,901,2),
     1   QV,QC(31),FVCS(901,2),FKCS(16,901,2),
     1   FVCSE(901,2),FKCSE(16,901,2),FVCSI(901,2),FKCSI(16,901,2),
     1   HCTC(901,2),FBC(901,2),POX(901,2)
C
C
	DOUBLE PRECISION CURC(901),OSMC(901),OSMS(901),
     1    TA(901),FTA(901),NAE(901),FNAE(901),BETA(10)
        CHARACTER*8 SOL(31)
C
C
	COMMON/PAR/ SOLS,T,X,CHOP,
     1   Z,ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS,LCS,HCS,
     1   LPCSE,SCSE,LCSE,HCSE,
     1   LPCSI,SCSI,LCSI,HCSI,
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
	COMMON/VAR/
     1   VC,PC,CC,VS,PS,CS,IMPC,IMPS,
     1   LCHC,LCHS,XC,XS,CCS,
     1   SAT,CCN,CCX,CCY,HCTC,FBC,POX,
     1   SC,AC,FVC,FKC,QC,QV,FVCS,FKCS,
     1   FVCSE,FKCSE,FVCSI,FKCSI
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
	SOL(11)=' NH4+'
	SOL(12)='  H+ '
	SOL(13)='HCO2-'
	SOL(14)='H2CO2'
	SOL(15)=' GLUC'
	SOL(16)=' PROT'
C
	SOL(17)='   ONH2 '
	SOL(18)='   NH2  '
	SOL(19)='   ONH3+'
	SOL(20)='   NH3+ '
	SOL(21)=' ONHCO2-'
	SOL(22)='  NHCO2-'
	SOL(23)=' OXH '
	SOL(24)='XH   '
	SOL(25)='OX-  '
	SOL(26)='X-   '
	SOL(27)='OYH  '
	SOL(28)='YH   '
	SOL(29)='OY-  '
	SOL(30)='Y-   '
	SOL(31)=' IMP '
C
C CALCULATION OF TITRATABLE ACID FOR THE 9 IMPORTANT BUFFER REACTIONS:
C  PHOSPHATE PLUS 8 FOR HGB:
C
C	17<-11      ONH2        KZO
C       18<-12      NH2         KZD
C	19<-13      ONH3+
C	20<-14      NH3+ 
C	21<-15      ONHCO2-     KCO
C	22<-16      NHCO2-      KCD
C	23<-17      OXH 
C	24<-18      XH   
C	25<-19      OX-         KXO
C	26<-20      X-          KXD
C	27<-21      OYH  
C	28<-22      YH   
C	29<-23      OY-         KYO
C	30<-24      Y-          KYD
C
 	BETA(1)=10**(7.4 - 6.8)
        BETA(2)=KXO*(10**7.40)
        BETA(3)=KXD*(10**7.40)
        BETA(4)=KYO*(10**7.40)
        BETA(5)=KYD*(10**7.40)
        BETA(6)=KZO*(10**7.40)
        BETA(7)=KZD*(10**7.40)
        BETA(8)=KCO*(10**7.40)
        BETA(9)=KCD*(10**7.40)
C
	DO 320 KX=X-CHOP,X
	CURC(KX)=0.D0
	OSMC(KX)=0.D0
	OSMS(KX)=0.D0
C	OSMC(KX)=IMPC(KX,T)/RT
 	TA(KX)=(BETA(1)*CC(8,KX,T) - CC(7,KX,T))/(BETA(1) + 1.D0)+
     1   (BETA(2)*CC(23,KX,T) - CC(25,KX,T))/(BETA(2) + 1.D0)+
     1   (BETA(3)*CC(24,KX,T) - CC(26,KX,T))/(BETA(3) + 1.D0)+
     1   (BETA(4)*CC(27,KX,T) - CC(29,KX,T))/(BETA(4) + 1.D0)+
     1   (BETA(5)*CC(28,KX,T) - CC(30,KX,T))/(BETA(5) + 1.D0)+
     1   (BETA(6)*CC(19,KX,T) - CC(17,KX,T))/(BETA(6) + 1.D0)+
     1   (BETA(7)*CC(20,KX,T) - CC(18,KX,T))/(BETA(7) + 1.D0)+
     1   (CC(6,KX,T)*BETA(8)*CC(17,KX,T) - CC(21,KX,T))/
     1         (CC(6,KX,T)*BETA(8) + 1.D0)+
     1   (CC(6,KX,T)*BETA(9)*CC(18,KX,T) - CC(22,KX,T))/
     1         (CC(6,KX,T)*BETA(9) + 1.D0)
C
        FTA(KX)=FVC(KX,T)*TA(KX)
        NAE(KX)=TA(KX)+CC(11,KX,T)-CC(4,KX,T)
        FNAE(KX)=FVC(KX,T)*NAE(KX)
	DO 320 I=1,SOLS
        CURC(KX)=CURC(KX)+F*Z(I)*FKCS(I,KX,T)
        OSMC(KX)=OSMC(KX)+CC(I,KX,T)
        OSMS(KX)=OSMS(KX)+CS(I,KX,T)
  320   CONTINUE
C
	WRITE (11,398) TIME
  398   FORMAT (1H1,//,5X,'TIME=',F15.4)
C
	WRITE (11,490) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15)
  490   FORMAT(//,5X,'DIST',8X,'PS',4X,'CS:',10(A5,6X),1X,
     1    'PHS',7X,'OSMS')
	DIST = DIST - CL/FLOAT(CHOP)
	DO 405 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  405   WRITE (11,495)  DIST, PS(KX,T), 
     1   (CS(I,KX,T),I=1,4),(CS(I,KX,T),I=7,9),
     1   CS(11,KX,T),CS(13,KX,T),CS(15,KX,T),LCHS(KX),OSMS(KX)
  495   FORMAT (3X,F8.6,13F11.6)
C
	WRITE (11,400) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15),SOL(16)
  400   FORMAT(//,5X,'DIST',7X,'PC',5X,'CC:', 11(A5,6X),'IMPC')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 410 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  410   WRITE (11,420) DIST,PC(KX,T),
     1   (CC(I,KX,T),I=1,4),(CC(I,KX,T),I=7,9),
     1   CC(11,KX,T),CC(13,KX,T),CC(15,KX,T),CC(16,KX,T),IMPC(KX,T)
  420   FORMAT (3X,F8.6,13F11.6)
C
	WRITE (11,430) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15)
  430   FORMAT(//,5X,'DIST',7X,'VC',5X,
     1  'XC:',10(A5,6X),1X,'OSMC')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 440 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  440   WRITE (11,450) DIST,VC(KX,T),
     1  (XC(I,KX)-XS(I,KX),I=1,4),(XC(I,KX)-XS(I,KX),I=7,9),
     1  XC(11,KX)-XS(11,KX),XC(13,KX)-XS(13,KX),
     1  XC(15,KX)-XS(15,KX),OSMC(KX)
  450   FORMAT (3X,F8.6,12F11.6)
C
C
	WRITE (11,460) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15)
  460   FORMAT(1H1,//,5X,'DIST',6X,'FVC',4X,'FKC:',
     1  10(A5,6X),' FBC ',7X,'HCTC')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 470 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  470   WRITE (11,480) DIST, 1.D6*FVC(KX,T),
     1   (1.D6*FKC(I,KX,T),I=1,4),(1.D6*FKC(I,KX,T),I=7,9),
     1   1.D6*FKC(11,KX,T),1.D6*FKC(13,KX,T),1.D6*FKC(15,KX,T),
     1   1.D6*FBC(KX,T),HCTC(KX,T)
  480   FORMAT (3X,F8.6,13F11.6)
C
	WRITE (11,464) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15)
  464   FORMAT(/,11X,'APR:',2X,'VOL',6X,10(A5,6X),'  TA ',6X,' NAE ')
	WRITE (11,484) 
     1  1.D6*(FVC(1,T)-FVC(CHOP+1,T)), 
     1  (1.D6*(FKC(I,1,T)-FKC(I,CHOP+1,T)), I=1,4),
     1  (1.D6*(FKC(I,1,T)-FKC(I,CHOP+1,T)), I=7,9),
     1  1.D6*(FKC(11,1,T)-FKC(11,CHOP+1,T)), 
     1  1.D6*(FKC(13,1,T)-FKC(13,CHOP+1,T)), 
     1  1.D6*(FKC(15,1,T)-FKC(15,CHOP+1,T)), 
     1  1.D6*(FTA(1)-FTA(CHOP+1)),1.D6*(FNAE(1)-FNAE(CHOP+1))
  484   FORMAT (11X,13F11.6)
C
	WRITE (11,462) (SOL(I),I=1,4),(SOL(I),I=7,9),
     1   SOL(11),SOL(13),SOL(15)
  462   FORMAT(//,5X,'DIST',6X,'FVCS',2X,'FKCS:',
     1  10(A5,6X),2X,'CUR')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 472 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  472   WRITE (11,482) DIST,1.D6*FVCS(KX,T),
     1  (1.D6*FKCS(I,KX,T), I=1,4),(1.D6*FKCS(I,KX,T), I=7,9),
     1  1.D6*FKCS(11,KX,T), 1.D6*FKCS(13,KX,T), 1.D6*FKCS(15,KX,T), 
     1  CURC(KX)*1.D6
  482   FORMAT (3X,F8.6,1P,12D11.3)
C
C
	WRITE (11,598) CCN(1,T),CCX(1,T),CCY(1,T)
  598   FORMAT (1H1,//,27X,'CCN=',F11.6,5X,'CCX=',F11.6,5X,'CCY=',F11.6)
C
	WRITE (11,600) SOL(12),(SOL(I),I=17,22)
  600   FORMAT(//,5X,'DIST',7X,'PHC',6X,A5,2X,'HGB:',6(A8,3X),
     1   12X,'  TA ',6X,' NAE ',6X,'O2-SAT')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 610 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  610   WRITE (11,620) DIST,LCHC(KX),CC(12,KX,T),
     1  (CC(I,KX,T),I=17,22),TA(KX),NAE(KX),SAT(KX,T)
  620   FORMAT (3X,F8.6,F10.6,D12.4,6F11.6,11X,3F11.6)
C
	WRITE (11,630) (SOL(I),I=SOLS+8,SOLS+15)
  630   FORMAT(//,5X,'DIST',4X,'HGB:',8(A5,6X),10X,
     1   ' FTA ',6X,' FNAE')
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 640 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  640   WRITE (11,650) DIST,(CC(I,KX,T),I=23,30),
     1   1.D6*FTA(KX),1.D6*FNAE(KX)
  650   FORMAT (3X,F8.6,8F11.6,11X,2F11.6)
C
C
	DIST = DIST - CL - CL/FLOAT(CHOP)
	DO 710 KX=1,CHOP+1
	DIST = DIST + CL/FLOAT(CHOP)
  710   WRITE (12,715) DIST,1.D6*FVC(KX,T),
     1   PS(KX,T)-PC(KX,T)+IMPC(KX,T),
     1   1.D3*(OSMC(KX)-OSMS(KX)),1.D6*FVCS(KX,1),
     1   1.D6*FVCSE(KX,1), 1.D6*FVCSI(KX,1)
  715   FORMAT(8D16.8)
	RETURN
	END
