	SUBROUTINE ARESUL
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
     1   LHP(3),XIHP(3),XHP(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
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
     1   LHP,XIHP,XHP,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
C
	REAL 
     1    HAME(12),HAMI(3,12),HAIE(3,12),HAES(12),HAIS(3,12),
     1    GME(12),GMI(3,12),GIE(3,12),GES(12),GIS(3,12),
     1    TGME,TGMI(3),TGIE(3),TGES,TGIS(3)
C
C 	G.(I)-	AREA ADJUSTED PARTIAL CONDUCTANCE (MHO/CM2)
C 	TG.-	AREA ADJUSTED TOTAL CONDUCTANCE (MHO/CM2)
C
	CHARACTER*5 SOL(12),CELL(3)
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
	CELL(1)='PRINC'
	CELL(2)='ALPHA'
	CELL(3)='BETA '
C
C THIS IS A STEADY STATE SOLUTION
C     OR THE FIRST PASS OF A TIMED EXPERIMENT
C
	TGME=0.D0
	TGES=0.D0
	DO 355 I=1,SOLS
	HAME(I)=HME(I)*AME
	HAES(I)=HES(I)*AE(T)
C
	GME(I)=1.D-6*((F*Z(I))**2)*HAME(I)*CME(I)/RTE
	GES(I)=1.D-6*((F*Z(I))**2)*HAES(I)*CES(I)/RTE
C
	TGME=TGME+GME(I)
	TGES=TGES+GES(I)
  355   CONTINUE
C
C
	DO 370 K=1,1
	TGMI(K)=0.D0
	TGIE(K)=0.D0
	TGIS(K)=0.D0
C
	DO 365 I=1,SOLS
	HAMI(K,I)=HMI(K,I)*AI0(K)
	HAIE(K,I)=HIS(K,I)*AIE(K)
	HAIS(K,I)=HIS(K,I)*AI0(K)
C
	GMI(K,I)=1.D-6*((F*Z(I))**2)*HAMI(K,I)*CMI(K,I)/RTE
	GIE(K,I)=1.D-6*((F*Z(I))**2)*HAIE(K,I)*CIE(K,I)/RTE
	GIS(K,I)=1.D-6*((F*Z(I))**2)*HAIS(K,I)*CIS(K,I)/RTE
C
	TGMI(K)=TGMI(K)+GMI(K,I)
	TGIE(K)=TGIE(K)+GIE(K,I)
	TGIS(K)=TGIS(K)+GIS(K,I)
  365   CONTINUE
  370   CONTINUE
C
	WRITE (21,372) RTE,F,PKC,PKP,PKN,KHY(4),KDHY(4)
	WRITE (21,374) L0,AME,AE0,MUA,CHVL0,MUV
	WRITE (21,376)
	DO 375 K=1,1
  375   WRITE (21,378) CELL(K),AIE(K),AI0(K),CLVL0(K),
     1  IMP0(K),ZIMP(K),TBUF(K),PKB(K),KHY(K),KDHY(K)
C
  372   FORMAT (1H1,//,5X,'RTE=',F4.2,7X,
     1  'F=',D10.2,5X,'PKC=',F4.2,7X,'PKP=',F4.2,7X,'PKN=',F4.2,7X,
     1  'KHY=',D11.3,5X,'KDHY=',D11.3)
  374   FORMAT (/,5X,'L0=',D10.2,5X,'AME=',D10.2,5X,'AE0=',D10.2,5X,
     1  'MUA=',D10.2,5X,'CHVL0=',D10.2,5X,'MUV=',D10.2)
  376   FORMAT (/,5X,'CELL',6X,'AIE',10X,'AI0',9X,'CLVL0',7X,
     1  'IMP0',6X,'ZIMP',6X,'TBUF',6X,'PKB',8X
     1  'KHY',7X,'KDHY')
  378   FORMAT (5X,A5,D12.4,3X,F6.4,4X,D12.4,3X,
     1  F6.4,4X,F6.3,4X,F6.4,4X,F6.4,2D11.3)
C
	WRITE (21,380)
  380   FORMAT(//)
C
	WRITE (21,381) LPME*AME,LPES*AE(T),
     1    (LPMI(K)*AI0(K),LPIS(K)*AIE(K),LPIS(K)*AI0(K),K=1,1)
	WRITE (21,384)
	WRITE (21,385) (SOL(I),
     1  SME(I),SES(I),(SMI(K,I),SIS(K,I),SIS(K,I),K=1,1),I=1,SOLS)
	WRITE (21,386)
	WRITE (21,385) (SOL(I),
     1  HAME(I),HAES(I),(HAMI(K,I),HAIE(K,I),HAIS(K,I),K=1,1),I=1,SOLS)
	WRITE (21,388)
	WRITE (21,385) (SOL(I),
     1  GME(I),GES(I),(GMI(K,I),GIE(K,I),GIS(K,I),K=1,1),I=1,SOLS)
	WRITE (21,389) TGME,TGES,(TGMI(K),TGIE(K),TGIS(K),K=1,1)
	WRITE (21,390) 1./TGME,1./TGES,
     1  (1./TGMI(K),1./TGIE(K),1./TGIS(K),K=1,1)
C
  381   FORMAT (/,14X,'LPME',7X,'LPES',
     1  6X,' LPMI',6X,' LPIE',6X,' LPIS',/,9X,5D11.4)
  384   FORMAT (//,15X,'SME',8X,'SES',
     1  7X,' SMI',7X,' SIE',7X,' SIS')
  385   FORMAT (3X,A5,1X,5D11.4)
  386   FORMAT (//,15X,'HME',8X,'HES',
     1  7X,' HMI',7X,' HIE',7X,' HIS')
  388   FORMAT (//,15X,'GME',8X,'GES',
     1  7X,' GMI',7X,' GIE',7X,' GIS')
  389   FORMAT (/,3X,'TOTAL',1X,5D11.4)
  390   FORMAT (/,3X,'RESIS',1X,5D11.4)
C
C
 	RETURN
	END
