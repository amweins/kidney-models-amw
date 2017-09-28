	SUBROUTINE PCTRESA
C
	INTEGER SOLS,X,T
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,L0,L(81,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(81,2),PM(81,2),CM(15,81,2),
     1   IMPM(81,2),LCHM(81),XM(15,81),
     1   VS(81,2),PS(81,2),CS(15,81,2),
     1   IMPS(81,2),LCHS(81),XS(15,81),
     1   TL,DX,RM0,MUM,ETA,
     1   SM(81,2),AM(81,2),FVM(81,2),FKM(16,81,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(81,2),MUA,CHVL0,CHVL(81,2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(81,2),PE(81,2),CE(15,81,2),LCHE(81),XE(15,81),
     1   FEVM(81,2),FEKM(15,81,2),FEVS(81,2),FEKS(15,81,2),CURE(81)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),CLVL0(3),IMP0(3),CLVL(3,81,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,81,2),HCBUF(3,81,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),ATMI(3,15,81),ATIS(3,15,81),
     1   VI(3,81,2),PI(3,81,2),CI(3,15,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,15,81),FIVM(3,81,2),FIKM(3,15,81,2),FIVS(3,81,2),
     1   FIKS(3,15,81,2),CURI(3,81),JV(3,81,2),JK(3,15,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2),QIAMM
C
        COMMON SOLS,T,X,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   TL,DX,RM0,MUM,ETA,
     1   SM,AM,FVM,FKM,
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
	REAL 
     1    PFME,PFMI(3),PFIE(3),PFES,PFIS(3),
     1    HAME(15),HAMI(3,15),HAIE(3,15),HAES(15),HAIS(3,15),
     1    GME(15),GMI(3,15),GIE(3,15),GES(15),GIS(3,15),
     1    TGME,TGMI(3),TGIE(3),TGES,TGIS(3)
C
C 	G.(I)-	AREA ADJUSTED PARTIAL CONDUCTANCE (MHO/CM2)
C 	TG.-	AREA ADJUSTED TOTAL CONDUCTANCE (MHO/CM2)
C
	CHARACTER*5 SOL(15),CELL(3)
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
	CELL(1)='PRINC'
	CELL(2)='ALPHA'
	CELL(3)='BETA '
	NCELLS=1
C
C THIS IS A STEADY STATE SOLUTION
C     OR THE FIRST PASS OF A TIMED EXPERIMENT
C
	DO 100 I=1,SOLS
C
	CME(I)=(CE(I,X,T)-CM(I,X,T))
	IF (CME(I)) 147,146,147
  146   CME(I)=CM(I,X,T)
	GO TO 148
  147   CME(I)=CME(I)/DLOG(CE(I,X,T)/CM(I,X,T))
  148   CES(I)=(CE(I,X,T)-CS(I,X,T))
	IF (CES(I)) 150,149,150
  149   CES(I)=CS(I,X,T)
	GO TO 151
  150   CES(I)=CES(I)/DLOG(CE(I,X,T)/CS(I,X,T))
  151   DO 160 K=1,1
	CMI(K,I)=(CI(K,I,X,T)-CM(I,X,T))
	IF (CMI(K,I)) 153,152,153
  152   CMI(K,I)=CM(I,X,T)
	GO TO 154
  153   CMI(K,I)=CMI(K,I)/DLOG(CI(K,I,X,T)/CM(I,X,T))
  154   CIS(K,I)=(CI(K,I,X,T)-CS(I,X,T))
	IF (CIS(K,I)) 156,155,156
  155   CIS(K,I)=CS(I,X,T)
	GO TO 157
  156   CIS(K,I)=CIS(K,I)/DLOG(CI(K,I,X,T)/CS(I,X,T))
  157   CIE(K,I)=(CI(K,I,X,T)-CE(I,X,T))
	IF (CIE(K,I)) 159,158,159
  158   CIE(K,I)=CI(K,I,X,T)
	GO TO 160
  159   CIE(K,I)=CIE(K,I)/DLOG(CI(K,I,X,T)/CE(I,X,T))
  160   CONTINUE
  100   CONTINUE
C
C
	PFME=LPME*AME/.018
	PFES=LPES*AE(X,T)/.018
	DO 350 K=1,1
	PFMI(K)=LPMI(K)*AMI(K)/.018
	PFIE(K)=LPIS(K)*AIE(K)/.018
	PFIS(K)=LPIS(K)*AIS(K)/.018
  350   CONTINUE
C
	TGME=0.D0
	TGES=0.D0
	DO 355 I=1,SOLS
	HAME(I)=HME(I)*AME
	HAES(I)=HES(I)*AE(X,T)
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
	HAMI(K,I)=HMI(K,I)*AMI(K)
	HAIE(K,I)=HIS(K,I)*AIE(K)
	HAIS(K,I)=HIS(K,I)*AIS(K)
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
C
	WRITE (21,372) RTE,F,PKC,PKP,PKN,KHY(4),KDHY(4)
	WRITE (21,374) L0,AME,AE0,MUA,CHVL0,MUV
	WRITE (21,376)
	DO 375 K=1,1
  375   WRITE (21,378) CELL(K),AIE(K),AMI(K),AIS(K),CLVL0(K),
     1  IMP0(K),ZIMP(K),TBUF(K),PKB(K),KHY(K),KDHY(K),QIAMM
C
  372   FORMAT (1H1,//,5X,'RTE=',F4.2,7X,
     1  'F=',D10.2,5X,'PKC=',F4.2,7X,'PKP=',F4.2,7X,'PKN=',F4.2,7X,
     1  'KHY=',D11.3,5X,'KDHY=',D11.3)
  374   FORMAT (/,5X,'L0=',D10.2,5X,'AME=',D10.2,5X,'AE0=',D10.2,5X,
     1  'MUA=',D10.2,5X,'CHVL0=',D10.2,5X,'MUV=',D10.2)
  376   FORMAT (/,5X,'CELL',6X,'AIE',10X,'AMI',9X,'AIS',9X,
     1  'CLVL0',7X,'IMP0',6X,'ZIMP',6X,'TBUF',6X,'PKB',7X
     1  'KHY',7X,'KDHY',6X,'QIAMM')
  378   FORMAT (5X,A5,4D12.4,4X,
     1  F6.4,4X,F6.3,4X,F6.4,4X,F6.4,3D11.3)
C
C
	WRITE (21,381) PFME,PFES,(PFMI(K),PFIE(K),PFIS(K),K=1,1)
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
  381   FORMAT (/,14X,'PFME',7X,'PFES',
     1  7X,'PFMI',7X,'PFIE',7X,'PFIS',
     1  /,9X,11D11.4)
  384   FORMAT (//,15X,'SME',8X,'SES',
     1  8X,'SMI',8X,'SIE',8X,'SIS')
C  381   FORMAT (/,14X,'PFME',7X,'PFES',
C     1  7X,'PFMI',7X,'PFIE',7X,'PFIS',
C     1  6X,'APFMI',6X,'APFIE',6X,'APFIS',
C     1  6X,'BPFMI',6X,'BPFIE',6X,'BPFIS',/,9X,11D11.4)
C  384   FORMAT (//,15X,'SME',8X,'SES',
C     1  8X,'SMI',8X,'SIE',8X,'SIS',
C     1  7X,'ASMI',7X,'ASIE',7X,'ASIS',
C     1  7X,'BSMI',7X,'BSIE',7X,'BSIS')
  385   FORMAT (3X,A5,1X,5D11.4)
  386   FORMAT (//,15X,'HME',8X,'HES',
     1  8X,'HMI',8X,'HIE',8X,'HIS')
  388   FORMAT (//,15X,'GME',8X,'GES',
     1  8X,'GMI',8X,'GIE',8X,'GIS')
C  386   FORMAT (//,15X,'HME',8X,'HES',
C     1  8X,'HMI',8X,'HIE',8X,'HIS',
C     1  7X,'AHMI',7X,'AHIE',7X,'AHIS',
C     1  7X,'BHMI',7X,'BHIE',7X,'BHIS')
C  388   FORMAT (//,15X,'GME',8X,'GES',
C     1  8X,'GMI',8X,'GIE',8X,'GIS',
C     1  7X,'AGMI',7X,'AGIE',7X,'AGIS',
C     1  7X,'BGMI',7X,'BGIE',7X,'BGIS')
  389   FORMAT (/,3X,'TOTAL',1X,11D11.4)
  390   FORMAT (/,3X,'RESIS',1X,11D11.4)
C
C
 	RETURN
	END
