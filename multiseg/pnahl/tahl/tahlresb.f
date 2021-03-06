	SUBROUTINE TAHLRESB
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
C
C DERIVED QUANTITIES ARE COMPUTED:OSMOLALITY,CURRENT FLOW,AND THE
C NET QUANTITY "Q(I)" OF SOLUTE INTRODUCED INTO THE COMPARTMENT
C
	OSMM=IMPM(X,T)
	OSMS=IMPS(X,T)
	DO 320 I=1,SOLS
	OSMM=OSMM+CM(I,X,T)
	OSMS=OSMS+CS(I,X,T)
  320   CONTINUE
C
C
	WRITE (91,396) TIME,DIST,AM(X,T),SM(X,T)
  396   FORMAT (1H1,/,5X,'TIME=',F15.4,7X,'DIST=',F8.6,
     1  7X,'AM=',D10.4,7X,'SM=',D10.4)
C
C
	WRITE (91,400) VM(X,T),VS(X,T)
	WRITE (91,402) PM(X,T),PS(X,T)
C
	WRITE (91,410) (SOL(I),
     1  CM(I,X,T),CS(I,X,T),XM(I,X)-XS(I,X),I=1,4)	
	WRITE (91,412) (SOL(I),
     1  CM(I,X,T),CS(I,X,T),XM(I,X)-XS(I,X),I=5,5)	
	WRITE (91,414) (SOL(I),
     1  CM(I,X,T),CS(I,X,T),XM(I,X)-XS(I,X),I=6,SOLS-3)	
	WRITE (91,412) (SOL(I),
     1  CM(I,X,T),CS(I,X,T),XM(I,X)-XS(I,X),I=SOLS-2,SOLS-2)	
	WRITE (91,414) (SOL(I),
     1  CM(I,X,T),CS(I,X,T),XM(I,X)-XS(I,X),I=SOLS-1,SOLS-1)	
	WRITE (91,412) SOL(SOLS),
     1  CM(SOLS,X,T),CS(SOLS,X,T),XM(SOLS,X)-XS(SOLS,X)
	WRITE (91,420) IMPM(X,T),IMPS(X,T)
	WRITE (91,430) OSMM,OSMS
	WRITE (91,440) LCHM(X),LCHS(X)
	WRITE (91,442) 
C
	WRITE (91,470) (FIVM(K,X,T),K=1,1),FEVM(X,T)
	WRITE (91,480) (SOL(I),(FIKM(K,I,X,T),K=1,1),
     1   FEKM(I,X,T),I=1,SOLS)
	WRITE (91,485) (CURI(K,X),K=1,1),CURE(X)
C
C
  400   FORMAT(//,16X,'VM',10X,'VS',/,8X,2F12.3)
  402   FORMAT(//,16X,'PM',10X,'PS',/,8X,2F12.3)
C
  410   FORMAT(//,16X,'CM',10X,'CS',9X,'XM-XS',/,(4X,A5,3F12.6))
  412   FORMAT(4X,A5,1X,2D12.4,F11.6)
  414   FORMAT(4X,A5,3F12.6)
  420   FORMAT(5X,'IMP',1X,2F12.6)
  430   FORMAT(5X,'OSM',1X,2F12.6)
  440   FORMAT(6X,'PH',1X,2F12.6)
  442   FORMAT(//)
C
  470   FORMAT(1H1,//,14X,'FIVM',7X,'FEVM',/,9X,2D12.4)
  480   FORMAT(//,14X,'FIKM',7X,'FEKM',/,(4X,A5,2D12.4))
  485   FORMAT(4X,'CUR',2X,2D12.4)
C
	RETURN
	END
