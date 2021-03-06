	SUBROUTINE CCTRESD
C
	INTEGER SOLS,T,X,CHOP
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
	DOUBLE PRECISION OSMM(1601),OSMS(1601),
     1    TA(1601),QTA0,QTA1,ATA0,ATA1,TTIME
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
C
	CHARACTER*5 SOL(14)
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
	SOL(13)=' IMPM'
	SOL(14)='  TA '
	IS=SOLS+1
 	ATA0=10**(7.4 - 6.10142)
 	ATA1=10**(7.4 - 6.8)
C
	TTIME=0.D0
	DO 310 KX=X-CHOP,X-1
  310   TTIME=TTIME + 
     1    DX*0.5*(AM(KX,T)/FVM(KX,T) + AM(KX+1,T)/FVM(KX+1,T))
C
	DO 320 KX=X-CHOP,X
	OSMM(KX)=IMPM(KX,T)
 	QTA0=ATA0*CM(6,KX,T) - CM(4,KX,T)
 	QTA1=(ATA1*CM(8,KX,T) - CM(7,KX,T))/(ATA1 + 1.D0)
 	TA(KX)= QTA1
	DO 320 I=1,SOLS
  320 	OSMM(KX)=OSMM(KX)+CM(I,KX,T)
C
	WRITE (10,398) TTIME
  398   FORMAT (1H1,//,5X,'TTIME=',F15.4)
C
	WRITE (10,400) (SOL(I),I=1,SOLS-1),SOL(13)
  400   FORMAT(//,2X,'DIST',4X,'VM',7X,'PM',3X,'CM:',
     1  13(A5,4X))
	DIST=DIST-TL-DX
	DO 410 KX=X-CHOP,X
	DIST=DIST+DX
  410   WRITE (10,420) DIST,
     1  VM(KX,T),PM(KX,T),(CM(I,KX,T),I=1,SOLS-1),IMPM(KX,T)
  420   FORMAT (F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,D9.3,2F9.6)
C
	WRITE (10,430) (SOL(I),I=1,SOLS)
  430   FORMAT(//,2X,'DIST',4X,'TA',7X,'PHM',2X,
     1  'XM:',12(A5,4X))
	DIST=DIST-TL-DX
	DO 440 KX=X-CHOP,X
	DIST=DIST+DX
  440   WRITE (10,450) DIST,TA(KX),LCHM(KX),
     1    (XM(I,KX)-XS(I,KX),I=1,SOLS)
  450   FORMAT (F7.4,F9.6,13F9.4)
C
	WRITE (10,460) (SOL(I),I=1,SOLS-1),SOL(13),SOL(14)
  460   FORMAT(//,2X,'DIST',4X,'OSMM',4X,'FVM',3X,'FKM:',
     1  13(A5,4X))
	DIST=DIST-TL-DX
	DO 470 KX=X-CHOP,X
	DIST=DIST+DX
  470   WRITE (10,480) DIST,OSMM(KX),
     1  1.D3*FVM(KX,T),(1.D3*FKM(I,KX,T),I=1,SOLS-1),
     1  1.D3*FKM(13,KX,T),1.D3*FVM(KX,T)*TA(KX)
  480   FORMAT (F7.4,6F9.6,D9.3,4F9.6,D9.3,3F9.6)
C
	DO 520 KX=X-CHOP,X
	OSMS(KX)=IMPS(KX,T)
	DO 520 I=1,SOLS
  520 	OSMS(KX)=OSMS(KX)+CS(I,KX,T)
C
	WRITE (11,598) TIME
  598   FORMAT (1H1,//,5X,'TIME=',F15.4)
C
	WRITE (11,600) (SOL(I),I=1,SOLS-1),SOL(13)
  600   FORMAT(//,2X,'DIST',4X,'VS',7X,'PS',3X,'CS:',
     1  13(A5,4X))
	DIST=DIST-TL-DX
	DO 610 KX=X-CHOP,X
	DIST=DIST+DX
  610   WRITE (11,620) DIST,
     1  VS(KX,T),PS(KX,T),(CS(I,KX,T),I=1,SOLS-1),IMPS(KX,T)
  620   FORMAT (F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,D9.3,2F9.6)
C
	WRITE (11,630) (SOL(I),I=1,SOLS)
  630   FORMAT(//,2X,'DIST',3X,'OSMS',6X,'PHS',2X,
     1  'XS:',12(A5,4X))
	DIST=DIST-TL-DX
	DO 640 KX=X-CHOP,X
	DIST=DIST+DX
  640   WRITE (11,650) DIST,OSMS(KX),LCHS(KX),(XS(I,KX),I=1,SOLS)
  650   FORMAT (F7.4,14F9.4)
C
	RETURN
	END
