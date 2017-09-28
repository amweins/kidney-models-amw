	SUBROUTINE DHLRESF(CHOP)
C A NEW SUBROUTINE THAT OUTPUTS SOLUTE FLUXES POINTWISE - SORTED BY SPECIES
C BASED ON FRESUL IN DHL DIRECTORY
C
	INTEGER SOLS,X,T,CHOP
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
     1   LMI(3,12,12),LIS(3,12,12),
     1   ATMI(3,12,81),ATIS(3,12,81),ATIE(3,12,81),
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2),
     1   NNKCC(3),MKCC(3),JNKCC(3,4,81,2),JKCC(3,3,81,2)
C
	DOUBLE PRECISION 
     1    PIKM(3,12,81),PJK(3,12,81),PIKS(3,12,81),
     1    CIKM(3,12,81),CJK(3,12,81),CIKS(3,12,81),
     1    GIKM(3,12,81),GJK(3,12,81),GIKS(3,12,81)
C
C	PIKM-	ACTIVE TRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	PJK-	ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	PIKS-	ACTIVE TRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	CIKM-	COTRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
C	CJK-	COTRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
C	CIKS-	COTRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
C	GIKM-	CHANNEL FLUX OF SOLUTE ACROSS APICAL MEMBRANE
C	GJK-	CHANNEL FLUX OF SOLUTE ACROSS LATERAL MEMBRANE	
C	GIKS-	CHANNEL FLUX OF SOLUTE ACROSS BASAL MEMBRANE
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION FLUXV(81),FLUXS(12,81),SCALE
C
C	FLUXV-	TOTAL EPITHELIAL VOLUME FLUX
C	FLUXS-	TOTAL EPITHELIAL SOLUTE FLUX
C	SCALE-  FACTOR TO TRANSLATE FROM MMOL/S*CM2 TO PMOL/MIN*MM
C		=PI*DIAM*(1.0E+09)*60/10
C			DIAM = 0.0030 cm  => SCALE = 5.65E+07
C			DIAM = 0.0016 cm  => SCALE = 3.02E+07
C
	DOUBLE PRECISION JNAHS(81),JCLB(81),JNAP(81),JNAB(81)
C
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
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,
     1   NNKCC,MKCC,JNKCC,JKCC
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
C SCALE FOR A 30 MICRON TUBULE DIAMETER
C	SCALE=5.65D+07
C SCALE FOR A 20 MICRON TUBULE DIAMETER
	SCALE=3.77D+07
C SCALE FOR A 16 MICRON TUBULE DIAMETER
C	SCALE=3.02D+07
C SCALE FOR A 15 MICRON TUBULE DIAMETER
C	SCALE=2.83D+07
C
	DO 570 X=1,CHOP+1
	FEVM(X,T)=FEVM(X,T)*SCALE
	DO 520 I=1,SOLS
  520   FEKM(I,X,T)=FEKM(I,X,T)*SCALE
	DO 530 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)*SCALE
	DO 530 I=1,SOLS
  530   FIKM(K,I,X,T)=FIKM(K,I,X,T)*SCALE
C
	FLUXV(X)=FEVM(X,T)
	DO 550 I=1,SOLS
  550   FLUXS(I,X)=FEKM(I,X,T)
	DO 560 K=1,1
	FLUXV(X)=FLUXV(X) + FIVM(K,X,T)
	DO 555 I=1,SOLS
  555   FLUXS(I,X)=FLUXS(I,X) + FIKM(K,I,X,T)
  560	CONTINUE
  570	CONTINUE
C
C
	DO 170 X=1,CHOP+1
C
C NET COTRANSPORTERS
	DO 110 K=1,1
	DO 110 I=1,SOLS
  110   CIKM(K,I,X)=0.
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
  120   CIKM(K,I,X)=CIKM(K,I,X)+
     1    LMI(K,I,J)*AI0(K)*(XM(J,X)-XI(K,J,X))*SCALE
C
C CONDUCTANCE FLUXES
	DO 150 K=1,1
	DO 150 I=1,SOLS
  150   GIKM(K,I,X)=FIKM(K,I,X,T)-CIKM(K,I,X)
C
  170   CONTINUE
C
C THE SOLUTE FLUXES ARE OUTPUT ACCORDING TO SOLUTE:
C   Na, K, Cl, NH4
C
	WRITE (53,300) 
  300   FORMAT(//,1X,'DIST',2X,'Na:',3X,'FLUX',6X,'FEKM',6X,'FIKM')
	DIST=DIST-TL-DX
	DO 310 X=1,CHOP+1
	DIST=DIST+DX
  310   WRITE (53,320) DIST,FLUXS(1,X),FEKM(1,X,T),FIKM(1,1,X,T)
  320   FORMAT (F7.4,2X,3F10.4)
C
	WRITE (53,330) 
  330   FORMAT(//,1X,'DIST',2X,'K:',5X,'FLUX',6X,'FEKM',6X,'FIKM')
	DIST=DIST-TL-DX
	DO 340 X=1,CHOP+1
	DIST=DIST+DX
  340   WRITE (53,350) DIST,FLUXS(2,X),FEKM(2,X,T),FIKM(1,2,X,T)
  350   FORMAT (F7.4,2X,3F10.4)
C
	WRITE (53,360) 
  360   FORMAT(//,1X,'DIST',2X,'Cl:',3X,'FLUX',6X,'FEKM',6X,'FIKM')
	DIST=DIST-TL-DX
	DO 370 X=1,CHOP+1
	DIST=DIST+DX
  370   WRITE (53,380) DIST,FLUXS(3,X),FEKM(3,X,T),FIKM(1,3,X,T)
  380   FORMAT (F7.4,2X,3F10.4)
C
	WRITE (53,390) 
  390   FORMAT(//,1X,'DIST',2X,'NH4:',3X,'FLUX',6X,'FEKM',6X,'FIKM',
     1   6X,'JNH3')
	DIST=DIST-TL-DX
	DO 400 X=1,CHOP+1
	DIST=DIST+DX
  400   WRITE (53,410) DIST,FLUXS(11,X),FEKM(11,X,T),FIKM(1,11,X,T),
     1  FIKM(1,10,X,T)
  410   FORMAT (F7.4,2X,4F10.4)
C
	DO 670 X=1,CHOP+1
	FEVM(X,T)=FEVM(X,T)/SCALE
	DO 620 I=1,SOLS
  620   FEKM(I,X,T)=FEKM(I,X,T)/SCALE
	DO 630 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)/SCALE
	DO 630 I=1,SOLS
  630   FIKM(K,I,X,T)=FIKM(K,I,X,T)/SCALE
  670	CONTINUE
C
  700	CONTINUE
	RETURN
	END
