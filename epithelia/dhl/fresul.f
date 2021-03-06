	SUBROUTINE FRESUL
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
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),ATIE(3,12),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3),JHK(3),JHP(3),JAE1(3),JTSC(3),JNHE3(3,3),
     1   NNKCC(3),NKCL(3),JNKCC(3,4),JKCC(3,3)
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION 
     1    PIKM(3,12),PJK(3,12),PIKS(3,12),
     1    CIKM(3,12),CJK(3,12),CIKS(3,12),
     1    GIKM(3,12),GJK(3,12),GIKS(3,12),
     1    SNKCC(8),VNKCC(4),SKCC(6),VKCCJK(3),VKCCIS(3),VKCC(3),
     1    VNHE3(3),SNHE3(3,2)
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
C
	DOUBLE PRECISION FLUXV,FLUXS(12),SCALE
C
C	FLUXV-	TOTAL EPITHELIAL VOLUME FLUX
C	FLUXS-	TOTAL EPITHELIAL SOLUTE FLUX
C	SCALE-  FACTOR TO TRANSLATE FROM MMOL/S*CM2 TO PMOL/MIN*MM
C		=PI*DIAM*(1.0E+09)*60/10
C			DIAM = 0.0030 cm  => SCALE = 5.65E+07
C			DIAM = 0.0016 cm  => SCALE = 3.02E+07
C
	DOUBLE PRECISION JKCL,JNAHM,JNAHS,JCLB,JNAP,JNAKCL,JNAB
	CHARACTER*12, JCOTR(12,12)
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
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,
     1   NNKCC,NKCL,JNKCC,JKCC
C
	DO 510 I=1,12
	DO 510 J=1,12
  510   JCOTR(I,J)='            '
C
	DO 700 II=1,2
	IF (II.EQ.1) SCALE=1.0
C SCALE FOR A 30 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=5.65D+07
C SCALE FOR A 20 MICRON TUBULE DIAMETER
	IF (II.EQ.2) SCALE=3.77D+07
C SCALE FOR A 16 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=3.02D+07
C SCALE FOR A 15 MICRON TUBULE DIAMETER
C	IF (II.EQ.2) SCALE=2.83D+07
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
	FEVM(T)=FEVM(T)*SCALE
	DO 520 I=1,SOLS
  520   FEKM(I,T)=FEKM(I,T)*SCALE
	DO 530 K=1,1
	FIVM(K,T)=FIVM(K,T)*SCALE
	DO 530 I=1,SOLS
  530   FIKM(K,I,T)=FIKM(K,I,T)*SCALE
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
C OUTPUT THE SUMMARY FLUXES
C
	WRITE (27,470) (FIVM(K,T),K=1,1),FEVM(T),FLUXV
	WRITE (27,475) (SOL(I), (FIKM(K,I,T),K=1,1),
     1   FEKM(I,T),FLUXS(I),I=1,SOLS)
C
  470   FORMAT(1H1,/,14X,'FIVMA',7X,'FEVM',8X,'FLUXV',/,9X,3D12.4)
  475   FORMAT(//,14X,'FIKMA',7X,
     1  'FEKM',8X,'FLUXS',/,(4X,A5,3D12.4))
C
C
C NET COTRANSPORTERS
C
	DO 110 K=1,1
	DO 110 I=1,SOLS
  110   CIKM(K,I)=0.
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
  120   CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AI0(K)*(XM(J)-XS(J))*SCALE
C
C
	DO 150 K=1,1
	DO 150 I=1,SOLS
  150   GIKM(K,I)=FIKM(K,I,T)-CIKM(K,I)
C
	WRITE (27,480) (SOL(I),(CIKM(K,I),K=1,1),I=1,SOLS)
	WRITE (27,495) (SOL(I), (GIKM(K,I),K=1,1),I=1,SOLS)
C
  480   FORMAT(/,14X,'CIKMA',/,(4X,A5,D12.4))
  495   FORMAT(/,14X,'GIKMA',/,(4X,A5,D12.4))
C
	GO TO 501
	IF (II.EQ.2) WRITE (23,500) 1.D3*CM(2),
     1    FLUXS(1),FEKM(1,T),FIKM(1,1,T),FLUXS(2),
     1    FLUXS(3),FEKM(3,T),FIKM(1,3,T)
  500   FORMAT (8D16.8)
  501   CONTINUE
C
	FEVM(T)=FEVM(T)/SCALE
	DO 620 I=1,SOLS
  620   FEKM(I,T)=FEKM(I,T)/SCALE
	DO 630 K=1,1
	FIVM(K,T)=FIVM(K,T)/SCALE
	DO 630 I=1,SOLS
  630   FIKM(K,I,T)=FIKM(K,I,T)/SCALE
C
  700	CONTINUE
	RETURN
	END
