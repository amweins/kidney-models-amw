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
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C
	CHARACTER*5 SOL(12)
C
	DOUBLE PRECISION 
     1    PIKM(3,12),PJK(3,12),PIKS(3,12),
     1    CIKM(3,12),CJK(3,12),CIKS(3,12),
     1    GIKM(3,12),GJK(3,12),GIKS(3,12)
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
C
C
	DOUBLE PRECISION JNACL,JKCL,JNAH,JCLB,JNAP,JNAKCL
	CHARACTER*12, JCOTR(12,6)
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
	DO 510 I=1,12
	DO 510 J=1,6
  510   JCOTR(I,J)='            '
C
	DO 700 II=1,2
	IF (II.EQ.1) SCALE=1.0
	IF (II.EQ.2) SCALE=4.71D+07
C SCALE FOR A 25 MICRON TUBULE DIAMETER
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
	FEVS(T)=FEVS(T)*SCALE
	DO 520 I=1,SOLS
	FEKM(I,T)=FEKM(I,T)*SCALE
  520   FEKS(I,T)=FEKS(I,T)*SCALE
	DO 530 K=1,1
	FIVM(K,T)=FIVM(K,T)*SCALE
	JV(K,T)=JV(K,T)*SCALE
	FIVS(K,T)=FIVS(K,T)*SCALE
	DO 530 I=1,SOLS
	FIKM(K,I,T)=FIKM(K,I,T)*SCALE
	JK(K,I,T)=JK(K,I,T)*SCALE
  530   FIKS(K,I,T)=FIKS(K,I,T)*SCALE
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
	WRITE (27,470) (FIVM(K,T),JV(K,T)+FIVS(K,T),K=1,1),
     1   FEVM(T),FEVS(T),FLUXV
	WRITE (27,475) (SOL(I),
     1  (FIKM(K,I,T),JK(K,I,T)+FIKS(K,I,T),K=1,1),
     1  FEKM(I,T),FEKS(I,T),FLUXS(I),I=1,SOLS)
C
  470   FORMAT(1H1,/,14X,'FIVMA',7X,'FIVSA',7X,
     1  'FEVM',8X,'FEVS',8X,'FLUXV',/,9X,5D12.4)
  475   FORMAT(//,14X,'FIKMA',7X,'FIKSA',7X,
     1  'FEKM',8X,'FEKS',8X,'FLUXS',/,(4X,A5,5D12.4))
C
C
C NET COTRANSPORTERS
C
	DO 110 K=1,1
	DO 110 I=1,SOLS
	CIKM(K,I)=0.
	CIKS(K,I)=0.
  110   CJK(K,I)=0.
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	CIKM(K,I)=CIKM(K,I)+LMI(K,I,J)*AI0(K)*(XM(J)-XI(K,J))*SCALE
	CIKS(K,I)=CIKS(K,I)+LIS(K,I,J)*AI0(K)*(XI(K,J)-XS(J))*SCALE
  120   CJK(K,I)=CJK(K,I)+LIS(K,I,J)*AIE(K)*(XI(K,J)-XE(J))*SCALE
C
C THIS PORTION OF THE CODE BACK-CALCULATES THE INDIVIDUAL COTRANSPORTERS OF THE CELL
C
	JNACL=CIKM(1,1)
	JNAH= CIKS(1,12)+CJK(1,12)
	JCLB= CIKS(1,4)+CJK(1,4)
	JNAP= CIKS(1,7)+CJK(1,7)
	JNAKCL=CIKS(1,1)+CJK(1,1) + JNAH - 2.*JNAP
	JKCL= CIKS(1,2)+CJK(1,2) - JNAKCL
C
  125   FORMAT(D12.4)
	WRITE (JCOTR(1,1),125) JNACL
	WRITE (JCOTR(3,1),125) JNACL
	WRITE (JCOTR(2,2),125) JKCL
	WRITE (JCOTR(3,2),125) JKCL
	WRITE (JCOTR(1,3),125) -JNAH
	WRITE (JCOTR(12,3),125) JNAH
	WRITE (JCOTR(3,4),125) -JCLB
	WRITE (JCOTR(4,4),125) JCLB
	WRITE (JCOTR(1,5),125) 2.*JNAP
	WRITE (JCOTR(7,5),125) JNAP
	WRITE (JCOTR(1,6),125) JNAKCL
	WRITE (JCOTR(2,6),125) JNAKCL
	WRITE (JCOTR(3,6),125) 2.*JNAKCL
C
C PUMPS
C
	DO 130 K=1,1
	DO 130 I=1,SOLS
	PIKM(K,I)=0.
	PIKS(K,I)=0.
  130   PJK(K,I)=0.
C
	DO 140 K=1,1
	DO 140 I=1,SOLS
	PIKM(K,I)=PIKM(K,I)+AI0(K)*ATMI(K,I)*SCALE
	PIKS(K,I)=PIKS(K,I)+AI0(K)*ATIS(K,I)*SCALE
  140   PJK(K,I)=PJK(K,I)+AIE(K)*ATIS(K,I)*SCALE
C
	DO 150 K=1,1
	DO 150 I=1,SOLS
	GIKM(K,I)=FIKM(K,I,T)-PIKM(K,I)-CIKM(K,I)
	GIKS(K,I)=FIKS(K,I,T)-PIKS(K,I)-CIKS(K,I)
  150   GJK(K,I)=JK(K,I,T)-PJK(K,I)-CJK(K,I)
C
	WRITE (27,480) (SOL(I),
     1  (CIKM(K,I),CJK(K,I)+CIKS(K,I),K=1,1),
     1  (JCOTR(I,J),J=1,6),I=1,SOLS)
	WRITE (27,490) (SOL(I),
     1  (PIKM(K,I),PJK(K,I)+PIKS(K,I),K=1,1),I=1,SOLS)
	WRITE (27,495) (SOL(I),
     1  (GIKM(K,I),GJK(K,I)+GIKS(K,I),K=1,1),I=1,SOLS)
C
  480   FORMAT(/,14X,'CIKMA',7X,'CIKSA',7X,
     1   'JNACL',7X,'JKCL',8X,'JNAH',8X,'JCLB',8X,'JNAP',8X,
     1   'JNAKCL',/,(4X,A5,2D12.4,6A12))
  490   FORMAT(/,14X,'PIKMA',7X,'PIKSA',7X,/,(4X,A5,2D12.4))
  495   FORMAT(/,14X,'GIKMA',7X,'GIKSA',7X,/,(4X,A5,2D12.4))
C
C
	FEVM(T)=FEVM(T)/SCALE
	FEVS(T)=FEVS(T)/SCALE
	DO 620 I=1,SOLS
	FEKM(I,T)=FEKM(I,T)/SCALE
  620   FEKS(I,T)=FEKS(I,T)/SCALE
	DO 630 K=1,1
	FIVM(K,T)=FIVM(K,T)/SCALE
	JV(K,T)=JV(K,T)/SCALE
	FIVS(K,T)=FIVS(K,T)/SCALE
	DO 630 I=1,SOLS
	FIKM(K,I,T)=FIKM(K,I,T)/SCALE
	JK(K,I,T)=JK(K,I,T)/SCALE
  630   FIKS(K,I,T)=FIKS(K,I,T)/SCALE
C
  700	CONTINUE
	RETURN
	END
