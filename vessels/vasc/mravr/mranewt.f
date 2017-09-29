	SUBROUTINE MRANEWT(COUNT,IVR,SWOUT)
C
	INTEGER SOLS,X,TAU,T,TLIM,COUNT,EXP,CHOP,SWOUT,SLICE
C
C PARAMETERS
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
C VARIABLES
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
	DOUBLE PRECISION RJ,TPS,TNS,TFS,QPS,QNS,QFS,POX0(20)
	DOUBLE PRECISION PHI,GAMMA,DGAM,PPHI,MAXPHI,PREMAX,
     1   DELGAM,NDERIV
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       DGAM-   VARIABLE INCREMENT TO COMPUTE NUMERICAL DERIVATIVE
C       PPHI-   ERROR VECTOR AT GAMMA + DGAM
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C       PREMAX- MAXPHI FOR THE PREVIOUS NEWTON ITERATE
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C
	DIMENSION PHI(40),GAMMA(40),PPHI(40),NDERIV(40,40),
     1   RPHI(40),DELGAM(40)
C
	DOUBLE PRECISION SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(40),PFL(40)
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
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
C
C  SOLUTE INDEX:
C	1-          NA+
C	2-          K+
C	3-          CL-
C	4-          HCO3-
C	5-          H2CO3
C	6-          CO2
C	7-          HPO4--
C	8-          H2PO4-
C	9-          UREA
C	10-         NH3
C	11-         NH4+
C	12-         H+
C	13-         HCO2-
C	14-         H2CO2
C	15-         GLUC
C	16-         PROT
C
C	17          ONH2 
C       18          NH2  
C	19          ONH3+
C	20          NH3+ 
C	21          ONHCO2-
C	22          NHCO2-
C	23          OXH 
C	24          XH   
C	25          OX-  
C	26          X-   
C	27          OYH  
C	28          YH   
C	29          OY-  
C	30          Y-   
C
	PTOL=1.D-20
	M1=40
	SLICE=1
C
C START THE VR CALCULATION
C
        IF(SWOUT.EQ.2) GO TO 57
	CALL MRAPSET(TAU,TLIM)
	DX=0.5D0*CL/DBLE(FLOAT(CHOP))
C
C BOUNDARY VALUES,ACTIVE TRANSPORT,AND THE TIME INCREMENT ARE
C READ. THE TIME STEP WILL CYCLE FROM THIS POINT.
C
	T=0
	TIME=0.D0
   40   T=T+1
	READ(252,50) DT,DUMFVC,DUMHCTC,(POX0(I),I=1,SLICE+1)
	READ(252,51) DUMVC,VS(1,T),VS(CHOP+1,T),
     1  DUMPC,PS(1,T),PS(CHOP+1,T),
     1  DUMCC,CS(16,1,T),CS(16,CHOP+1,T),
     1  (DUMCC,CS(I,1,T),CS(I,CHOP+1,T),I=1,11),
     1  (DUMCC,CS(I,1,T),CS(I,CHOP+1,T),I=13,15)
   50   FORMAT (2D12.4,F8.4,/,9F8.4)
   51   FORMAT (3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
	DO 45 IS=1,SLICE
	ICS=CHOP/SLICE
	DO 44 J=1,ICS+1
	RJ=DBLE(FLOAT(J-1)/FLOAT(ICS))
	POX((IS-1)*ICS+J,T)=(1.D0-RJ)*POX0(IS) + RJ*POX0(IS+1)
   44   CONTINUE
   45   CONTINUE
C
	DO 54 J=1,CHOP+1
	RJ=DBLE(FLOAT(J-1)/FLOAT(CHOP))
	VS(J,T)=(1.D0-RJ)*VS(1,T) + RJ*VS(CHOP+1,T)
	PS(J,T)=(1.D0-RJ)*PS(1,T) + RJ*PS(CHOP+1,T)
	CS(16,J,T)=(1.D0-RJ)*CS(16,1,T) + RJ*CS(16,CHOP+1,T)
	DO 52 I=1,SOLS
   52   CS(I,J,T)=(1.-RJ)*CS(I,1,T) + RJ*CS(I,CHOP+1,T)
	LCHS(J)=PKC + DLOG10(CS(4,J,T)/CS(5,J,T))
	CS(12,J,T)=10.**(-LCHS(J))
C
	TPS=CS(7,J,T)+CS(8,J,T)
	TNS=CS(10,J,T)+CS(11,J,T)
	TFS=CS(13,J,T)+CS(14,J,T)
C
	QPS=10.**(LCHS(J)-PKP)
	QNS=10.**(LCHS(J)-PKN)
	QFS=10.**(LCHS(J)-PKF)
C
	CS(8,J,T)=TPS/(1.+QPS)
	CS(7,J,T)=TPS-CS(8,J,T)
	CS(11,J,T)=TNS/(1.+QNS)
	CS(10,J,T)=TNS-CS(11,J,T)
	CS(14,J,T)=TFS/(1.+QFS)
	CS(13,J,T)=TFS-CS(14,J,T)
C
	CS(3,J,T)= CS(1,J,T)*Z(1) + CS(2,J,T)*Z(2)
	DO 53 I=4,SOLS
   53   CS(3,J,T)=CS(3,J,T)+CS(I,J,T)*Z(I)
   54   CONTINUE
C
C
	CCN(1,T)=.0081D0*HCTC(1,T)/0.45D0
	CCX(1,T)=.0307D0*HCTC(1,T)/0.45D0
	CCY(1,T)=.0453D0*HCTC(1,T)/0.45D0
	FVC(1,T)=FVC(1,T)*SCALE
	FBC(1,T)=FVC(1,T)/(1.D0-HCTC(1,T))
	TIME=TIME+FLOAT(T-1)*DT
	RTAU=2.D0*DBLE(FLOAT(TAU))/DT
C
	IF (COUNT.GT.1) GO TO 57
C
C  SET INITIAL GUESSES FOR THE HGB CONCENTRATIONS
C
	ESAT=HILL*(DLOG10(POX(1,T)/P50))
	ESAT=10**ESAT
	SAT(1,T)=ESAT/(1+ESAT)
C
	CC(17,1,T)=SAT(1,T)*CCN(1,T)*.3
	CC(18,1,T)=(1.D0-SAT(1,T))*CCN(1,T)*.3
	CC(19,1,T)=SAT(1,T)*CCN(1,T)*.3
	CC(20,1,T)=(1.D0-SAT(1,T))*CCN(1,T)*.3
	CC(21,1,T)=SAT(1,T)*CCN(1,T)*.3
	CC(22,1,T)=(1.D0-SAT(1,T))*CCN(1,T)*.3
	CC(23,1,T)=SAT(1,T)*CCX(1,T)*.5
	CC(24,1,T)=(1.D0-SAT(1,T))*CCX(1,T)*.5
	CC(25,1,T)=SAT(1,T)*CCX(1,T)*.5
	CC(26,1,T)=(1.D0-SAT(1,T))*CCX(1,T)*.5
	CC(27,1,T)=SAT(1,T)*CCY(1,T)*.5
	CC(28,1,T)=(1.D0-SAT(1,T))*CCY(1,T)*.5
	CC(29,1,T)=SAT(1,T)*CCY(1,T)*.5
	CC(30,1,T)=(1.D0-SAT(1,T))*CCY(1,T)*.5
C
C USE IVR TO CALCULATE CHOP AND CL (NO CHANGE IN DX):
C   IVR = 1,2 FOR OM          IVR = 3,7 FOR IM           IVR = 8 FOR MR
C
	IF (IVR.LE.2) THEN
	CHOP = IVR*(CHOP/2)
	CL = FLOAT(IVR)*CL/2.
	ELSE IF (IVR.LE.7) THEN
	CHOP = (IVR-2)*(CHOP/5)
	CL = FLOAT(IVR-2)*CL/5.
	ENDIF
C
C   START THE SPATIAL ITERATION.
C
   57   X=0
	CCN(1,T)=.0081D0*HCTC(1,T)/0.45D0
	CCX(1,T)=.0307D0*HCTC(1,T)/0.45D0
	CCY(1,T)=.0453D0*HCTC(1,T)/0.45D0
	FBC(1,T)=FVC(1,T)/(1.D0-HCTC(1,T))
	DO 400 KX=1,CHOP+1
	DIST= DIST + CL*FLOAT(KX-1)/FLOAT(CHOP)
C	PRINT 55, TIME,DIST,KX-1
   55   FORMAT(1X,'TIME=',F7.3,4X,'DIST=',F7.4,4X,'CHOP=',I2)
        IF (KX.EQ.1) THEN
        NN=15
        ELSE
	NN=15 + (SOLS+2)
        ENDIF
   60   X=X+1
   70   ITER=0
	PREMAX=0.D0
   80   ITER=ITER+1
C	PRINT 82, ITER
   82   FORMAT(1X,45X,'ITER=',I5,$)
	IF (ITER.GT.1) PREMAX=MAXPHI
	CALL MRAFLUX
	CALL MRAERR(PHI)
	MAXPHI=DABS(PHI(1))
	IMXPHI=1
	DO 230 J=2,NN
	IF (MAXPHI.GT.DABS(PHI(J))) GO TO 230
	MAXPHI=DABS(PHI(J))
	IMXPHI=J
  230   CONTINUE
C	PRINT 232,MAXPHI,IMXPHI
  232   FORMAT ('+',D12.4,I5)
	IF(MAXPHI.LT.EPSI) GO TO 558
C
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
C
	CALL MRAGAM(2,GAMMA)
	DO  900 K=1,NN
	DGAM=GAMMA(K)*1.D-2
	IF (DGAM.EQ.0.D0) DGAM=1.D-4
	GAMMA(K)=GAMMA(K)+DGAM
	CALL MRAGAM(1,GAMMA)
	CALL MRAFLUX
	CALL MRAERR(PPHI)
	DO 825 J=1,NN
  825   NDERIV(J,K)=(PPHI(J)-PHI(J))/DGAM
	GAMMA(K)=GAMMA(K)-DGAM
  900   CONTINUE
C
C
	CALL LES(NDERIV,PHI,PR,PFL,NN,SW,PTOL,DELGAM,M1)
	IF (SW) 930,600,930
C
  930   DO 940 J=1,NN
  940   GAMMA(J)=GAMMA(J)-DELGAM(J)
C
C
	DO 970 J=1,14
  970   IF (GAMMA(J).LT.0.D0) GAMMA(J)=1.D-10
	DO 975 J=18,32
  975   IF (GAMMA(J).LT.0.D0) GAMMA(J)=1.D-10
C
C RESET THE GUESSES IN THE NAMES BY WHICH WE KNOW THEM.
C
	CALL MRAGAM(1,GAMMA)
	IF (ITER.LT.35) GO TO 80
C IF THE SPATIAL ITERATION FAILS TO CONVERGE WITHIN 25 PASSES
C THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS')
	GO TO 600
C
C
C THE SPATIAL STEP CONVERGED,AND
C THE PROGRAM PICKS UP AT THIS POINT:
C
  558   IF (KX.EQ.CHOP+1) GO TO 400
C
C IF THIS IS NOT THE END OF THE CAPILLARY THE VARIABLES NEED RESETTING
	CALL MRARESET(1)
C
  400   CONTINUE
C
C IF THERE IS NO NEED FOR OUTPUT GO TO THE END
        IF ((SWOUT.EQ.0).OR.(SWOUT.EQ.2)) GO TO 600
C
C  IF THIS IS THE LAST SPATIAL CHOP, STOP
C
	IF (T.EQ.1) CALL MRARESA
	CALL MRARESB
	IF (TAU.EQ.1) GO TO 570
C
C  THUS IN THE STEADY STATE CASE WE ARE DONE.
	PRINT 560, COUNT
  560   FORMAT (' PROBLEM ',I2,' SOLVED FOR THE MRAVR')
	WRITE(254,565) DT,FVC(CHOP+1,T),HCTC(CHOP+1,T),SAT(CHOP+1,T),
     1  VC(CHOP+1,T), PC(CHOP+1,T), CC(16,CHOP+1,T),
     1  (CC(I,CHOP+1,T),I=1,11),
     1  (CC(I,CHOP+1,T),I=13,15)
  565   FORMAT (2D12.4,2F8.4,/,F8.4,/,F8.4,/,F8.4,/,(F14.9))
	GO TO 600
C
C  FOR THE TRANSIENT EXPERIMENT,THE VARIABLES AT T=1 ARE SET;
C  IF NEEDED, THE GUESS FOR T=2 IS GENERATED.
C
  570   CALL MRARESET(2)
C
C  IF TIME IS UP,STOP.
  585   TLIM=TLIM-1
	IF (TLIM.GT.0) GO TO 40
	PRINT 590
  590   FORMAT (' TIME IS UP')
C
C
  600   RETURN
C
	END
