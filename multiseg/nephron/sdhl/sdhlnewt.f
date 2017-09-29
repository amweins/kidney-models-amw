	SUBROUTINE SDHLNEWT(COUNT,SWOUT)
C
	INTEGER SOLS,TAU,T,TLIM,COUNT,NPAR,X,CHOP,SWOUT
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
	DOUBLE PRECISION RJ,TPS,TNS,QPS,QNS
	DOUBLE PRECISION PHI,GAMMA,DGAM,PPHI,MAXPHI,
     1    DELGAM,NDERIV,RPHI
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       DGAM-   VARIABLE INCREMENT TO COMPUTE NUMERICAL DERIVATIVE
C       PPHI-   ERROR VECTOR AT GAMMA + DGAM
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C
	DIMENSION PHI(80),GAMMA(80),PPHI(80),NDERIV(80,80),
     1    RPHI(80),DELGAM(80)
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(80),PFL(80)
C
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
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
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
	DGAM=1.D-4
	PTOL=1.D-10
	M1=80
	IERR=0
C
C START THE DHL CALCULATION
C
	X=1
	IF (T.GT.1) GO TO 49
	NN=1
	REWIND 124
	READ (124,24) (GAMMA(J),J=1,NN)
   24   FORMAT (D30.20)
	CALL SDHLGAM(1,GAMMA)
	IF (SWOUT.EQ.2) GO TO 60
C
C PARAMETERS AND BOUNDARY VALUES ARE READ IN:
C IT IS ASSUMED THAT THE PERITUBULAR PROFILE IS A LINEAR FUNCTION OF X.
C ASSUMING EQUILIBRIUM BETWEEN DISSOLVED CO2 AND H2CO3
C
   49   CALL SDHLPSET
	READ(122,50) DUMDT,DUMFVM,
     1  DUMVM,VS(X,T),VS(CHOP+X,T),
     1  DUMPM,PS(X,T),PS(CHOP+X,T),
     1  DUMIMPM,IMPS(X,T),IMPS(CHOP+X,T),
     1  (DUMCM,CS(I,X,T),CS(I,CHOP+X,T),I=1,SOLS-1)
   50   FORMAT (2D12.4,/,3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
	DO 54 J=X,CHOP+X
	RJ=DBLE(FLOAT(J-X)/FLOAT(CHOP))
	VS(J,T)=(1.D0-RJ)*VS(X,T) + RJ*VS(CHOP+X,T)
	PS(J,T)=(1.D0-RJ)*PS(X,T) + RJ*PS(CHOP+X,T)
	IMPS(J,T)=(1.D0-RJ)*IMPS(X,T) + RJ*IMPS(CHOP+X,T)
	DO 52 I=1,SOLS-1
   52   CS(I,J,T)=(1.-RJ)*CS(I,X,T) + RJ*CS(I,CHOP+X,T)
	LCHS(J)=PKC + DLOG10(CS(4,J,T)/CS(5,J,T))
	CS(SOLS,J,T)=10.**(-LCHS(J))
C
	TPS=CS(7,J,T)+CS(8,J,T)
	TNS=CS(10,J,T)+CS(11,J,T)
C
	QPS=10.**(LCHS(J)-PKP)
	QNS=10.**(LCHS(J)-PKN)
C
	CS(8,J,T)=TPS/(1.+QPS)
	CS(7,J,T)=TPS-CS(8,J,T)
	CS(11,J,T)=TNS/(1.+QNS)
	CS(10,J,T)=TNS-CS(11,J,T)
C
	CS(3,J,T)=IMPS(J,T)*ZIMPS + CS(1,J,T)*Z(1) + CS(2,J,T)*Z(2)
	DO 53 I=4,SOLS
   53   CS(3,J,T)=CS(3,J,T)+CS(I,J,T)*Z(I)
   54   CONTINUE
C
	TIME=TIME+DT*FLOAT(T-1)
	RTAU=DBLE(FLOAT(TAU))/DT
	IF (TAU.EQ.1) CALL SDHLPTIM(NPAR)
C
C   START THE SPATIAL ITERATION.
C
   60   CONTINUE
	X=0
	DX=TL/DBLE(FLOAT(CHOP))
	DIST=0.0 - DX
C
	DO 400 KX=0,CHOP
	X=X+1
	DIST=DIST + DX
C	PRINT 55, COUNT,TIME,DIST,X-1
   55   FORMAT(1X,'EXP=',I3,4X,'TIME=',F7.3,4X,
     1   'DIST=',F7.4,4X,'CHOP=',I3)
	IF (KX.EQ.0) THEN
	NN=1
	ELSE 
	NN=(2+SOLS)+2
	ENDIF
   70   ITER=0
   80   ITER=ITER+1
C	PRINT 82, ITER
   82   FORMAT(46X,'ITER=',I5,$)
	CALL SDHLFLUX
	CALL SDHLERR(PHI)
	MAXPHI=DABS(PHI(1))
	DO 230 J=2,NN
  230   IF (MAXPHI.LT.DABS(PHI(J))) MAXPHI=DABS(PHI(J))
C	PRINT 232,MAXPHI
  232   FORMAT ('+',D12.4)
	IF(MAXPHI.LT.EPSI) GO TO 558
C
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
C
	CALL SDHLGAM(2,GAMMA)
C
	DO  900 K=1,NN
	GAMMA(K)=GAMMA(K)+DGAM
C
	CALL SDHLGAM(1,GAMMA)
C
	CALL SDHLFLUX
	CALL SDHLERR(PPHI)
C
	DO 825 J=1,NN
  825   NDERIV(J,K)=(PPHI(J)-PHI(J))/DGAM
C
	GAMMA(K)=GAMMA(K)-DGAM
  900   CONTINUE
C
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C
  923   DO 924 J=1,NN
  924   RPHI(J)=PHI(J)
C
	CALL LES8(NDERIV,RPHI,PR,PFL,NN,SW,PTOL,DELGAM,M1)
	IF (SW) 930,600,930
C
  930   DO 940 J=1,NN
  940   GAMMA(J)=GAMMA(J)-DELGAM(J)
C
C
C RESET THE GUESSES IN THE NAMES BY WHICH WE KNOW THEM.
C
	CALL SDHLGAM(1,GAMMA)
C
	IF (ITER.LT.50) GO TO 80
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS IN SDHL')
	WRITE (19,982) COUNT
  982   FORMAT ('SDHL ERROR IN EXP. ', I3)
	GO TO 600
C
C THE SPATIAL STEP CONVERGED,AND
C THE PROGRAM PICKS UP AT THIS POINT:
  558   CONTINUE
C
C IF THIS IS NOT THE END OF THE TUBULE THE VARIABLES NEED RESETTING
	IF (X.EQ.CHOP+1) GO TO 400
	CALL SDHLRESET(1)
  400   CONTINUE
C
	X=1
	CALL SDHLGAM(2,GAMMA)
	X=1+CHOP
	REWIND 125
	WRITE (125,624) (GAMMA(J),J=1,NN)
  624   FORMAT(D30.20)
C
C IF THERE IS NO NEED FOR OUTPUT GO TO THE END
	IF ((SWOUT.EQ.0).OR.(SWOUT.EQ.2)) GO TO 600
C
C  IF THIS IS THE LAST SPATIAL CHOP, STOP
C
	IF(T.EQ.2) GO TO 420
	DIST=DIST - TL*1.5
	X=X-3*(CHOP/2)
	DO 410 LL=1,3
	DIST=DIST + TL*0.5
	X=X+(CHOP/2)
	IF(LL.EQ.1) CALL SDHLRESA
	CALL SDHLRESB
  410   CALL SDHLRESC
  420   CALL SDHLRESD
	CALL SDHLRESF
C
	PRINT 560
  560   FORMAT (' PROBLEM SOLVED FOR THE DHL')
C
C  FOR THE TRANSIENT EXPERIMENT,THE VARIABLES AT T=1 ARE SET;
C  IF NEEDED, THE GUESS FOR T=2 IS GENERATED.
C
	IF (TAU.EQ.1) CALL SDHLRESET(2)
C
  600   RETURN
	END
