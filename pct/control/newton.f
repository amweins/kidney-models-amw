	SUBROUTINE NEWTON(NN,GAMMA)
C
	INTEGER SOLS,TAU,T,TLIM,COUNT,EXP,SW1,NPAR
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,L0,L(2),
     1   PKC,PKF,PKN,PKP,KHY(4),KDHY(4)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(15),IMPM,LCHM,XM(15),
     1   VS,PS,CS(15),IMPS,LCHS,XS(15)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(2),PE(2),CE(15,2),LCHE,XE(15),
     1   FEVM(2),FEKM(15,2),FEVS(2),FEKS(15,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),ATMI(3,15),ATIS(3,15),
     1   VI(3,2),PI(3,2),CI(3,15,2),IMP(3),LCHI(3),XI(3,15),
     1   FIVM(3,2),FIKM(3,15,2),FIVS(3,2),FIKS(3,15,2),CURI(3),
     1   JV(3,2),JK(3,15,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3),JHK(3),JHP(3),JAE1(3),JTSC(3),JNHE3(3,3),
     1   QIAMM
C COTRANSPORTER INFORMATION
	INTEGER  
     1   NMI,ORDMI(9),IMI(9,9),JMI(9,9),MMI(9,9),
     1   NIS,ORDIS(9),IIS(9,9),JIS(9,9),MIS(9,9)
	DOUBLE PRECISION RMI(9),RIS(9)
C CONTROL VARIABLES
	DOUBLE PRECISION 
     1   HCON(10),HCON0(10),DHCON(10,4),CON(4),CON0(10,4)
C
C
	DOUBLE PRECISION PHI,GAMMA,DGAM,PPHI,MAXPHI,PREMAX,
     1    DELGAM,NDERIV,RPHI
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
	DIMENSION PHI(60),GAMMA(60),PPHI(60),NDERIV(60,60),
     1    RPHI(60),DELGAM(60)
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(60),PFL(60)
C
C
        COMMON SOLS,T,TAU,TLIM,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
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
	COMMON/COTR/
     1   NMI,ORDMI,IMI,JMI,MMI,
     1   NIS,ORDIS,IIS,JIS,MIS,RMI,RIS,
     1   HCON,HCON0,DHCON,CON,CON0
C
C
	DGAM=1.D-10
	PTOL=1.D-8
	M1=60
	IERR=0
C
   70   ITER=0
   80   ITER=ITER+1
	PRINT 82, ITER
   82   FORMAT(I5,$)
	IF (ITER.GT.1) PREMAX=MAXPHI
	CALL ERRVEC(PHI)
	MAXPHI=DABS(PHI(1))
	MAXJ=1
	DO 230 J=2,NN
	IF (MAXPHI.LT.DABS(PHI(J))) MAXJ=J
  230   IF (MAXPHI.LT.DABS(PHI(J))) MAXPHI=DABS(PHI(J))
	PRINT 232, maxj, MAXPHI
  232   FORMAT ('+',I5,D12.4)
	IF(MAXPHI.LT.EPSI) RETURN
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
	CALL GAMSET(2,GAMMA)
C
	DO  900 K=1,NN
	GAMMA(K)=GAMMA(K)+DGAM
C
	CALL GAMSET(1,GAMMA)
C
	CALL ERRVEC(PPHI)
C
	DO 825 J=1,NN
  825   NDERIV(J,K)=(PPHI(J)-PHI(J))/DGAM
C  825   NDERIV(J,K)=SNGL((PPHI(J)-PHI(J))/DGAM)
C
	GAMMA(K)=GAMMA(K)-DGAM
  900   CONTINUE
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C
  923   DO 924 J=1,NN
  924   RPHI(J)=PHI(J)
C
	CALL LES8(NDERIV,RPHI,PR,PFL,NN,SW,PTOL,DELGAM,M1)
	IF (SW) 930,990,930
C
  930   DO 940 J=1,NN
  940   GAMMA(J)=GAMMA(J)-DELGAM(J)
C
C
C RESET THE GUESSES IN THE NAMES BY WHICH WE KNOW THEM.
C
	CALL GAMSET(1,GAMMA)
C
	IF (ITER.LT.45) GO TO 80
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE WITHIN 15 PASSES
C THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS')
	STOP
  990   PRINT 995
  995   FORMAT (' SMALL PIVOT DIAGNOSTIC')
	STOP
	END
