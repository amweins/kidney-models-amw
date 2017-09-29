	SUBROUTINE OIAERR(PHI)
C
	INTEGER SOLS,T,X,CHOP
C
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
	DOUBLE PRECISION PHI
	DIMENSION PHI(40)
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
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C  AT THE INLET (X=0) THERE ARE 15 EQUATIONS TO BE SATISFIED
C
C THE OPEN CIRCUIT CONDITION
C
	PHI(15)=0.D0
	DO 115 I=1,SOLS
  115   PHI(15)=PHI(15)+F*Z(I)*FKCS(I,X,T)
C
C THE EIGHT BUFFER RELATIONS
C
C	17<-11      ONH2        KZO
C       18<-12      NH2         KZD
C	19<-13      ONH3+
C	20<-14      NH3+ 
C	21<-15      ONHCO2-     KCO
C	22<-16      NHCO2-      KCD
C	23<-17      OXH 
C	24<-18      XH   
C	25<-19      OX-         KXO
C	26<-20      X-          KXD
C	27<-21      OYH  
C	28<-22      YH   
C	29<-23      OY-         KYO
C	30<-24      Y-          KYD
C
	PHI(1)=PKC+
     1   DLOG10((KXO*CC(23,X,T)*CC(4,X,T))/(CC(25,X,T)*CC(5,X,T)))
	PHI(2)=PKC+
     1   DLOG10((KXD*CC(24,X,T)*CC(4,X,T))/(CC(26,X,T)*CC(5,X,T)))
	PHI(3)=PKC+
     1   DLOG10((KYO*CC(27,X,T)*CC(4,X,T))/(CC(29,X,T)*CC(5,X,T)))
	PHI(4)=PKC+
     1   DLOG10((KYD*CC(28,X,T)*CC(4,X,T))/(CC(30,X,T)*CC(5,X,T)))
	PHI(5)=PKC+
     1   DLOG10((KZO*CC(19,X,T)*CC(4,X,T))/(CC(17,X,T)*CC(5,X,T)))
	PHI(6)=PKC+
     1   DLOG10((KZD*CC(20,X,T)*CC(4,X,T))/(CC(18,X,T)*CC(5,X,T)))
	PHI(7)=PKC+ DLOG10
     1   ((KCO*CC(17,X,T)*CC(6,X,T)*CC(4,X,T))/(CC(21,X,T)*CC(5,X,T)))
	PHI(8)=PKC+ DLOG10
     1   ((KCD*CC(18,X,T)*CC(6,X,T)*CC(4,X,T))/(CC(22,X,T)*CC(5,X,T)))
C
C FOR OTHER THAN THE INITIAL STEP CHOOSE A DIFFERENCE SCHEME
C
	IF (X.GT.1) GO TO 300
C
C THE SIX INITIAL HGB SPECIES CONCENTRATIONS
C
	PHI(9)=CCN(X,T)- (CC(17,X,T)+CC(18,X,T)+CC(19,X,T)
     1  +CC(20,X,T)+CC(21,X,T)+CC(22,X,T))
	PHI(10)=CCX(X,T)- (CC(23,X,T)+CC(24,X,T)+CC(25,X,T)+CC(26,X,T))
	PHI(11)=CCY(X,T)- (CC(27,X,T)+CC(28,X,T)+CC(29,X,T)+CC(30,X,T))
	PHI(12)=SAT(X,T)*CCN(X,T)- (CC(17,X,T)+CC(19,X,T)+CC(21,X,T))
	PHI(13)=SAT(X,T)*CCX(X,T)- (CC(23,X,T)+CC(25,X,T))
	PHI(14)=SAT(X,T)*CCY(X,T)- (CC(27,X,T)+CC(29,X,T))
C
	RETURN
C
C FOR OTHER THAN THE INITIAL SPATIAL STEP (X=2) CONTINUE HERE
C   HERE IS A SCHEME CENTERED IN SPACE AND TIME
C
  200   IF (T.EQ.1) GO TO 220
C
C   MASS GENERATION IN THE TIME-DEPENDENT CASE
C
C
	QV=(FVC(X,T)-FVC(X-1,T))/DX
     1  +SC(X,T)*FVCS(X,T)
     1  +SC(X-1,T)*FVCS(X-1,T)
     1            +(FVC(X,T-1)-FVC(X-1,T-1))/DX
     1  +SC(X,T-1)*FVCS(X,T-1)
     1  +SC(X-1,T-1)*FVCS(X-1,T-1)
     1       +RTAU*(AC(X,T)-AC(X,T-1))
     1       +RTAU*(AC(X-1,T)-AC(X-1,T-1))
	DO 210 I=1,SOLS
	QC(I)=(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1  +SC(X,T)*FKCS(I,X,T)
     1  +SC(X-1,T)*FKCS(I,X-1,T)
     1            +(FKC(I,X,T-1)-FKC(I,X-1,T-1))/DX
     1  +SC(X,T-1)*FKCS(I,X,T-1)
     1  +SC(X-1,T-1)*FKCS(I,X-1,T-1)
     1   +RTAU*(AC(X,T)*CC(I,X,T)-AC(X,T-1)*CC(I,X,T-1))
     1   +RTAU*(AC(X-1,T)*CC(I,X-1,T)-AC(X-1,T-1)*CC(I,X-1,T-1))
  210   CONTINUE
	DO 215 I=SOLS+1,SOLS+15
	QC(I)=(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1            +(FKC(I,X,T-1)-FKC(I,X-1,T-1))/DX
     1   +RTAU*(AC(X,T)*CC(I,X,T)-AC(X,T-1)*CC(I,X,T-1))
     1   +RTAU*(AC(X-1,T)*CC(I,X-1,T)-AC(X-1,T-1)*CC(I,X-1,T-1))
  215   CONTINUE
C
	GO TO 240
C
C   MASS GENERATION IN THE STEADY-STATE CASE
C
  220   QV=(FVC(X,T)-FVC(X-1,T))/DX
     1  +SC(X,T)*FVCS(X,T)
     1  +SC(X-1,T)*FVCS(X-1,T)
	DO 230 I=1,SOLS
	QC(I)=(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1  +SC(X,T)*FKCS(I,X,T)
     1  +SC(X-1,T)*FKCS(I,X-1,T)
  230   CONTINUE
	DO 235 I=SOLS+1,SOLS+15
	QC(I)=(FKC(I,X,T)-FKC(I,X-1,T))/DX
  235   CONTINUE
C
  240   GO TO 400
C
C FOR OTHER THAN THE INITIAL SPATIAL STEP (X=2) CONTINUE HERE
C   HERE IS A SCHEME BACKWARD IN SPACE AND CENTERED IN TIME
C
  300   IF (T.EQ.1) GO TO 320
C
C   MASS GENERATION IN THE TIME-DEPENDENT CASE
C
C
	QV=.5*(FVC(X,T)-FVC(X-1,T))/DX
     1  +SC(X,T)*FVCS(X,T)
     1            +.5*(FVC(X,T-1)-FVC(X-1,T-1))/DX
     1  +SC(X,T-1)*FVCS(X,T-1)
     1       +RTAU*(AC(X,T)-AC(X,T-1))
	DO 310 I=1,SOLS
	QC(I)=.5*(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1  +SC(X,T)*FKCS(I,X,T)
     1            +.5*(FKC(I,X,T-1)-FKC(I,X-1,T-1))/DX
     1  +SC(X,T-1)*FKCS(I,X,T-1)
     1   +RTAU*(AC(X,T)*CC(I,X,T)-AC(X,T-1)*CC(I,X,T-1))
  310   CONTINUE
	DO 315 I=SOLS+1,SOLS+15
	QC(I)=.5*(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1            +.5*(FKC(I,X,T-1)-FKC(I,X-1,T-1))/DX
     1   +RTAU*(AC(X,T)*CC(I,X,T)-AC(X,T-1)*CC(I,X,T-1))
  315   CONTINUE
C
	GO TO 340
C
C   MASS GENERATION IN THE STEADY-STATE CASE
C
  320   QV=.5*(FVC(X,T)-FVC(X-1,T))/DX
     1  +SC(X,T)*FVCS(X,T)
	DO 330 I=1,SOLS
	QC(I)=.5*(FKC(I,X,T)-FKC(I,X-1,T))/DX
     1  +SC(X,T)*FKCS(I,X,T)
  330   CONTINUE
	DO 335 I=SOLS+1,SOLS+15
	QC(I)=.5*(FKC(I,X,T)-FKC(I,X-1,T))/DX
  335   CONTINUE
C
  340   CONTINUE
C
C
C HGB SPECIES CONSERVATION
C
  400   PHI(9)=QC(17)+QC(18)+QC(19)+QC(20)+QC(21)+QC(22)
	PHI(10)=QC(23)+QC(24)+QC(25)+QC(26)
	PHI(11)=QC(27)+QC(28)+QC(29)+QC(30)
C
	PHI(12)=SAT(X,T)*CCN(X,T)- (CC(17,X,T)+CC(19,X,T)+CC(21,X,T))
	PHI(13)=SAT(X,T)*CCX(X,T)- (CC(23,X,T)+CC(25,X,T))
	PHI(14)=SAT(X,T)*CCY(X,T)- (CC(27,X,T)+CC(29,X,T))
C
C   POISEUILLE FLOW
C
	PHI(16)=(PC(X,T)-PC(X-1,T))/DX
     1  +0.5*(25.1*ETA*FBC(X,T)/(AC(X,T)**2)
     1  +25.1*ETA*FBC(X-1,T)/(AC(X-1,T)**2))
C
C
C   MASS CONSERVATION
C
	PHI(17)=QV
	PHI(18)=QC(16)
	PHI(19)=QC(1)
	PHI(20)=QC(2)
	PHI(21)=QC(3)
	PHI(22)=QC(4)+QC(5)+QC(6)+QC(21)+QC(22)
	PHI(23)=QC(6)+QC(21)+QC(22)+
     1   AC(X,T)*(KHY*CC(6,X,T)-KDHY*CC(5,X,T))
	PHI(24)=QC(9)
	PHI(25)=QC(15)
	PHI(26)=QC(7)+QC(8)
	PHI(27)=QC(10)+QC(11)
	PHI(28)=QC(13)+QC(14)
	PHI(29)=QC(12)-QC(4)-QC(7)+QC(11)-QC(13)
     1  +QC(19)+QC(20)-QC(21)-QC(22)-QC(25)-QC(26)-QC(29)-QC(30)
C
C
C   PH EQUILIBRIA:
C
	PHI(30)=PKP-PKC
     1     -DLOG10((CC(4,X,T)*CC(8,X,T))/(CC(5,X,T)*CC(7,X,T)))
	PHI(31)=PKN-PKC
     1     -DLOG10((CC(4,X,T)*CC(11,X,T))/(CC(5,X,T)*CC(10,X,T)))
	PHI(32)=PKF-PKC
     1     -DLOG10((CC(4,X,T)*CC(14,X,T))/(CC(5,X,T)*CC(13,X,T)))
C
	RETURN
	END
C