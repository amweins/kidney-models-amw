	SUBROUTINE OMPCTRESET(SWR,CHOP)
C
C PROGRAM TO RESET VARIABLES EITHER AT THE END OF A SPATIAL INTEGRATION
C OR AT THE BEGINNING OF A NEW TIME STEP
C
	INTEGER SOLS,T,X,CHOP,SWR
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
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12,81),ATIS(3,12,81),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
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
     1   LMI,LIS,ATMI,ATIS,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
C
	GO TO (401,501) SWR
C
C IF NOT THE END OF THE TUBULE A NEW GUESS IS GENERATED FOR SPATIAL VARIABLES
C
  401   L(X+1,T)=L(X,T)
	AE(X+1,T)=AE(X,T)
	CHVL(X+1,T)=CHVL(X,T)
	VE(X+1,T)=VE(X,T)
	PE(X+1,T)=PE(X,T)
	FEVM(X+1,T)=FEVM(X,T)
	FEVS(X+1,T)=FEVS(X,T)
	DO 405 I=1,SOLS
	CE(I,X+1,T)=CE(I,X,T)
	FEKM(I,X+1,T)=FEKM(I,X,T)
  405   FEKS(I,X+1,T)=FEKS(I,X,T)
	DO 409 K=1,2
	VI(K,X+1,T)=VI(K,X,T)
	PI(K,X+1,T)=PI(K,X,T)
	IMP(K,X+1)=IMP(K,X)
	FIVM(K,X+1,T)=FIVM(K,X,T)
	FIVS(K,X+1,T)=FIVS(K,X,T)
	JV(K,X+1,T)=JV(K,X,T)
	CLVL(K,X+1,T)=CLVL(K,X,T)
	CBUF(K,X+1,T)=CBUF(K,X,T)
	HCBUF(K,X+1,T)=HCBUF(K,X,T)
	DO 409 I=1,SOLS
	CI(K,I,X+1,T)=CI(K,I,X,T)
	FIKM(K,I,X+1,T)=FIKM(K,I,X,T)
	FIKS(K,I,X+1,T)=FIKS(K,I,X,T)
	JK(K,I,X+1,T)=JK(K,I,X,T)
  409   CONTINUE
	VM(X+1,T)=VM(X,T)
	PM(X+1,T)=PM(X,T)
	DO 413 I=1,SOLS
  413   CM(I,X+1,T)=CM(I,X,T)
	IMPM(X+1,T)=IMPM(X,T)
	FVM(X+1,T)=FVM(X,T)
	RETURN
C
C
C  FOR THE TRANSIENT EXPERIMENT
C  IF T=1, THE GUESS FOR T=2 IS GENERATED.
C  IF T=2, THE VALUES AT T=2 ARE MOVED TO T=1.
C
  501   GO TO (561,571) T
  561   X=X-CHOP
	L(X,2)=L(X,1)
	AE(X,2)=AE(X,1)
	CHVL(X,2)=CHVL(X,1)
	VE(X,2)=VE(X,1)
	PE(X,2)=PE(X,1)
	FEVM(X,2)=FEVM(X,1)
	FEVS(X,2)=FEVS(X,1)
	DO 565 I=1,SOLS
	CE(I,X,2)=CE(I,X,1)
	FEKM(I,X,2)=FEKM(I,X,1)
  565   FEKS(I,X,2)=FEKS(I,X,1)
	DO 569 K=1,2
	VI(K,X,2)=VI(K,X,1)
	PI(K,X,2)=PI(K,X,1)
	FIVM(K,X,2)=FIVM(K,X,1)
	FIVS(K,X,2)=FIVS(K,X,1)
	JV(K,X,2)=JV(K,X,1)
	CLVL(K,X,2)=CLVL(K,X,1)
	CBUF(K,X,2)=CBUF(K,X,1)
	HCBUF(K,X,2)=HCBUF(K,X,1)
	DO 569 I=1,SOLS
	CI(K,I,X,2)=CI(K,I,X,1)
	FIKM(K,I,X,2)=FIKM(K,I,X,1)
	FIKS(K,I,X,2)=FIKS(K,I,X,1)
	JK(K,I,X,2)=JK(K,I,X,1)
  569   CONTINUE
	VM(X,2)=VM(X,1)
	X=X+CHOP
	RETURN
C
  571   DO 590 X=1,CHOP+1
	L(X,1)=L(X,2)
	AE(X,1)=AE(X,2)
	CHVL(X,1)=CHVL(X,2)
	VE(X,1)=VE(X,2)
	PE(X,1)=PE(X,2)
	FEVM(X,1)=FEVM(X,2)
	FEVS(X,1)=FEVS(X,2)
	DO 575 I=1,SOLS
	CE(I,X,1)=CE(I,X,2)
	FEKM(I,X,1)=FEKM(I,X,2)
  575   FEKS(I,X,1)=FEKS(I,X,2)
	DO 579 K=1,2
	VI(K,X,1)=VI(K,X,2)
	PI(K,X,1)=PI(K,X,2)
	FIVM(K,X,1)=FIVM(K,X,2)
	FIVS(K,X,1)=FIVS(K,X,2)
	JV(K,X,1)=JV(K,X,2)
	CLVL(K,X,1)=CLVL(K,X,2)
	CBUF(K,X,1)=CBUF(K,X,2)
	HCBUF(K,X,1)=HCBUF(K,X,2)
	DO 579 I=1,SOLS
	CI(K,I,X,1)=CI(K,I,X,2)
	FIKM(K,I,X,1)=FIKM(K,I,X,2)
	FIKS(K,I,X,1)=FIKS(K,I,X,2)
	JK(K,I,X,1)=JK(K,I,X,2)
  579   CONTINUE
	VM(X,1)=VM(X,2)
	PM(X,1)=PM(X,2)
	DO 583 I=1,SOLS
  583   CM(I,X,1)=CM(I,X,2)
	IMPM(X,1)=IMPM(X,2)
	FVM(X,1)=FVM(X,2)
  590   CONTINUE
	RETURN
C
	END
