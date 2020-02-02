	SUBROUTINE MRDRESET(SWR)
C
C PROGRAM TO RESET VARIABLES EITHER AT THE END OF A SPATIAL INTEGRATION
C OR AT THE BEGINNING OF A NEW TIME STEP
C
C
	INTEGER SOLS,T,X,CHOP,SWR
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
C
	COMMON/PARVR/ SOLS,T,X,CHOP,
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
	COMMON/VARVR/
     1   VC,PC,CC,VS,PS,CS,IMPC,IMPS,
     1   LCHC,LCHS,XC,XS,CCS,
     1   SAT,CCN,CCX,CCY,HCTC,FBC,POX,
     1   SC,AC,FVC,FKC,QC,QV,FVCS,FKCS,
     1   FVCSE,FKCSE,FVCSI,FKCSI
C
C
	GO TO (401,501) SWR
C
C IF THIS IS NOT THE END OF THE TUBULE THE VARIABLES NEED RESETTING
C
  401   VC(X+1,T)=VC(X,T)
	PC(X+1,T)=PC(X,T)
	DO 413 I=1,SOLS+15
  413   CC(I,X+1,T)=CC(I,X,T)
	FVC(X+1,T)=FVC(X,T)
	RETURN
C
C
C  FOR THE TRANSIENT EXPERIMENT
C  IF T=1, THE GUESS FOR T=2 IS GENERATED.
C  IF T=2, THE VALUES AT T=2 ARE MOVED TO T=1.
C
  501   GO TO (561,571) T
  561   VC(X,2)=VC(X,1)
	PC(X,2)=PC(X,1)
	DO 563 I=1,SOLS
  563   CC(I,X,2)=CC(I,X,1)
	FVC(X,2)=FVC(X,1)
	RETURN
C
  571   DO 590 X=1,CHOP+1
	VC(X,1)=VC(X,2)
	PC(X,1)=PC(X,2)
	DO 583 I=1,SOLS
  583   CC(I,X,1)=CC(I,X,2)
	FVC(X,1)=FVC(X,2)
  590   CONTINUE
	RETURN
C
	END
