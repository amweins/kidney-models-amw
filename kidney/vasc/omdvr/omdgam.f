	SUBROUTINE OMDGAM(SWG,GAMMA)
C PROGRAM TO SET AND RESET THE THE GAMMA ARRAY OF NN INDEPENDENT
C VARIABLES.  FOR SWG=1 THE GAMMA ARE DEFINED; FOR SWG=2, THE
C NAMED VARIABLES ARE RESET.
C
	INTEGER SOLS,T,SWG,X,CHOP
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
	DOUBLE PRECISION GAMMA
	DIMENSION GAMMA(40)
C
C
	GO TO (100,200) SWG
C
C
  100   DO 34 I=1,14
   34   CC(SOLS+1+I,X,T)=GAMMA(I)
	VC(X,T)=GAMMA(15)
	PC(X,T)=GAMMA(16)
	FVC(X,T)=GAMMA(17)
	CC(SOLS+1,X,T)=GAMMA(18)
	DO 36 I=1,11
   36   CC(I,X,T)=GAMMA(18+I)
	DO 37 I=12,14
   37   CC(I+1,X,T)=GAMMA(18+I)
C SCALING
        CC(5,X,T)=GAMMA(18+5)/1.D3
        CC(10,X,T)=GAMMA(18+10)/1.D3
        CC(15,X,T)=GAMMA(18+14)/1.D3
	RETURN
C
C
  200   DO 44 I=1,14
   44   GAMMA(I)=CC(SOLS+1+I,X,T)
	GAMMA(15)=VC(X,T)
	GAMMA(16)=PC(X,T)
	GAMMA(17)=FVC(X,T)
	GAMMA(18)=CC(SOLS+1,X,T)
	DO 46 I=1,11
   46   GAMMA(18+I)=CC(I,X,T)
	DO 47 I=12,14
   47   GAMMA(18+I)=CC(I+1,X,T)
C SCALING
        GAMMA(18+5)=CC(5,X,T)*1.D3
        GAMMA(18+10)=CC(10,X,T)*1.D3
        GAMMA(18+14)=CC(15,X,T)*1.D3
	RETURN
C
	END
