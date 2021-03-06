	SUBROUTINE OMAFLUX
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
     1   FVCSE(901,2),FKCSE(16,901,2),FVCSI(901,2),FKCSI(16,901,2),
     1   QV,QC(31),FVCS(901,2),FKCS(16,901,2),
     1   HCTC(901,2),FBC(901,2),POX(901,2)
C
        DOUBLE PRECISION ZCS(15),ESAT
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
C
	SC(X,T)=6.28D0*RC0*(1.D0+MUC*(PC(X,T)-PS(X,T)))*SCALE
	AC(X,T)=((SC(X,T)**2)/12.56D0)/SCALE
	IMPC(X,T)=A0*CC(SOLS+1,X,T) + A1*CC(SOLS+1,X,T)**2
     1     +A2*CC(SOLS+1,X,T)**3
	IMPS(X,T)=A0*CS(SOLS+1,X,T) + A1*CS(SOLS+1,X,T)**2
     1     +A2*CS(SOLS+1,X,T)**3
C
	LCHC(X)=PKC+DLOG10(CC(4,X,T)/CC(5,X,T))
	LCHS(X)=PKC+DLOG10(CS(4,X,T)/CS(5,X,T))
	CC(12,X,T)=(10.D0**(-PKC))*CC(5,X,T)/CC(4,X,T)
	CS(12,X,T)=(10.**(-PKC))*CS(5,X,T)/CS(4,X,T)
C
	ESAT=HILL*(DLOG10(POX(X,T)/P50)+BOHR*(7.40-LCHC(X)))
	ESAT=10**ESAT
	SAT(X,T)=ESAT/(1+ESAT)
C
	IF (X.GT.1) THEN
	CCN(X,T)= (CC(17,X,T)+CC(18,X,T)+CC(19,X,T)
     1     +CC(20,X,T)+CC(21,X,T)+CC(22,X,T))
	CCX(X,T)= (CC(23,X,T)+CC(24,X,T)+CC(25,X,T)+CC(26,X,T))
	CCY(X,T)= (CC(27,X,T)+CC(28,X,T)+CC(29,X,T)+CC(30,X,T))
	ELSE
	ENDIF
C
C COMPUTE THE FLUXES
C
	FVCSE(X,T)=LPCSE*(IMPS(X,T)-IMPC(X,T)+PC(X,T)-PS(X,T))/RT
	FVCSI(X,T)=LPCSI*(IMPS(X,T)-IMPC(X,T)+PC(X,T)-PS(X,T))/RT
	DO 95 I=1,SOLS
        FVCSE(X,T)=FVCSE(X,T)+LPCSE*SCSE(I)*(CS(I,X,T)-CC(I,X,T))
        FVCSI(X,T)=FVCSI(X,T)+LPCSI*SCSI(I)*(CS(I,X,T)-CC(I,X,T))
   95   CONTINUE
        FVCS(X,T) = FVCSE(X,T) + FVCSI(X,T)
C
	DO 100 I=1,SOLS
C
  354   CCS(I)=(CC(I,X,T)-CS(I,X,T))
	IF (CCS(I)) 356,355,356
  355   CCS(I)=CS(I,X,T)
	GO TO 360
  356   CCS(I)=CCS(I)/DLOG(CC(I,X,T)/CS(I,X,T))
C
  360   XC(I,X)=RTE*DLOG(CC(I,X,T))+Z(I)*F*VC(X,T)*1.D-6
	XS(I,X)=RTE*DLOG(CS(I,X,T))+Z(I)*F*VS(X,T)*1.D-6
C
C CONVECTIVE FLUX
	FKCSE(I,X,T)=FVCSE(X,T)*(1.D0-SCSE(I))*CCS(I)
	FKCSI(I,X,T)=FVCSI(X,T)*(1.D0-SCSI(I))*CCS(I)
C
C GOLDMAN FLUX
	ZCS(I)=Z(I)*F*(VC(X,T)-VS(X,T))*1.D-6/RTE
	IF (DABS(ZCS(I)).GT.0.1D-15) THEN
	FKCSE(I,X,T)=FKCSE(I,X,T)+HCSE(I)*ZCS(I)*
     1   (CC(I,X,T) - CS(I,X,T)*DEXP(-ZCS(I)))/(1.D0 - DEXP(-ZCS(I)))
	FKCSI(I,X,T)=FKCSI(I,X,T)+HCSI(I)*ZCS(I)*
     1   (CC(I,X,T) - CS(I,X,T)*DEXP(-ZCS(I)))/(1.D0 - DEXP(-ZCS(I)))
        ELSE
	FKCSE(I,X,T)=FKCSE(I,X,T)+HCSE(I)*(CC(I,X,T) - CS(I,X,T))
	FKCSI(I,X,T)=FKCSI(I,X,T)+HCSI(I)*(CC(I,X,T) - CS(I,X,T))
        ENDIF
        FKCS(I,X,T) = FKCSE(I,X,T) + FKCSI(I,X,T)
  100   CONTINUE
        FKCS(SOLS+1,X,T)=0.D0
C
C COMPUTE THE FLOWS
C
	DO 106 I=1,SOLS+1
  106   FKC(I,X,T)=FVC(X,T)*CC(I,X,T)
C
	IF (X.EQ.1) GO TO 108
	FBC(X,T)=FBC(X-1,T) + (FVC(X,T)-FVC(X-1,T))
	HCTC(X,T)=(HCTC(X-1,T)*FBC(X-1,T))/FBC(X,T)
C
  108   DO 110 I=SOLS+2,SOLS+15
  110   FKC(I,X,T)=FBC(X,T)*CC(I,X,T)
C
	RETURN
	END
