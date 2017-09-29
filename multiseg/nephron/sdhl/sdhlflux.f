	SUBROUTINE SDHLFLUX
C
	INTEGER SOLS,T,X,CHOP
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
	DOUBLE PRECISION VHK,SHK(4),HMI11,
     1   VNHE3(3),SNHE3(3,2),VTSC,STSC(4),
     1   ZME(12),ZMI(3,12),ZIE(3,12),ZES(12),ZIS(3,12),
     1   SNKCC(8),VNKCC(4),SKCC(6),VKCC(3)
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
C COMPUTE THE RELEVANT VOLUMES
C
	SM(X,T)=6.28D0*RM0*MUM
	AM(X,T)=3.14D0*(RM0**2)*MUM
C	SM(X,T)=6.28D0*RM0*(1.D0+MUM*(PM(X,T)-PS(X,T)))
C	AM(X,T)=((SM(X,T)**2)/12.56D0)
C
C COMPUTE THE MEMBRANE FLUXES
C
	FEVM(X,T)=AME*LPME*
     1   (PM(X,T)-PS(X,T)+RT*IMPS(X,T)-RT*IMPM(X,T))/RT
	DO 30 K=1,1
	FIVM(K,X,T)=AI0(K)*LPMI(K)*
     1   (PM(X,T)-PS(X,T)+RT*IMPS(X,T)-RT*IMPM(X,T))/RT
   30   JV(K,X,T)=LPIS(K)*AIE(K)*(PI(K,X,T)-PE(X,T)-RT*IMP(K,X))/RT
	DO 35 I=1,SOLS
	FEVM(X,T)=FEVM(X,T)+AME*LPME*SME(I)*(CS(I,X,T)-CM(I,X,T))
	DO 35 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)+AI0(K)*LPMI(K)*
     1   SMI(K,I)*(CS(I,X,T)-CM(I,X,T))
   35   CONTINUE
C
C
C THE CIVAN AND BOOKMAN MODIFICATION TO HMI(1,1)
	HMI11=HMI(1,1)
C	HMI(1,1)=HMI11*((1.+CM(1,X,T)/.03)**-1)*
C     1      (1.-CI(1,1,X,T)/.05)*1.666667
C
	DO 100 I=1,SOLS
C
	CME(I)=(CS(I,X,T)-CM(I,X,T))
	IF (CME(I)) 347,346,347
  346   CME(I)=CM(I,X,T)
	GO TO 351
  347   CME(I)=CME(I)/DLOG(CS(I,X,T)/CM(I,X,T))
  351   DO 360 K=1,1
	CMI(K,I)=(CS(I,X,T)-CM(I,X,T))
	IF (CMI(K,I)) 353,352,353
  352   CMI(K,I)=CM(I,X,T)
	GO TO 360
  353   CMI(K,I)=CMI(K,I)/DLOG(CS(I,X,T)/CM(I,X,T))
  360   CONTINUE
C
C
	XM(I,X)=RTE*DLOG(CM(I,X,T))+Z(I)*F*VM(X,T)*1.D-6
	XS(I,X)=RTE*DLOG(CS(I,X,T))+Z(I)*F*VS(X,T)*1.D-6
C
	ZME(I)=Z(I)*F*(VM(X,T)-VS(X,T))*1.D-6/RTE
	DO 380 K=1,1
  380   ZMI(K,I)=Z(I)*F*(VM(X,T)-VS(X,T))*1.D-6/RTE
C
C CONVECTIVE FLUXES
C
	FEKM(I,X,T)=FEVM(X,T)*(1.D0-SME(I))*CME(I)
	DO 390 K=1,1
  390   FIKM(K,I,X,T)=FIVM(K,X,T)*(1.D0-SMI(K,I))*CMI(K,I)
C
C GOLDMAN FLUXES
C
	IF (DABS(Z(I)).LT.0.1) GO TO 90
	IF (DABS(ZME(I)).GT.1.D-10) THEN
	FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*ZME(I)*
     1   (CM(I,X,T) - CS(I,X,T)*DEXP(-ZME(I)))/(1.D0 - DEXP(-ZME(I)))
	ELSE
	FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*(CM(I,X,T) - CS(I,X,T))
	ENDIF
	DO 85 K=1,1
	IF (DABS(ZMI(K,I)).GT.1.D-10) THEN
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*AI0(K)*ZMI(K,I)*
     1   (CM(I,X,T) - CS(I,X,T)*DEXP(-ZMI(K,I)))/
     1     (1.D0 - DEXP(-ZMI(K,I)))
	ELSE
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*AI0(K)*
     1   (CM(I,X,T) - CS(I,X,T))
	ENDIF
   85   CONTINUE
	GO TO 100
C
   90   FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*(CM(I,X,T) - CS(I,X,T))
	DO 95 K=1,1
   95   FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*
     1   AI0(K)*(CM(I,X,T) - CS(I,X,T))
  100   CONTINUE
C
C
C NET COTRANSPORTERS
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
  120   FIKM(K,I,X,T)=FIKM(K,I,X,T)+LMI(K,I,J)*AI0(K)*(XM(J,X)-XS(J,X))
C
	DO 110 I=1,SOLS
  110   FKM(I,X,T)=FVM(X,T)*CM(I,X,T)
	FKM(13,X,T)=FVM(X,T)*IMPM(X,T)
C
	RETURN
	END
