	SUBROUTINE OMCTFLUX
C
	INTEGER SOLS,T,X
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
	DOUBLE PRECISION VHK,SS(4),VAE1,
     1   ZME,ZMI,ZIE,ZES,ZIS
	DIMENSION
     1   ZME(12),ZMI(3,12),ZIE(3,12),ZES(12),ZIS(3,12)
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
C
C COMPUTE THE RELEVANT VOLUMES
	AE(X,T)=AE0*(1.0+MUA*(PE(X,T)-PM(X,T)))
	IF((AE(X,T).LT.AE0)) AE(X,T)=AE0
	CHVL(X,T)=CHVL0*(1.0+MUV*(PE(X,T)-PM(X,T)))
	L(X,T)=CHVL(X,T)
	DO 20 K=2,2
	CLVL(K,X,T)=CLVL0(K)
C	CLVL(K,X,T)=CLVL0(K)*IMP0(K)/IMP(K,X)
	L(X,T)=L(X,T)+CLVL(K,X,T)
	PI(K,X,T)=PM(X,T)
C
	HCBUF(K,X,T)=(TBUF(K)*CLVL0(K)/CLVL(K,X,T))/
     1   (1.0 + (10.**(-PKB(K)))/CI(K,12,X,T))
	CBUF(K,X,T)=(TBUF(K)*CLVL0(K)/CLVL(K,X,T))-HCBUF(K,X,T)
   20   CONTINUE
	SM(X,T)=6.28D0*RM0*MUM
	AM(X,T)=3.14D0*(RM0**2)*MUM
C	SM(X,T)=6.28D0*RM0*(1.D0+MUM*(PM(X,T)-PS(X,T)))
C	AM(X,T)=((SM(X,T)**2)/12.56D0)
C
C
C COMPUTE THE MEMBRANE FLUXES
C
	FEVM(X,T)=AME*LPME*(PM(X,T)-PE(X,T)-RT*IMPM(X,T))/RT
	FEVS(X,T)=AE(X,T)*LPES*(RT*IMPS(X,T)+PE(X,T)-PS(X,T))/RT
	DO 30 K=2,2
	FIVM(K,X,T)=AI0(K)*LPMI(K)*
     1   (RT*IMP(K,X)-RT*IMPM(X,T)+PM(X,T)-PI(K,X,T))/RT
	FIVS(K,X,T)=AI0(K)*LPIS(K)*
     1   (PI(K,X,T)-PS(X,T)+RT*IMPS(X,T)-RT*IMP(K,X))/RT
   30   JV(K,X,T)=LPIS(K)*AIE(K)*(PI(K,X,T)-PE(X,T)-RT*IMP(K,X))/RT
	DO 35 I=1,SOLS
	FEVM(X,T)=FEVM(X,T)+AME*LPME*SME(I)*(CE(I,X,T)-CM(I,X,T))
	FEVS(X,T)=FEVS(X,T)+AE(X,T)*LPES*SES(I)*(CS(I,X,T)-CE(I,X,T))
	DO 35 K=2,2
	FIVM(K,X,T)=FIVM(K,X,T)+AI0(K)*LPMI(K)*
     1   SMI(K,I)*(CI(K,I,X,T)-CM(I,X,T))
	FIVS(K,X,T)=FIVS(K,X,T)+AI0(K)*LPIS(K)*
     1   SIS(K,I)*(CS(I,X,T)-CI(K,I,X,T))
	JV(K,X,T)=JV(K,X,T)+LPIS(K)*AIE(K)*
     1   SIS(K,I)*(CE(I,X,T)-CI(K,I,X,T))
   35   CONTINUE
C
C
	DO 100 I=1,SOLS
C
	CME(I)=(CE(I,X,T)-CM(I,X,T))
	IF (CME(I)) 347,346,347
  346   CME(I)=CM(I,X,T)
	GO TO 348
  347   CME(I)=CME(I)/DLOG(CE(I,X,T)/CM(I,X,T))
  348   CES(I)=(CE(I,X,T)-CS(I,X,T))
	IF (CES(I)) 350,349,350
  349   CES(I)=CS(I,X,T)
	GO TO 351
  350   CES(I)=CES(I)/DLOG(CE(I,X,T)/CS(I,X,T))
  351   DO 360 K=2,2
	CMI(K,I)=(CI(K,I,X,T)-CM(I,X,T))
	IF (CMI(K,I)) 353,352,353
  352   CMI(K,I)=CM(I,X,T)
	GO TO 354
  353   CMI(K,I)=CMI(K,I)/DLOG(CI(K,I,X,T)/CM(I,X,T))
  354   CIS(K,I)=(CI(K,I,X,T)-CS(I,X,T))
	IF (CIS(K,I)) 356,355,356
  355   CIS(K,I)=CS(I,X,T)
	GO TO 357
  356   CIS(K,I)=CIS(K,I)/DLOG(CI(K,I,X,T)/CS(I,X,T))
  357   CIE(K,I)=(CI(K,I,X,T)-CE(I,X,T))
	IF (CIE(K,I)) 359,358,359
  358   CIE(K,I)=CI(K,I,X,T)
	GO TO 360
  359   CIE(K,I)=CIE(K,I)/DLOG(CI(K,I,X,T)/CE(I,X,T))
  360   CONTINUE
C
C
	XE(I,X)=RTE*DLOG(CE(I,X,T))+Z(I)*F*VE(X,T)*1.D-6
	XM(I,X)=RTE*DLOG(CM(I,X,T))+Z(I)*F*VM(X,T)*1.D-6
	XS(I,X)=RTE*DLOG(CS(I,X,T))+Z(I)*F*VS(X,T)*1.D-6
	DO 370 K=2,2
  370   XI(K,I,X)=RTE*DLOG(CI(K,I,X,T))+Z(I)*F*VI(K,X,T)*1.D-6
C
	ZME(I)=Z(I)*F*(VM(X,T)-VE(X,T))*1.D-6/RTE
	ZES(I)=Z(I)*F*(VE(X,T)-VS(X,T))*1.D-6/RTE
	DO 380 K=2,2
	ZMI(K,I)=Z(I)*F*(VM(X,T)-VI(K,X,T))*1.D-6/RTE
	ZIE(K,I)=Z(I)*F*(VI(K,X,T)-VE(X,T))*1.D-6/RTE
  380   ZIS(K,I)=Z(I)*F*(VI(K,X,T)-VS(X,T))*1.D-6/RTE
C
C CONVECTIVE FLUXES
C
	FEKM(I,X,T)=FEVM(X,T)*(1.D0-SME(I))*CME(I)
	FEKS(I,X,T)=FEVS(X,T)*(1.D0-SES(I))*CES(I)
	DO 390 K=2,2
	FIKM(K,I,X,T)=FIVM(K,X,T)*(1.D0-SMI(K,I))*CMI(K,I)
	FIKS(K,I,X,T)=FIVS(K,X,T)*(1.D0-SIS(K,I))*CIS(K,I)
  390   JK(K,I,X,T)=JV(K,X,T)*(1.D0-SIS(K,I))*CIE(K,I)
C
C GOLDMAN FLUXES
C
	IF (DABS(Z(I)).LT.0.1) GO TO 90
	FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*ZME(I)*
     1   (CM(I,X,T) - CE(I,X,T)*DEXP(-ZME(I)))/(1.D0 - DEXP(-ZME(I)))
	FEKS(I,X,T)=FEKS(I,X,T)+HES(I)*AE(X,T)*ZES(I)*
     1   (CE(I,X,T) - CS(I,X,T)*DEXP(-ZES(I)))/(1.D0 - DEXP(-ZES(I)))
	DO 85 K=2,2
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*AI0(K)*ZMI(K,I)*
     1   (CM(I,X,T) - CI(K,I,X,T)*DEXP(-ZMI(K,I)))/
     1     (1.D0 - DEXP(-ZMI(K,I)))
	JK(K,I,X,T)=JK(K,I,X,T)+HIS(K,I)*AIE(K)*ZIE(K,I)*
     1   (CI(K,I,X,T) - CE(I,X,T)*DEXP(-ZIE(K,I)))/
     1     (1.D0 - DEXP(-ZIE(K,I)))
   85   FIKS(K,I,X,T)=FIKS(K,I,X,T)+HIS(K,I)*AI0(K)*ZIS(K,I)*
     1   (CI(K,I,X,T) - CS(I,X,T)*DEXP(-ZIS(K,I)))/
     1     (1.D0 - DEXP(-ZIS(K,I)))
	GO TO 100
C
   90   FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*(CM(I,X,T) - CE(I,X,T))
	FEKS(I,X,T)=FEKS(I,X,T)+HES(I)*AE(X,T)*(CE(I,X,T) - CS(I,X,T))
	DO 95 K=2,2
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*
     1   AI0(K)*(CM(I,X,T) - CI(K,I,X,T))
	JK(K,I,X,T)=JK(K,I,X,T)+
     1   HIS(K,I)*AIE(K)*(CI(K,I,X,T) - CE(I,X,T))
   95   FIKS(K,I,X,T)=FIKS(K,I,X,T)+
     1   HIS(K,I)*AI0(K)*(CI(K,I,X,T) - CS(I,X,T))
  100   CONTINUE
C
C NET COTRANSPORTERS
C
	DO 120 K=2,2
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+LMI(K,I,J)*AI0(K)*(XM(J,X)-XI(K,J,X))
	FIKS(K,I,X,T)=FIKS(K,I,X,T)+LIS(K,I,J)*AI0(K)*(XI(K,J,X)-XS(J,X))
  120   JK(K,I,X,T)=JK(K,I,X,T)+LIS(K,I,J)*AIE(K)*(XI(K,J,X)-XE(J,X))
C
C
C AE1 cotransporter
C
	SS(1)=CI(2,3,X,T)
	SS(2)=CI(2,4,X,T)
	SS(3)=CE(3,X,T)
	SS(4)=CE(4,X,T)
	CALL AE1(VAE1,SS)
	JK(2,3,X,T)=JK(2,3,X,T)-NAE1(2)*AIE(2)*VAE1
	JK(2,4,X,T)=JK(2,4,X,T)+NAE1(2)*AIE(2)*VAE1
	SS(3)=CS(3,X,T)
	SS(4)=CS(4,X,T)
	CALL AE1(VAE1,SS)
	FIKS(2,3,X,T)=FIKS(2,3,X,T)-NAE1(2)*AI0(2)*VAE1
	FIKS(2,4,X,T)=FIKS(2,4,X,T)+NAE1(2)*AI0(2)*VAE1
C
C
C PUMPS
C
C  SODIUM PUMPS
C
	DO 125 K=2,2
        KNPN(K)=0.0002*(1.D0+CI(K,2,X,T)/.00833)
        KNPK(K)=0.0001*(1.D0+CS(1,X,T)/.0185)
        ATIS(K,1,X)=NP(K)*(CI(K,1,X,T)/(KNPN(K)+CI(K,1,X,T)))**3*
     1    (CS(2,X,T)/(KNPK(K)+CS(2,X,T)))**2
C  ALLOW FOR COMPETITION BETWEEN K+ AND NH4+ 
	ATIS(K,2,X)=-ATIS(K,1,X)*0.667*CE(2,X,T)/
     1    (CE(2,X,T)+CE(11,X,T)/KNH4(K))
	ATIS(K,11,X)=-ATIS(K,1,X)*0.667*CE(11,X,T)/
     1    (KNH4(K)*CE(2,X,T)+CE(11,X,T))
C
	JK(K,1,X,T)=JK(K,1,X,T)+AIE(K)*ATIS(K,1,X)
	JK(K,2,X,T)=JK(K,2,X,T)+AIE(K)*ATIS(K,2,X)
	JK(K,11,X,T)=JK(K,11,X,T)+AIE(K)*ATIS(K,11,X)
	FIKS(K,1,X,T)=FIKS(K,1,X,T)+AI0(K)*ATIS(K,1,X)
	FIKS(K,2,X,T)=FIKS(K,2,X,T)+AI0(K)*ATIS(K,2,X)
	FIKS(K,11,X,T)=FIKS(K,11,X,T)+AI0(K)*ATIS(K,11,X)
  125   CONTINUE
C
C  PROTON PUMPS
C
	ATMI(2,2,X)=0.D0
 	ATMI(2,12,X)=-LHP(2)/
     1   (1.+EXP(XIHP(2)*(XM(12,X)-XI(2,12,X)-XHP(2))))
C
	SS(1)=CI(2,2,X,T)
	SS(2)=CI(2,12,X,T)
	SS(3)=CM(2,X,T)
	SS(4)=CM(12,X,T)
	CALL HK(VHK,SS)
	ATMI(2,2,X)=ATMI(2,2,X) + 2.D0*NPHK(2)*VHK
	ATMI(2,12,X)=ATMI(2,12,X) - 2.D0*NPHK(2)*VHK
C
	FIKM(2,2,X,T)=FIKM(2,2,X,T)+AI0(2)*ATMI(2,2,X)
	FIKM(2,12,X,T)=FIKM(2,12,X,T)+AI0(2)*ATMI(2,12,X)
C
	DO 110 I=1,SOLS
  110   FKM(I,X,T)=FVM(X,T)*CM(I,X,T)
	FKM(13,X,T)=FVM(X,T)*IMPM(X,T)
C
	RETURN
	END
