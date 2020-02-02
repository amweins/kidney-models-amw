	SUBROUTINE SFPCTFLUX
C
	INTEGER SOLS,T,X,CHOP
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,DIST,
     1   EPSI,DT,RTAU,TIME,L0,L(1601,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5),BCO2(3)
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
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   TQM(81,2),RM(81,2),MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0(3),NP0(3),NNHE30(3),LPMI0(3),LPIS0(3),
     1   HMI0(3,15),HIS0(3,15),LMI0(3,15,15),LIS0(3,15,15),
     1   QIAMM0,MVL,MVD,RMT0
C
	DOUBLE PRECISION VNHE3(3),SNHE3(3,2),
     1   ZME(15),ZMI(3,15),ZIE(3,15),ZES(15),ZIS(3,15)
C
        COMMON SOLS,T,X,CHOP,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,BCO2,L0,L,
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
	COMMON/TORQUE/
     1   TQM,RM,MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0,NP0,NNHE30,LPMI0,LPIS0,
     1   HMI0,HIS0,LMI0,LIS0,QIAMM0,MVL,MVD,RMT0
C
C
C COMPUTE THE RELEVANT VOLUMES
	AE(X,T)=AE0*(1.0+MUA*(PE(X,T)-PM(X,T)))
	IF((AE(X,T).LT.AE0)) AE(X,T)=AE0
	CHVL(X,T)=CHVL0*(1.0+MUV*(PE(X,T)-PM(X,T)))
	IF((CHVL(X,T).LT.CHVL0)) CHVL(X,T)=CHVL0
	L(X,T)=CHVL(X,T)
	DO 20 K=1,1
	CLVL(K,X,T)=CLVL0(K)*IMP0(K)/IMP(K,X)
	L(X,T)=L(X,T)+CLVL(K,X,T)
	PI(K,X,T)=PM(X,T)
C
	HCBUF(K,X,T)=(TBUF(K)*CLVL0(K)/CLVL(K,X,T))/
     1   (1.0 + (10.**(-PKB(K)))/CI(K,12,X,T))
	CBUF(K,X,T)=(TBUF(K)*CLVL0(K)/CLVL(K,X,T))-HCBUF(K,X,T)
   20   CONTINUE
	SM(X,T)=6.28D0*RMT0*MUM
	AM(1,T)=3.14D0*(RM(1,T)**2)*MUM
C
C
C COMPUTE THE MEMBRANE FLUXES
C
	FEVM(X,T)=AME*LPME*(PM(X,T)-PE(X,T)-RT*IMPM(X,T))/RT
	FEVS(X,T)=AE(X,T)*LPES*(RT*IMPS(X,T)+PE(X,T)-PS(X,T))/RT
	DO 30 K=1,1
	FIVM(K,X,T)=AMI(K)*LPMI(K)*
     1   (RT*IMP(K,X)-RT*IMPM(X,T)+PM(X,T)-PI(K,X,T))/RT
	FIVS(K,X,T)=AIS(K)*LPIS(K)*
     1   (PI(K,X,T)-PS(X,T)+RT*IMPS(X,T)-RT*IMP(K,X))/RT
   30   JV(K,X,T)=LPIS(K)*AIE(K)*(PI(K,X,T)-PE(X,T)-RT*IMP(K,X))/RT
	DO 35 I=1,SOLS
	FEVM(X,T)=FEVM(X,T)+AME*LPME*SME(I)*(CE(I,X,T)-CM(I,X,T))
	FEVS(X,T)=FEVS(X,T)+AE(X,T)*LPES*SES(I)*(CS(I,X,T)-CE(I,X,T))
	DO 35 K=1,1
	FIVM(K,X,T)=FIVM(K,X,T)+AMI(K)*LPMI(K)*
     1   SMI(K,I)*(CI(K,I,X,T)-CM(I,X,T))
	FIVS(K,X,T)=FIVS(K,X,T)+AIS(K)*LPIS(K)*
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
  351   DO 360 K=1,1
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
	DO 370 K=1,1
  370   XI(K,I,X)=RTE*DLOG(CI(K,I,X,T))+Z(I)*F*VI(K,X,T)*1.D-6
C
	ZME(I)=Z(I)*F*(VM(X,T)-VE(X,T))*1.D-6/RTE
	ZES(I)=Z(I)*F*(VE(X,T)-VS(X,T))*1.D-6/RTE
	DO 380 K=1,1
	ZMI(K,I)=Z(I)*F*(VM(X,T)-VI(K,X,T))*1.D-6/RTE
	ZIE(K,I)=Z(I)*F*(VI(K,X,T)-VE(X,T))*1.D-6/RTE
  380   ZIS(K,I)=Z(I)*F*(VI(K,X,T)-VS(X,T))*1.D-6/RTE
C
C CONVECTIVE FLUXES
C
	FEKM(I,X,T)=FEVM(X,T)*(1.D0-SME(I))*CME(I)
	FEKS(I,X,T)=FEVS(X,T)*(1.D0-SES(I))*CES(I)
	DO 390 K=1,1
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
	DO 85 K=1,1
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*AMI(K)*ZMI(K,I)*
     1   (CM(I,X,T) - CI(K,I,X,T)*DEXP(-ZMI(K,I)))/
     1     (1.D0 - DEXP(-ZMI(K,I)))
	JK(K,I,X,T)=JK(K,I,X,T)+HIS(K,I)*AIE(K)*ZIE(K,I)*
     1   (CI(K,I,X,T) - CE(I,X,T)*DEXP(-ZIE(K,I)))/
     1     (1.D0 - DEXP(-ZIE(K,I)))
   85   FIKS(K,I,X,T)=FIKS(K,I,X,T)+HIS(K,I)*AIS(K)*ZIS(K,I)*
     1   (CI(K,I,X,T) - CS(I,X,T)*DEXP(-ZIS(K,I)))/
     1     (1.D0 - DEXP(-ZIS(K,I)))
	GO TO 100
C
   90   FEKM(I,X,T)=FEKM(I,X,T)+HME(I)*AME*(CM(I,X,T) - CE(I,X,T))
	FEKS(I,X,T)=FEKS(I,X,T)+HES(I)*AE(X,T)*(CE(I,X,T) - CS(I,X,T))
	DO 95 K=1,1
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+HMI(K,I)*
     1   AMI(K)*(CM(I,X,T) - CI(K,I,X,T))
	JK(K,I,X,T)=JK(K,I,X,T)+
     1   HIS(K,I)*AIE(K)*(CI(K,I,X,T) - CE(I,X,T))
   95   FIKS(K,I,X,T)=FIKS(K,I,X,T)+
     1   HIS(K,I)*AIS(K)*(CI(K,I,X,T) - CS(I,X,T))
  100   CONTINUE
C
C
C NET COTRANSPORTERS
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	FIKM(K,I,X,T)=FIKM(K,I,X,T)+LMI(K,I,J)*AMI(K)*(XM(J,X)-XI(K,J,X))
	FIKS(K,I,X,T)=FIKS(K,I,X,T)+LIS(K,I,J)*AIS(K)*(XI(K,J,X)-XS(J,X))
  120   JK(K,I,X,T)=JK(K,I,X,T)+LIS(K,I,J)*AIE(K)*(XI(K,J,X)-XE(J,X))
C
C THE NAH EXCHANGER
C   TRANSLATE CONCENTRATIONS TO THE NAH MODEL
C
	SNHE3(1,1)=CM(1,X,T)
	SNHE3(2,1)=CM(12,X,T)
	SNHE3(3,1)=CM(11,X,T)
	SNHE3(1,2)=CI(1,1,X,T)
	SNHE3(2,2)=CI(1,12,X,T)
	SNHE3(3,2)=CI(1,11,X,T)
C
	CALL NHE3(VNHE3,SNHE3)
C
	JNHE3(1,1,X,T)= NNHE3(1)*AMI(1)*VNHE3(1)
	JNHE3(1,2,X,T)= NNHE3(1)*AMI(1)*VNHE3(2)
	JNHE3(1,3,X,T)= NNHE3(1)*AMI(1)*VNHE3(3)
C
	FIKM(1,1,X,T)=FIKM(1,1,X,T) + JNHE3(1,1,X,T)
	FIKM(1,12,X,T)=FIKM(1,12,X,T) + JNHE3(1,2,X,T)
	FIKM(1,11,X,T)=FIKM(1,11,X,T) + JNHE3(1,3,X,T)
C
C PUMPS
C
C  SODIUM PUMPS
C
	DO 125 K=1,1
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
	FIKS(K,1,X,T)=FIKS(K,1,X,T)+AIS(K)*ATIS(K,1,X)
	FIKS(K,2,X,T)=FIKS(K,2,X,T)+AIS(K)*ATIS(K,2,X)
	FIKS(K,11,X,T)=FIKS(K,11,X,T)+AIS(K)*ATIS(K,11,X)
  125   CONTINUE
C
C  PROTON PUMPS
C
 	ATMI(1,12,X)= - LHP(1)
     1       /(1.+EXP(XIHP(1)*(XM(12,X)-XI(1,12,X)-XHP(1))))
	JHP(1,X,T)= -AMI(1)*LHP(1)
     1       /(1.+EXP(XIHP(1)*(XM(12,X)-XI(1,12,X)-XHP(1))))
C
 	FIKM(1,12,X,T)=FIKM(1,12,X,T)+AMI(1)*ATMI(1,12,X)
C
C
	DO 110 I=1,SOLS
  110   FKM(I,X,T)=FVM(X,T)*CM(I,X,T)
	FKM(16,X,T)=FVM(X,T)*IMPM(X,T)
C
	RETURN
	END
