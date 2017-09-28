	SUBROUTINE ERRVEC(PHI)
C
	INTEGER SOLS,T
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY(4),KDHY(4),L0,L(2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(12),IMPM,LCHM,XM(12),
     1   VS,PS,CS(12),IMPS,LCHS,XS(12)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(2),PE(2),CE(12,2),LCHE,XE(12),
     1   FEVM(2),FEKM(12,2),FEVS(2),FEKS(12,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C
	DOUBLE PRECISION PHI(60),VHK,SS(4),VAE1,HMI11,
     1   ZME(12),ZMI(3,12),ZIE(3,12),ZES(12),ZIS(3,12)
C
C
        COMMON SOLS,T,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
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
	AE(T)=AE0*(1.0+MUA*(PE(T)-PM))
	IF((AE(T).LT.AE0)) AE(T)=AE0
	CHVL(T)=CHVL0*(1.0+MUV*(PE(T)-PM))
	IF((CHVL(T).LT.CHVL0)) CHVL(T)=CHVL0
	L(T)=CHVL(T)
	DO 20 K=1,2
	L(T)=L(T)+CLVL(K,T)
	IMP(K)=CLVL0(K)*IMP0(K)/CLVL(K,T)
	PI(K,T)=PM
C
	HCBUF(K,T)=(TBUF(K)*CLVL0(K)/CLVL(K,T))/
     1   (1.0 + (10.**(-PKB(K)))/CI(K,12,T))
	CBUF(K,T)=(TBUF(K)*CLVL0(K)/CLVL(K,T))-HCBUF(K,T)
   20   CONTINUE
C
C
C COMPUTE THE MEMBRANE FLUXES
C
	FEVM(T)=AME*LPME*(PM-PE(T)-RT*IMPM)/RT
	FEVS(T)=AE(T)*LPES*(RT*IMPS+PE(T)-PS)/RT
	DO 30 K=1,2
	FIVM(K,T)=AI0(K)*LPMI(K)*(RT*IMP(K)-RT*IMPM+PM-PI(K,T))/RT
	FIVS(K,T)=AI0(K)*LPIS(K)*(PI(K,T)-PS+RT*IMPS-RT*IMP(K))/RT
   30   JV(K,T)=LPIS(K)*AIE(K)*(PI(K,T)-PE(T)-RT*IMP(K))/RT
	DO 35 I=1,SOLS
	FEVM(T)=FEVM(T)+AME*LPME*SME(I)*(CE(I,T)-CM(I))
	FEVS(T)=FEVS(T)+AE(T)*LPES*SES(I)*(CS(I)-CE(I,T))
	DO 35 K=1,2
	FIVM(K,T)=FIVM(K,T)+AI0(K)*LPMI(K)*SMI(K,I)*(CI(K,I,T)-CM(I))
	FIVS(K,T)=FIVS(K,T)+AI0(K)*LPIS(K)*SIS(K,I)*(CS(I)-CI(K,I,T))
	JV(K,T)=JV(K,T)+LPIS(K)*AIE(K)*SIS(K,I)*(CE(I,T)-CI(K,I,T))
   35   CONTINUE
C
C THE CIVAN AND BOOKMAN MODIFICATION TO HMI(1,1)
	HMI11=HMI(1,1)
	HMI(1,1)=HMI11*((1.+CM(1)/.03)**(-1))*
     1      (1.-CI(1,1,T)/.05)*1.666667
C
	DO 100 I=1,SOLS
C
C
	CME(I)=(CE(I,T)-CM(I))
	IF (CME(I)) 347,346,347
  346   CME(I)=CM(I)
	GO TO 348
  347   CME(I)=CME(I)/DLOG(CE(I,T)/CM(I))
  348   CES(I)=(CE(I,T)-CS(I))
	IF (CES(I)) 350,349,350
  349   CES(I)=CS(I)
	GO TO 351
  350   CES(I)=CES(I)/DLOG(CE(I,T)/CS(I))
  351   DO 360 K=1,2
	CMI(K,I)=(CI(K,I,T)-CM(I))
	IF (CMI(K,I)) 353,352,353
  352   CMI(K,I)=CM(I)
	GO TO 354
  353   CMI(K,I)=CMI(K,I)/DLOG(CI(K,I,T)/CM(I))
  354   CIS(K,I)=(CI(K,I,T)-CS(I))
	IF (CIS(K,I)) 356,355,356
  355   CIS(K,I)=CS(I)
	GO TO 357
  356   CIS(K,I)=CIS(K,I)/DLOG(CI(K,I,T)/CS(I))
  357   CIE(K,I)=(CI(K,I,T)-CE(I,T))
	IF (CIE(K,I)) 359,358,359
  358   CIE(K,I)=CI(K,I,T)
	GO TO 360
  359   CIE(K,I)=CIE(K,I)/DLOG(CI(K,I,T)/CE(I,T))
  360   CONTINUE
C
C
	XE(I)=RTE*DLOG(CE(I,T))+Z(I)*F*VE(T)*1.D-6
	XM(I)=RTE*DLOG(CM(I))+Z(I)*F*VM*1.D-6
	XS(I)=RTE*DLOG(CS(I))+Z(I)*F*VS*1.D-6
	DO 370 K=1,2
  370   XI(K,I)=RTE*DLOG(CI(K,I,T))+Z(I)*F*VI(K,T)*1.D-6
C
	ZME(I)=Z(I)*F*(VM-VE(T))*1.D-6/RTE
	ZES(I)=Z(I)*F*(VE(T)-VS)*1.D-6/RTE
	DO 380 K=1,2
	ZMI(K,I)=Z(I)*F*(VM-VI(K,T))*1.D-6/RTE
	ZIE(K,I)=Z(I)*F*(VI(K,T)-VE(T))*1.D-6/RTE
  380   ZIS(K,I)=Z(I)*F*(VI(K,T)-VS)*1.D-6/RTE
C
C CONVECTIVE FLUXES
C
	FEKM(I,T)=FEVM(T)*(1.D0-SME(I))*CME(I)
	FEKS(I,T)=FEVS(T)*(1.D0-SES(I))*CES(I)
	DO 390 K=1,2
	FIKM(K,I,T)=FIVM(K,T)*(1.D0-SMI(K,I))*CMI(K,I)
	FIKS(K,I,T)=FIVS(K,T)*(1.D0-SIS(K,I))*CIS(K,I)
  390   JK(K,I,T)=JV(K,T)*(1.D0-SIS(K,I))*CIE(K,I)
C
C GOLDMAN FLUXES
C
	IF (DABS(Z(I)).LT.0.1) GO TO 90
	FEKM(I,T)=FEKM(I,T)+HME(I)*AME*ZME(I)*
     1   (CM(I) - CE(I,T)*DEXP(-ZME(I)))/(1.D0 - DEXP(-ZME(I)))
	FEKS(I,T)=FEKS(I,T)+HES(I)*AE(T)*ZES(I)*
     1   (CE(I,T) - CS(I)*DEXP(-ZES(I)))/(1.D0 - DEXP(-ZES(I)))
	DO 85 K=1,2
	FIKM(K,I,T)=FIKM(K,I,T)+HMI(K,I)*AI0(K)*ZMI(K,I)*
     1   (CM(I) - CI(K,I,T)*DEXP(-ZMI(K,I)))/(1.D0 - DEXP(-ZMI(K,I)))
	JK(K,I,T)=JK(K,I,T)+HIS(K,I)*AIE(K)*ZIE(K,I)*
     1   (CI(K,I,T) - CE(I,T)*DEXP(-ZIE(K,I)))/(1.D0 - DEXP(-ZIE(K,I)))
   85   FIKS(K,I,T)=FIKS(K,I,T)+HIS(K,I)*AI0(K)*ZIS(K,I)*
     1   (CI(K,I,T) - CS(I)*DEXP(-ZIS(K,I)))/(1.D0 - DEXP(-ZIS(K,I)))
	GO TO 100
C
   90   FEKM(I,T)=FEKM(I,T)+HME(I)*AME*(CM(I) - CE(I,T))
	FEKS(I,T)=FEKS(I,T)+HES(I)*AE(T)*(CE(I,T) - CS(I))
	DO 95 K=1,2
	FIKM(K,I,T)=FIKM(K,I,T)+HMI(K,I)*AI0(K)*(CM(I) - CI(K,I,T))
	JK(K,I,T)=JK(K,I,T)+HIS(K,I)*AIE(K)*(CI(K,I,T) - CE(I,T))
   95   FIKS(K,I,T)=FIKS(K,I,T)+HIS(K,I)*AI0(K)*(CI(K,I,T) - CS(I))
  100   CONTINUE
C
C THE CIVAN AND BOOKMAN MODIFICATION TO HMI(1,1)
	HMI(1,1)=HMI11
C
C
C NET COTRANSPORTERS
C
	DO 120 K=1,2
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
	FIKM(K,I,T)=FIKM(K,I,T)+LMI(K,I,J)*AI0(K)*(XM(J)-XI(K,J))
	FIKS(K,I,T)=FIKS(K,I,T)+LIS(K,I,J)*AI0(K)*(XI(K,J)-XS(J))
  120   JK(K,I,T)=JK(K,I,T)+LIS(K,I,J)*AIE(K)*(XI(K,J)-XE(J))
C
C AE1 COTRANSPORTERS
C
	SS(1)=CI(2,3,T)
	SS(2)=CI(2,4,T)
	SS(3)=CE(3,T)
	SS(4)=CE(4,T)
	CALL AE1(VAE1,SS)
	JK(2,3,T)=JK(2,3,T)-NAE1(2)*AIE(2)*VAE1
	JK(2,4,T)=JK(2,4,T)+NAE1(2)*AIE(2)*VAE1
	SS(3)=CS(3)
	SS(4)=CS(4)
	CALL AE1(VAE1,SS)
	FIKS(2,3,T)=FIKS(2,3,T)-NAE1(2)*AI0(2)*VAE1
	FIKS(2,4,T)=FIKS(2,4,T)+NAE1(2)*AI0(2)*VAE1
C
C
C PUMPS
C
C  SODIUM PUMPS
C
	DO 125 K=1,2
        KNPN(K)=0.0002*(1.D0+CI(K,2,T)/.00833)
        KNPK(K)=0.0001*(1.D0+CS(1)/.0185)
        ATIS(K,1)=NP(K)*(CI(K,1,T)/(KNPN(K)+CI(K,1,T)))**3*
     1    (CS(2)/(KNPK(K)+CS(2)))**2
C  ALLOW FOR COMPETITION BETWEEN K+ AND NH4+ 
	ATIS(K,2)=-ATIS(K,1)*0.667*CE(2,T)/(CE(2,T)+CE(11,T)/KNH4(K))
	ATIS(K,11)=-ATIS(K,1)*0.667*CE(11,T)/(KNH4(K)*CE(2,T)+CE(11,T))
C
	JK(K,1,T)=JK(K,1,T)+AIE(K)*ATIS(K,1)
	JK(K,2,T)=JK(K,2,T)+AIE(K)*ATIS(K,2)
	JK(K,11,T)=JK(K,11,T)+AIE(K)*ATIS(K,11)
	FIKS(K,1,T)=FIKS(K,1,T)+AI0(K)*ATIS(K,1)
	FIKS(K,2,T)=FIKS(K,2,T)+AI0(K)*ATIS(K,2)
	FIKS(K,11,T)=FIKS(K,11,T)+AI0(K)*ATIS(K,11)
  125   CONTINUE
C
C  PROTON PUMPS
C
	SS(1)=CI(3,2,T)
	SS(2)=CI(3,12,T)
	SS(3)=CM(2)
	SS(4)=CM(12)
	CALL HK(VHK,SS)
	ATMI(3,2)=  2.D0*NPHK(3)*VHK
	ATMI(3,12)=-2.D0*NPHK(3)*VHK
	ATIS(3,12)=LHP(3)/(1.+EXP(XIHP(3)*(XI(3,12)-XS(12)-XHP(3))))
C
	JK(3,12,T)=JK(3,12,T)+AIE(3)*ATIS(3,12)
	FIKS(3,12,T)=FIKS(3,12,T)+AI0(3)*ATIS(3,12)
	FIKM(3,2,T)=FIKM(3,2,T)+AI0(3)*ATMI(3,2)
	FIKM(3,12,T)=FIKM(3,12,T)+AI0(3)*ATMI(3,12)
C
C
	ATMI(2,2)=0.D0
 	ATMI(2,12)=-LHP(2)/(1.+EXP(XIHP(2)*(XM(12)-XI(2,12)-XHP(2))))
C
	SS(1)=CI(2,2,T)
	SS(2)=CI(2,12,T)
	SS(3)=CM(2)
	SS(4)=CM(12)
	CALL HK(VHK,SS)
 	ATMI(2,2)=ATMI(2,2) + 2.D0*NPHK(2)*VHK
	ATMI(2,12)=ATMI(2,12) - 2.D0*NPHK(2)*VHK
C
	FIKM(2,2,T)=FIKM(2,2,T)+AI0(2)*ATMI(2,2)
	FIKM(2,12,T)=FIKM(2,12,T)+AI0(2)*ATMI(2,12)
C
C
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C
C FIRST FOR THE INTERSPACE
C   ELECTRONEUTRALITY
	PHI(1)=0.D0
	DO 130 I=1,SOLS
  130   PHI(1)=PHI(1)+Z(I)*CE(I,T)
C
	IF (T.EQ.1) GO TO 140
C
C   MASS CONSERVATION IN THE TIME-DEPENDENT CASE
	PHI(2)=FEVS(T)-FEVM(T)
     1   +RTAU*(CHVL(T)-CHVL(T-1))
	DO 135 I=1,SOLS
	PHI(2+I)=FEKS(I,T)-FEKM(I,T)
     1    +RTAU*(CE(I,T)*CHVL(T)-CE(I,T-1)*CHVL(T-1))
  135   CONTINUE
	DO 138 K=1,2
	PHI(2)=PHI(2) -JV(K,T)
	DO 136 I=1,SOLS
	PHI(2+I)=PHI(2+I) -JK(K,I,T)
  136   CONTINUE
  138   CONTINUE
	GO TO 149
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
  140   PHI(2)=FEVS(T)-FEVM(T)
	DO 145 I=1,SOLS
  145   PHI(2+I)=FEKS(I,T)-FEKM(I,T)
	DO 148 K=1,2
	PHI(2)=PHI(2)-JV(K,T)
	DO 146 I=1,SOLS
  146   PHI(2+I)=PHI(2+I)-JK(K,I,T)
  148   CONTINUE
  149   CONTINUE
C
	DO 200 K=1,2
	KK=K*(2+SOLS)
C
C THEN FOR THE CELLS
C   ELECTRONEUTRALITY
	PHI(1+KK)=IMP(K)*ZIMP(K)-CBUF(K,T)
	DO 150 I=1,SOLS
  150   PHI(1+KK)=PHI(1+KK)+Z(I)*CI(K,I,T)
C
	IF (T.EQ.1) GO TO 160
C
C   MASS CONSERVATION IN THE TIME-DEPENDENT CASE
	PHI(2+KK)=FIVS(K,T)-FIVM(K,T)
     1  +JV(K,T)+RTAU*(CLVL(K,T)-CLVL(K,T-1))
	DO 155 I=1,SOLS
	PHI(2+KK+I)=FIKS(K,I,T)-
     1   FIKM(K,I,T)+JK(K,I,T)+
     1     RTAU*(CI(K,I,T)*CLVL(K,T)-CI(K,I,T-1)*CLVL(K,T-1))
  155   CONTINUE
C
C   THE PROTON FLUX MUST INCLUDE THE CELLULAR BUFFERS
C
	PHI(2+KK+SOLS)=PHI(2+KK+SOLS)-
     1    RTAU*(CBUF(K,T)*CLVL(K,T)-CBUF(K,T-1)*CLVL(K,T-1))
C
	GO TO 169
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
  160   PHI(2+KK)=FIVS(K,T)-FIVM(K,T)+JV(K,T)
	DO 165 I=1,SOLS
  165   PHI(2+KK+I)=FIKS(K,I,T)-FIKM(K,I,T)+JK(K,I,T)
  169   CONTINUE
C
C
  200   CONTINUE
	DO 210 I=1,3*(SOLS+2)
  210   PHI(I)=1.D6*PHI(I)
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
C   THE ERROR TERM FOR PROTON GENERATION IS REPLACED BY CONSERVATION
C   OF CHARGE IN THE BUFFER REACTIONS FOR SOLUTES 4,7,10 AND 12.
C
	PHI(2+SOLS)=PHI(2+4)+PHI(2+7)-PHI(2+11)-PHI(2+SOLS)
	DO 220 K=1,2
	KK=K*(2+SOLS)
  220   PHI(2+SOLS+KK)=PHI(2+4+KK)+PHI(2+7+KK)-
     1    PHI(2+11+KK)-PHI(2+SOLS+KK)
C
C   TOTAL CO2, PHOSPHATE, AND AMMONIA CONTENT:
C
	PHI(2+4)=PHI(2+4)+PHI(2+5)+PHI(2+6)
	PHI(2+7)=PHI(2+7)+PHI(2+8)
	PHI(2+10)=PHI(2+10)+PHI(2+11)
	DO 230 K=1,2
	KK=K*(2+SOLS)
	PHI(2+4+KK)=PHI(2+4+KK)+PHI(2+5+KK)+PHI(2+6+KK)
	PHI(2+7+KK)=PHI(2+7+KK)+PHI(2+8+KK)
	PHI(2+10+KK)=PHI(2+10+KK)+PHI(2+11+KK)
  230   CONTINUE
C
C   CO2, PHOSPHATE, AND AMMMONIA PH EQUILIBRIUM: 
C
	LCHE=-1.*DLOG10(CE(SOLS,T))
	PHI(2+5)=LCHE-PKC-DLOG10(CE(4,T)/CE(5,T))
	PHI(2+8)=LCHE-PKP-DLOG10(CE(7,T)/CE(8,T))
	PHI(2+11)=LCHE-PKN-DLOG10(CE(10,T)/CE(11,T))
	DO 240 K=1,2
	KK=K*(2+SOLS)
	LCHI(K)=-1.*DLOG10(CI(K,SOLS,T))
	PHI(2+5+KK)=LCHI(K)-PKC-DLOG10(CI(K,4,T)/CI(K,5,T))
	PHI(2+8+KK)=LCHI(K)-PKP-DLOG10(CI(K,7,T)/CI(K,8,T))
	PHI(2+11+KK)=LCHI(K)-PKN-DLOG10(CI(K,10,T)/CI(K,11,T))
  240   CONTINUE
C
C   HYDRATION AND DHYDRATION OF CO2
C
	IF (T.EQ.1) GO TO 251
C
	PHI(2+6)=PHI(2+6)+
     1   1.D6*CHVL(T)*(KHY(4)*CE(6,T)-KDHY(4)*CE(5,T))
	DO 250 K=1,2
	KK=K*(2+SOLS)
  250   PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.D6*CLVL(K,T)*(KHY(K)*CI(K,6,T)-KDHY(K)*CI(K,5,T))
	GO TO 253
C
  251   PHI(2+6)=PHI(2+6)+
     1   1.D6*CHVL(T)*(KHY(4)*CE(6,T)-KDHY(4)*CE(5,T))
	DO 252 K=1,2
	KK=K*(2+SOLS)
  252   PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.D6*CLVL(K,T)*(KHY(K)*CI(K,6,T)-KDHY(K)*CI(K,5,T))
C
  253   CONTINUE
C
C   THE OPEN CIRCUIT CONDITION
	CURE=0.D0
	DO 255 I=1,SOLS
  255   CURE=CURE+F*Z(I)*FEKM(I,T)
	PHI(1+3*(2+SOLS))=CURE
	DO 265 K=1,2
 	CURI(K)=0.D0
	DO 260 I=1,SOLS
  260   CURI(K)=CURI(K)+F*Z(I)*FIKM(K,I,T)
  265   PHI(1+3*(2+SOLS))=PHI(1+3*(2+SOLS))+CURI(K)
	RETURN
	END
