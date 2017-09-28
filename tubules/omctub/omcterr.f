	SUBROUTINE OMCTERR(PHI)
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
	DOUBLE PRECISION PHI
	DIMENSION PHI(80)
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
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C
C FIRST FOR THE INTERSPACE
C   ELECTRONEUTRALITY
	PHI(1)=0.D0
	DO 130 I=1,SOLS
  130   PHI(1)=PHI(1)+Z(I)*CE(I,X,T)
C
	IF (T.EQ.1) GO TO 140
C
C   MASS CONSERVATION IN THE TIME-DEPENDENT CASE
	PHI(2)=FEVS(X,T)-FEVM(X,T)
     1   +RTAU*(CHVL(X,T)-CHVL(X,T-1))
	DO 135 I=1,SOLS
	PHI(2+I)=FEKS(I,X,T)-FEKM(I,X,T)
     1    +RTAU*(CE(I,X,T)*CHVL(X,T)-CE(I,X,T-1)*CHVL(X,T-1))
  135   CONTINUE
	DO 138 K=2,2
	PHI(2)=PHI(2) -JV(K,X,T)
	DO 136 I=1,SOLS
	PHI(2+I)=PHI(2+I) -JK(K,I,X,T)
  136   CONTINUE
  138   CONTINUE
	GO TO 149
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
  140   PHI(2)=FEVS(X,T)-FEVM(X,T)
	DO 145 I=1,SOLS
  145   PHI(2+I)=FEKS(I,X,T)-FEKM(I,X,T)
	DO 148 K=2,2
	PHI(2)=PHI(2)-JV(K,X,T)
	DO 146 I=1,SOLS
  146   PHI(2+I)=PHI(2+I)-JK(K,I,X,T)
  148   CONTINUE
  149   CONTINUE
C
	DO 200 K=2,2
	KK=1*(2+SOLS)
C
C THEN FOR THE CELLS
C   ELECTRONEUTRALITY
	PHI(1+KK)=IMP(K,X)*ZIMP(K)-CBUF(K,X,T)
	DO 150 I=1,SOLS
  150   PHI(1+KK)=PHI(1+KK)+Z(I)*CI(K,I,X,T)
C
	IF (T.EQ.1) GO TO 160
C
C   MASS CONSERVATION IN THE TIME-DEPENDENT CASE
	PHI(2+KK)=FIVS(K,X,T)-FIVM(K,X,T)
     1  +JV(K,X,T)+RTAU*(CLVL(K,X,T)-CLVL(K,X,T-1))
	DO 155 I=1,SOLS
	PHI(2+KK+I)=FIKS(K,I,X,T)-
     1   FIKM(K,I,X,T)+JK(K,I,X,T)+
     1     RTAU*(CI(K,I,X,T)*CLVL(K,X,T)-CI(K,I,X,T-1)*CLVL(K,X,T-1))
  155   CONTINUE
C
C   THE PROTON FLUX MUST INCLUDE THE CELLULAR BUFFERS
C
	PHI(2+KK+SOLS)=PHI(2+KK+SOLS)-
     1    RTAU*(CBUF(K,X,T)*CLVL(K,X,T)-CBUF(K,X,T-1)*CLVL(K,X,T-1))
C
	GO TO 169
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
  160   PHI(2+KK)=FIVS(K,X,T)-FIVM(K,X,T)+JV(K,X,T)
	DO 165 I=1,SOLS
  165   PHI(2+KK+I)=FIKS(K,I,X,T)-FIKM(K,I,X,T)+JK(K,I,X,T)
  169   CONTINUE
C
C
  200   CONTINUE
	DO 210 I=1,4*(SOLS+2)
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
	DO 220 K=2,2
	KK=1*(2+SOLS)
  220   PHI(2+SOLS+KK)=PHI(2+4+KK)+PHI(2+7+KK)-
     1    PHI(2+11+KK)-PHI(2+SOLS+KK)
C
C   TOTAL CO2, PHOSPHATE, AND AMMONIA CONTENT:
C
	PHI(2+4)=PHI(2+4)+PHI(2+5)+PHI(2+6)
	PHI(2+7)=PHI(2+7)+PHI(2+8)
	PHI(2+10)=PHI(2+10)+PHI(2+11)
	DO 230 K=2,2
	KK=1*(2+SOLS)
	PHI(2+4+KK)=PHI(2+4+KK)+PHI(2+5+KK)+PHI(2+6+KK)
	PHI(2+7+KK)=PHI(2+7+KK)+PHI(2+8+KK)
	PHI(2+10+KK)=PHI(2+10+KK)+PHI(2+11+KK)
  230   CONTINUE
C
C   CO2, PHOSPHATE, AND AMMMONIA PH EQUILIBRIUM: 
C
	LCHE(X)=-1.*DLOG10(CE(SOLS,X,T))
	PHI(2+5)=LCHE(X)-PKC-DLOG10(CE(4,X,T)/CE(5,X,T))
	PHI(2+8)=LCHE(X)-PKP-DLOG10(CE(7,X,T)/CE(8,X,T))
	PHI(2+11)=LCHE(X)-PKN-DLOG10(CE(10,X,T)/CE(11,X,T))
	DO 240 K=2,2
	KK=1*(2+SOLS)
	LCHI(K,X)=-1.*DLOG10(CI(K,SOLS,X,T))
	PHI(2+5+KK)=LCHI(K,X)-PKC-DLOG10(CI(K,4,X,T)/CI(K,5,X,T))
	PHI(2+8+KK)=LCHI(K,X)-PKP-DLOG10(CI(K,7,X,T)/CI(K,8,X,T))
	PHI(2+11+KK)=LCHI(K,X)-PKN-DLOG10(CI(K,10,X,T)/CI(K,11,X,T))
  240   CONTINUE
C
C   HYDRATION AND DHYDRATION OF CO2
C
	IF (T.EQ.1) GO TO 251
C
	PHI(2+6)=PHI(2+6)+
     1   1.D6*CHVL(X,T)*(KHY(4)*CE(6,X,T)-KDHY(4)*CE(5,X,T))
	DO 250 K=2,2
	KK=1*(2+SOLS)
  250   PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.D6*CLVL(K,X,T)*(KHY(K)*CI(K,6,X,T)-KDHY(K)*CI(K,5,X,T))
	GO TO 253
C
  251   PHI(2+6)=PHI(2+6)+
     1   1.D6*CHVL(X,T)*(KHY(4)*CE(6,X,T)-KDHY(4)*CE(5,X,T))
	DO 252 K=2,2
	KK=1*(2+SOLS)
  252   PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.D6*CLVL(K,X,T)*(KHY(K)*CI(K,6,X,T)-KDHY(K)*CI(K,5,X,T))
C
  253   CONTINUE
C
C   THE OPEN CIRCUIT CONDITION
	KK=2*(2+SOLS)
	CURE(X)=0.D0
	DO 255 I=1,SOLS
  255   CURE(X)=CURE(X)+F*Z(I)*FEKM(I,X,T)
	PHI(1+KK)=CURE(X)
	DO 259 K=2,2
 	CURI(K,X)=0.D0
	DO 257 I=1,SOLS
  257   CURI(K,X)=CURI(K,X)+F*Z(I)*FIKM(K,I,X,T)
  259   PHI(1+KK)=PHI(1+KK)+CURI(K,X)
C
	IF (X.EQ.1) RETURN
C
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C  FOR THE LUMEN VARIABLES
C
	IF (T.EQ.1) GO TO 270
C
C   MASS CONSERVATION IN THE TIME-DEPENDENT CASE
	PHI(2+KK)=(FVM(X,T)-FVM(X-1,T))/DX
     1  +SM(X,T)*FEVM(X,T)+
     1            +(FVM(X,T-1)-FVM(X-1,T-1))/DX
     1  +SM(X,T-1)*FEVM(X,T-1)+
     1       +2.D0*RTAU*(AM(X,T)-AM(X,T-1))
	DO 260 I=1,SOLS
	PHI(2+I+KK)=(FKM(I,X,T)-FKM(I,X-1,T))/DX
     1  +SM(X,T)*FEKM(I,X,T)+
     1            +(FKM(I,X,T-1)-FKM(I,X-1,T-1))/DX
     1  +SM(X,T-1)*FEKM(I,X,T-1)+
     1   +2.D0*RTAU*(AM(X,T)*CM(I,X,T)-AM(X,T-1)*CM(I,X,T-1))
  260   CONTINUE
	DO 262 K=2,2
	PHI(2+KK)=PHI(2+KK)+
     1   (SM(X,T)*FIVM(K,X,T)+
     1   SM(X,T-1)*FIVM(K,X,T-1))
	DO 262 I=1,SOLS
	PHI(2+I+KK)=PHI(2+I+KK)+
     1   (SM(X,T)*FIKM(K,I,X,T)+
     1   SM(X,T-1)*FIKM(K,I,X,T-1))
  262   CONTINUE
C
C   AND FOR THE IMPERMEANT
C
	PHI(3+SOLS+KK)=(FKM(13,X,T)-FKM(13,X-1,T))/DX
     1            +(FKM(13,X,T-1)-FKM(13,X-1,T-1))/DX
     1 +2.D0*RTAU*(AM(X,T)*IMPM(X,T)-AM(X,T-1)*IMPM(X,T-1))
C
	GO TO 290
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
  270   PHI(2+KK)=(FVM(X,T)-FVM(X-1,T))/DX
     1  +SM(X,T)*FEVM(X,T)
	DO 280 I=1,SOLS+1
	PHI(2+I+KK)=(FKM(I,X,T)-FKM(I,X-1,T))/DX
     1  +SM(X,T)*FEKM(I,X,T)
  280   CONTINUE
	DO 282 K=2,2
	PHI(2+KK)=PHI(2+KK)+
     1   (SM(X,T)*FIVM(K,X,T))
	DO 282 I=1,SOLS
	PHI(2+I+KK)=PHI(2+I+KK)+
     1   (SM(X,T)*FIKM(K,I,X,T))
  282   CONTINUE
C
C   AND FOR THE IMPERMEANT
C
	PHI(3+SOLS+KK)=(FKM(13,X,T)-FKM(13,X-1,T))/DX
C
  290   CONTINUE
	DO 295 I=2+KK,3+SOLS+KK
  295   PHI(I)=1.D6*PHI(I)
C
C
C   THE ERROR TERM FOR PROTON GENERATION IS REPLACED BY CONSERVATION
C   OF CHARGE IN THE BUFFER REACTIONS FOR SOLUTES 4,7,10 AND 12.
C
	PHI(2+SOLS+KK)=PHI(2+4+KK)+PHI(2+7+KK)-
     1    PHI(2+11+KK)-PHI(2+SOLS+KK)
C
C   TOTAL CO2, PHOSPHATE, AND AMMONIA CONTENT:
C
	PHI(2+4+KK)=PHI(2+4+KK)+PHI(2+5+KK)+PHI(2+6+KK)
	PHI(2+7+KK)=PHI(2+7+KK)+PHI(2+8+KK)
	PHI(2+10+KK)=PHI(2+10+KK)+PHI(2+11+KK)
C
C   CO2, PHOSPHATE, AND AMMMONIA PH EQUILIBRIUM: 
C
	LCHM(X)=-1.*DLOG10(CM(SOLS,X,T))
	PHI(2+5+KK)=LCHM(X)-PKC-DLOG10(CM(4,X,T)/CM(5,X,T))
	PHI(2+8+KK)=LCHM(X)-PKP-DLOG10(CM(7,X,T)/CM(8,X,T))
	PHI(2+11+KK)=LCHM(X)-PKN-DLOG10(CM(10,X,T)/CM(11,X,T))
C
C   HYDRATION AND DHYDRATION OF CO2
C
	IF (T.EQ.2) THEN
	PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.0D6*AM(X,T)*(KHY(5)*CM(6,X,T)-KDHY(5)*CM(5,X,T))+
     1   1.0D6*AM(X,T-1)*(KHY(5)*CM(6,X,T-1)-KDHY(5)*CM(5,X,T-1))
	ELSE
	PHI(2+6+KK)=PHI(2+6+KK)+
     1   1.0D6*AM(X,T)*(KHY(5)*CM(6,X,T)-KDHY(5)*CM(5,X,T))
	ENDIF
C
C   POISEIULLE FLOW
C
	PHI(4+SOLS+KK)=(PM(X,T)-PM(X-1,T))/DX
     1  +50.2*ETA*MUM*FVM(X,T)/(AM(X,T)**2)
C
	RETURN
	END
