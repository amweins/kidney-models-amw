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
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),ATIE(3,12),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3),JHK(3),JHP(3),JAE1(3),JTSC(3),JNHE3(3,3),
     1   NNKCC(3),NKCL(3),JNKCC(3,4),JKCC(3,3)
C
	DOUBLE PRECISION PHI(60),VHK,SHK(4),HMI11,
     1   VNHE3(3),SNHE3(3,2),VTSC,STSC(4),
     1   ZME(12),ZMI(3,12),ZIE(3,12),ZES(12),ZIS(3,12),
     1   SNKCC(8),VNKCC(4),SKCC(6),VKCC(3)
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
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,
     1   NNKCC,NKCL,JNKCC,JKCC
C
C
C COMPUTE THE RELEVANT VOLUMES
C
C COMPUTE THE MEMBRANE FLUXES
C
	FEVM(T)=AME*LPME*(RT*IMPS-RT*IMPM+PM-PS)/RT
	DO 30 K=1,1
	FIVM(K,T)=AI0(K)*LPMI(K)*(RT*IMPS-RT*IMPM+PM-PS)/RT
   30   CONTINUE
	DO 35 I=1,SOLS
	FEVM(T)=FEVM(T)+AME*LPME*SME(I)*(CS(I)-CM(I))
	DO 35 K=1,1
	FIVM(K,T)=FIVM(K,T)+AI0(K)*LPMI(K)*SMI(K,I)*(CS(I)-CM(I))
   35   CONTINUE
C
	DO 100 I=1,SOLS
C
	CME(I)=(CS(I)-CM(I))
	IF (CME(I)) 347,346,347
  346   CME(I)=CM(I)
	GO TO 351
  347   CME(I)=CME(I)/DLOG(CS(I)/CM(I))
  351   DO 360 K=1,1
	CMI(K,I)=(CS(I)-CM(I))
	IF (CMI(K,I)) 353,352,353
  352   CMI(K,I)=CM(I)
	GO TO 360
  353   CMI(K,I)=CMI(K,I)/DLOG(CS(I)/CM(I))
  360   CONTINUE
C
C
	XM(I)=RTE*DLOG(CM(I))+Z(I)*F*VM*1.D-6
	XS(I)=RTE*DLOG(CS(I))+Z(I)*F*VS*1.D-6
C
	ZME(I)=Z(I)*F*(VM-VS)*1.D-6/RTE
	DO 380 K=1,1
  380   ZMI(K,I)=Z(I)*F*(VM-VS)*1.D-6/RTE
C
C CONVECTIVE FLUXES
C
	FEKM(I,T)=FEVM(T)*(1.D0-SME(I))*CME(I)
	DO 390 K=1,1
  390   FIKM(K,I,T)=FIVM(K,T)*(1.D0-SMI(K,I))*CMI(K,I)
C
C GOLDMAN FLUXES
C
	IF (DABS(Z(I)).LT.0.1) GO TO 90
	IF (DABS(ZME(I)).GT.1.D-10) THEN
	FEKM(I,T)=FEKM(I,T)+HME(I)*AME*ZME(I)*
     1   (CM(I) - CS(I)*DEXP(-ZME(I)))/(1.D0 - DEXP(-ZME(I)))
	ELSE	
	FEKM(I,T)=FEKM(I,T)+HME(I)*AME*(CM(I) - CS(I))
	ENDIF
C
	DO 85 K=1,1
	IF (DABS(ZMI(K,I)).GT.1.D-10) THEN
	FIKM(K,I,T)=FIKM(K,I,T)+HMI(K,I)*AI0(K)*ZMI(K,I)*
     1   (CM(I) - CS(I)*DEXP(-ZMI(K,I)))/(1.D0 - DEXP(-ZMI(K,I)))
	ELSE
	FIKM(K,I,T)=FIKM(K,I,T)+HMI(K,I)*AI0(K)*(CM(I) - CS(I))
	ENDIF
   85   CONTINUE
	GO TO 100
C
   90   FEKM(I,T)=FEKM(I,T)+HME(I)*AME*(CM(I) - CS(I))
	DO 95 K=1,1
   95   FIKM(K,I,T)=FIKM(K,I,T)+HMI(K,I)*AI0(K)*(CM(I) - CS(I))
  100   CONTINUE
C
C
C NET COTRANSPORTERS
C
	DO 120 K=1,1
	DO 120 I=1,SOLS
	DO 120 J=1,SOLS
  120   FIKM(K,I,T)=FIKM(K,I,T)+LMI(K,I,J)*AI0(K)*(XM(J)-XS(J))
C
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C
C   THE OPEN CIRCUIT CONDITION
	CURE=0.D0
	DO 255 I=1,SOLS
  255   CURE=CURE+F*Z(I)*FEKM(I,T)
	PHI(1)=CURE
	DO 265 K=1,1
 	CURI(K)=0.D0
	DO 260 I=1,SOLS
  260   CURI(K)=CURI(K)+F*Z(I)*FIKM(K,I,T)
  265   PHI(1)=PHI(1)+CURI(K)
	RETURN
	END
