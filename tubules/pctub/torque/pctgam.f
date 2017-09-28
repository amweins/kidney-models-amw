	SUBROUTINE PCTGAM(SWG,GAMMA)
C PROGRAM TO SET AND RESET THE THE GAMMA ARRAY OF NN INDEPENDENT
C VARIABLES.  FOR SWG=1 THE GAMMA ARE DEFINED; FOR SWG=2, THE
C NAMED VARIABLES ARE RESET.
C
	INTEGER SOLS,T,SWG,X
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,L0,L(81,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(81,2),PM(81,2),CM(15,81,2),
     1   IMPM(81,2),LCHM(81),XM(15,81),
     1   VS(81,2),PS(81,2),CS(15,81,2),
     1   IMPS(81,2),LCHS(81),XS(15,81),
     1   TL,DX,RM0,MUM,ETA,
     1   SM(81,2),AM(81,2),FVM(81,2),FKM(16,81,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(81,2),MUA,CHVL0,CHVL(81,2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(81,2),PE(81,2),CE(15,81,2),LCHE(81),XE(15,81),
     1   FEVM(81,2),FEKM(15,81,2),FEVS(81,2),FEKS(15,81,2),CURE(81)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),CLVL0(3),IMP0(3),CLVL(3,81,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,81,2),HCBUF(3,81,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),ATMI(3,15,81),ATIS(3,15,81),
     1   VI(3,81,2),PI(3,81,2),CI(3,15,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,15,81),FIVM(3,81,2),FIKM(3,15,81,2),FIVS(3,81,2),
     1   FIKS(3,15,81,2),CURI(3,81),JV(3,81,2),JK(3,15,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2),QIAMM
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   TQM(81,2),RM(81,2),MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0(3),NP0(3),NNHE30(3),LPMI0(3),LPIS0(3),
     1   HMI0(3,15),HIS0(3,15),LMI0(3,15,15),LIS0(3,15,15),
     1   QIAMM0,MVL,MVD,RMT0
C
C
	DOUBLE PRECISION GAMMA
	DIMENSION GAMMA(80)
C
        COMMON SOLS,T,X,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   TL,DX,RM0,MUM,ETA,
     1   SM,AM,FVM,FKM,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AMI,AIS,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM
	COMMON/TORQUE/
     1   TQM,RM,MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0,NP0,NNHE30,LPMI0,LPIS0,
     1   HMI0,HIS0,LMI0,LIS0,QIAMM0,MVL,MVD,RMT0
C
C
	NCELLS=1
	GO TO (100,200) SWG
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
C	13-	HCO2-
C	14-	H2CO2
C	15-	GLUC
C
C
  100   VE(X,T)=GAMMA(1)
	PE(X,T)=GAMMA(2)
	DO 136 I=1,SOLS
  136   CE(I,X,T)=GAMMA(2+I)
C RESCALE
	CE(5,X,T)=CE(5,X,T)*1.D-6
	CE(10,X,T)=CE(10,X,T)*1.D-6
	CE(11,X,T)=CE(11,X,T)*1.D-3
	CE(12,X,T)=CE(12,X,T)*1.D-9
	CE(14,X,T)=CE(14,X,T)*1.D-6
C
	DO 145 K=1,1
	KK=K*(2+SOLS)
	VI(K,X,T)=GAMMA(KK+1)
	IMP(K,X)=GAMMA(KK+2)
	DO 140 I=1,SOLS
  140   CI(K,I,X,T)=GAMMA(KK+2+I)
C RESCALE
	CI(K,5,X,T)=CI(K,5,X,T)*1.D-6
	CI(K,10,X,T)=CI(K,10,X,T)*1.D-6
	CI(K,11,X,T)=CI(K,11,X,T)*1.D-3
	CI(K,12,X,T)=CI(K,12,X,T)*1.D-9
	CI(K,14,X,T)=CI(K,14,X,T)*1.D-6
  145   CONTINUE
	KK=(1+NCELLS)*(2+SOLS)
	VM(X,T)=GAMMA(KK+1)
	PM(X,T)=GAMMA(KK+2)
	DO 150 I=1,SOLS
  150   CM(I,X,T)=GAMMA(KK+2+I)
C RESCALE
	CM(5,X,T)=CM(5,X,T)*1.D-6
	CM(10,X,T)=CM(10,X,T)*1.D-6
	CM(11,X,T)=CM(11,X,T)*1.D-3
	CM(12,X,T)=CM(12,X,T)*1.D-9
	CM(14,X,T)=CM(14,X,T)*1.D-6
C
	FVM(X,T)=GAMMA(3*(2+SOLS)+1)*1.D-6*MUM
	IMPM(X,T)=GAMMA(3*(2+SOLS)+2)
	RETURN
C
C
  200   GAMMA(1)=VE(X,T)
	GAMMA(2)=PE(X,T)
	DO 236 I=1,SOLS
  236   GAMMA(2+I)=CE(I,X,T)
C SCALE
	GAMMA(2+5)=GAMMA(2+5)*1.D6
	GAMMA(2+10)=GAMMA(2+10)*1.D6
	GAMMA(2+11)=GAMMA(2+11)*1.D3
	GAMMA(2+12)=GAMMA(2+12)*1.D9
	GAMMA(2+14)=GAMMA(2+14)*1.D6
	DO 245 K=1,1
	KK=K*(2+SOLS)
	GAMMA(KK+1)=VI(K,X,T)
	GAMMA(KK+2)=IMP(K,X)
	DO 240 I=1,SOLS
  240   GAMMA(KK+2+I)=CI(K,I,X,T)
C SCALE
	GAMMA(KK+2+5)=GAMMA(KK+2+5)*1.D6
	GAMMA(KK+2+10)=GAMMA(KK+2+10)*1.D6
	GAMMA(KK+2+11)=GAMMA(KK+2+11)*1.D3
	GAMMA(KK+2+12)=GAMMA(KK+2+12)*1.D9
	GAMMA(KK+2+14)=GAMMA(KK+2+14)*1.D6
  245   CONTINUE
	KK=(1+NCELLS)*(2+SOLS)
	GAMMA(KK+1)=VM(X,T)
	GAMMA(KK+2)=PM(X,T)
	DO 250 I=1,SOLS
  250   GAMMA(KK+2+I)=CM(I,X,T)
C SCALE
	GAMMA(KK+2+5)=GAMMA(KK+2+5)*1.D6
	GAMMA(KK+2+10)=GAMMA(KK+2+10)*1.D6
	GAMMA(KK+2+11)=GAMMA(KK+2+11)*1.D3
	GAMMA(KK+2+12)=GAMMA(KK+2+12)*1.D9
	GAMMA(KK+2+14)=GAMMA(KK+2+14)*1.D6
C
	GAMMA(3*(2+SOLS)+1)=FVM(X,T)*1.D+6/MUM
	GAMMA(3*(2+SOLS)+2)=IMPM(X,T)
	RETURN
C
	END

