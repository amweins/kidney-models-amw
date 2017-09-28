	SUBROUTINE DCTORQUE(SWT,CHOP)
C
C SUBROUTINE TO SET COMPLIANCE AND TORQUE VARIABLES (SWT=1)
C  TO COMPUTE TUBULE RADIUS AND TORQUE, AND MODULATE PARAMETERS (SWT=2)
C   TO OUTPUT TORQUE-RELATED VARIABLES (SWT=3)
C
	INTEGER SOLS,TAU,T,TLIM,X,CHOP,SWT
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
     1   VI(3,81,2),PI(3,81,2),CI(3,12,81,2),IMP(3,81),LCHI(3,81),
     1   XI(3,12,81),FIVM(3,81,2),FIKM(3,12,81,2),FIVS(3,81,2),
     1   FIKS(3,12,81,2),CURI(3,81),JV(3,81,2),JK(3,12,81,2)
C SPECIAL TRANSPORTERS
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,81,2),JHK(3,81,2),JHP(3,81,2),
     1   JAE1(3,81,2),JTSC(3,81,2),JNHE3(3,3,81,2)
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   TQM(81,2),RM(81,2),MUR,VM0,AM0,TQM0,RMT0,
     1   SCALMI(3),SCALIS(3),LHP0(3),NP0(3),NNHE30(3),NTSC0(3),
     1   LPMI0(3),LPIS0(3),HMI0(3,15),HIS0(3,15),
     1   LMI0(3,15,15),LIS0(3,15,15)
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
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3
	COMMON/TORQUE/
     1   TQM,RM,MUR,SCALMI,SCALIS,VM0,AM0,TQM0,
     1   LHP0,NP0,NNHE30,NTSC0,
     1   LPMI0,LPIS0,HMI0,HIS0,LMI0,LIS0,RMT0
C
C
	GO TO (100,200,300) SWT
C
C FIRST PASS OF THE PROGRAM: INITIAL PARAMETERS ARE SET.
C
  100   CONTINUE
	DO 80 K=1,1
	LHP0(K)=LHP(K)
	NP0(K)=NP(K)
	NNHE30(K)=NNHE3(K)
	NTSC0(K)=NTSC(K)
	LPMI0(K)=LPMI(K)
	LPIS0(K)=LPIS(K)
	DO 72 I=1,SOLS
	HMI0(K,I)=HMI(K,I)
	HIS0(K,I)=HIS(K,I)
   72   CONTINUE
	DO 75 I=1,SOLS
	DO 75 J=1,SOLS
	LMI0(K,I,J)=LMI(K,I,J)
	LIS0(K,I,J)=LIS(K,I,J)
   75   CONTINUE
C
   80   CONTINUE
C
	SM(1,T)=6.28D0*RMT0*MUM
	RM(1,T)=RM0*(1.0 + MUR*(PM(1,T)-PS(1,T)))
	AM0=3.14D0*(RMT0**2)*MUM
	AM(1,T)=3.14D0*(RM(1,T)**2)*MUM
C
C THE REFERENCE TORQUE IS COMPUTED WITH RESPECT TO RMT0
C WHICH IS A DIFFERENT RADIUS FROM RM0, USED IN THE COMPLIANCE CALCULATION
C
	TQM0=VM0*MUM
	TQM(1,T)=(FVM(1,T)/AM(1,T))*MUM
C
	RETURN
C
C MDULATION OF MODEL PARAMTERS WITH TORQUE
C
  200   CONTINUE
C
	SM(X,T)=6.28D0*RMT0*MUM
	RM(X,T)=RM0*(1.0 + MUR*(PM(X,T)-PS(X,T)))
	AM(X,T)=3.14D0*(RM(X,T)**2)*MUM
	TQM(X,T)=(FVM(X,T)/AM(X,T))*MUM
C
	DO 180 K=1,1
	LHP(K)=LHP0(K)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	NP(K)=NP0(K)*(1.0 + SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
	NNHE3(K)=NNHE30(K)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	NTSC(K)=NTSC0(K)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	LPMI(K)=LPMI0(K)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	LPIS(K)=LPIS0(K)*(1.0 + SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
	DO 172 I=1,SOLS
	HMI(K,I)=HMI0(K,I)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	HIS(K,I)=HIS0(K,I)*(1.0 + SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
  172   CONTINUE
	DO 175 I=1,SOLS
	DO 175 J=1,SOLS
	LMI(K,I,J)=LMI0(K,I,J)*(1.0 + SCALMI(K)*(TQM(X,T)/TQM0 -1.0))
	LIS(K,I,J)=LIS0(K,I,J)*(1.0 + SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
  175   CONTINUE
C
C Independent modulation of either TSC or NHE:
C	NNHE3(K)=NNHE30(K)*(1.0 + 2.0*SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
C	NTSC(K)=NTSC0(K)*(1.0 + 2.0*SCALIS(K)*(TQM(X,T)/TQM0 -1.0))
C
  180   CONTINUE
C
	RETURN
C
C   TO OUTPUT TORQUE-RELATED VARIABLES (SWT=3)
C
  300   CONTINUE
C
	DO 210 K=1,CHOP+1
  210	WRITE(13,215) FLOAT(K-1)*DX,
     1   PM(K,T),1.D4*RM(K,T),AM(K,T)/AM0,
     1   FVM(K,T)/FVM(1,T),FVM(K,T)/AM(K,T),TQM(K,T)/TQM0
  215   FORMAT(8D16.8)
C
	RETURN
	END