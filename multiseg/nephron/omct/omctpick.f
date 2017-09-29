	SUBROUTINE OMCTPICK
C SUBROUTINE TO SELECT INTRAEPITHELIAL MODEL VARIABLES FOR OUTPUT
C
	INTEGER SOLS,X,T,CHOP
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
	DOUBLE PRECISION CO2GEN,CO2FLX
C
	DOUBLE PRECISION TCO2,AEQ,BEQ,CEQ,PAEQ,PBEQ,
     1    ALPH,RHO,RHON,RHOP,NU,PHI,DPHDNU,PHEQ,
     1    OAEQ,OBEQ,OCEQ,OPAEQ,OPBEQ,OPHEQ
C
C	TCO2-	TOTAL CO2 IN DISQUILIBRIUM (AND EQUILIBRIUM) SOLUTIONS
C	AEQ-	EQUILIBRIUM CONCENTRATION OF H2CO3
C	BEQ-	EQUILIBRIUM CONCENTRATION OF HCO3
C	CEQ-	EQUILIBRIUM CONENTRATION OF CO2
C	PAEQ-	EQUILIBRIUM CONCENTRATION OF ACID PHOSPHATE
C	PBEQ-	EQUILIBRIUM CONCENTRATION OF BASIC PHOSPHATE
C	NU-	EXPONENTIATED DISEQUILIBRIUM PH
C	O...-	EQUILIBRIUM VALUES FOR AN OPEN SYSTEM
C
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
C
	CO2GEN = (KDHY(5)*CM(5,X,T) - KHY(5)*CM(6,X,T))*AM(X,T)
	CO2FLX = (FEKM(6,X,T) + FIKM(2,6,X,T))*SM(X,T)
C
C
C COMPUTE THE EQUILIBRIUM PH AND EQUILIBRIUM CO2
C
	ALPH = KHY(5)/KDHY(5)
	RHO  = 10.**(LCHM(X)-PKC)
	RHON = 10.**(LCHM(X)-PKN)
	RHOP = 10.**(LCHM(X)-PKP)
	TCO2 =CM(4,X,T) + CM(5,X,T) + CM(6,X,T)    
	NU = 1.0
C
C
   10   ITER=0
   15   ITER=ITER+1
	PHI = ALPH*NU*RHO*TCO2/(1. + ALPH + ALPH*NU*RHO)
     1   - CM(4,X,T) + (NU - 1.)*CM(7,X,T)/(1. + NU*RHOP)
     1    + (NU - 1.)*CM(10,X,T)/(1. + NU*RHON)
	IF(DABS(PHI).LT.1.D-16) GO TO 50
C
C PHI IS TOO LARGE.
C   WE DETERMINE THE DERIVATIVE OF PHI WITH RESPECT TO NU.
C
	DPHDNU = ALPH*(1.+ALPH)*RHO*TCO2/((1. + ALPH + ALPH*NU*RHO)**2)
     1    + (1. + RHOP)*CM(7,X,T)/((1. + NU*RHOP)**2)
     1    + (1. + RHON)*CM(10,X,T)/((1. + NU*RHON)**2)
C
C IT REMAINS TO COMPUTE C NDERIV-1(PHI) AT THE CORRECTION TO BE 
C  SUBTRACTED FROM NU.
C
	NU = NU - PHI/DPHDNU
	IF (ITER.LT.50) GO TO 15
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE
	PRINT 40
   40   FORMAT (' TOO MANY ITERATIONS IN OMCT:EQLIB')
	STOP
C
C THE SPATIAL STEP CONVERGED,AND
C THE PROGRAM PICKS UP AT THIS POINT:
C
   50   BEQ = ALPH*NU*RHO*TCO2/(1. + ALPH + ALPH*NU*RHO)
	AEQ = BEQ/(NU*RHO)
	CEQ = AEQ/ALPH
	PBEQ = NU*(1. + RHOP)*CM(7,X,T)/(1. + NU*RHOP) 
	PAEQ = PBEQ/(NU*RHOP)
	PHEQ = DLOG10(NU)
C
C
C NOW COMPUTE THE VALUES FOR AN OPEN SYSTEM
C
	OCEQ = CM(6,X,T)
	OAEQ = ALPH*OCEQ
C
  110   ITER=0
  115   ITER=ITER+1
	PHI = ALPH*NU*RHO*OCEQ  - CM(4,X,T) 
     1    + (NU - 1.)*CM(7,X,T)/(1. + NU*RHOP)
     1    + (NU - 1.)*CM(10,X,T)/(1. + NU*RHON)
	IF(DABS(PHI).LT.1.D-16) GO TO 150
C
C PHI IS TOO LARGE.
C   WE DETERMINE THE DERIVATIVE OF PHI WITH RESPECT TO NU.
C
	DPHDNU = ALPH*RHO*OCEQ
     1    + (1. + RHOP)*CM(7,X,T)/((1. + NU*RHOP)**2)
     1    + (1. + RHON)*CM(10,X,T)/((1. + NU*RHON)**2)
C
C IT REMAINS TO COMPUTE C NDERIV-1(PHI) AT THE CORRECTION TO BE 
C  SUBTRACTED FROM NU.
C
	NU = NU - PHI/DPHDNU
	IF (ITER.LT.50) GO TO 115
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE
	PRINT 140
  140   FORMAT (' TOO MANY ITERATIONS IN OMCT:EQLIB')
	STOP
C
C THE SPATIAL STEP CONVERGED,AND
C THE PROGRAM PICKS UP AT THIS POINT:
C
  150   OBEQ = ALPH*NU*RHO*OCEQ
	OPBEQ = NU*(1. + RHOP)*CM(7,X,T)/(1. + NU*RHOP) 
	OPAEQ = OPBEQ/(NU*RHOP)
	ONH3EQ = NU*(1. + RHON)*CM(10,X,T)/(1. + NU*RHON) 
	ONH4EQ = ONH3EQ/(NU*RHON)
	OPHEQ = DLOG10(NU)
C
C
  470   WRITE (12,480) DIST, 6.D4*CO2GEN, 
     1    6.D4*CO2GEN*DX,-PHEQ, 1.D3*CEQ,-OPHEQ
  480   FORMAT (8E16.8)
C
	RETURN
	END
