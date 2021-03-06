	SUBROUTINE NKCC(ISOFM,SS,VK,VNH4)
C
C IN THIS VERSION SUBSTRATE BINDING IS ASSUMED RAPID 
C   RELATIVE TO TRANSLOCATION
C CODE IS DERIVED FROM TSC.F
C
	INTEGER SOLS,SUBS,VELS
C
C       SOLS-   NUMBER OF SOLUTES
C       SUBS-   NUMBER OF SUBSTRATES
C       VELS-   NUMBER OF REACTIONS
C
	CHARACTER*1 ISOFM
C
C	ISOFM-	IDENTIFIES THE NKCC ISOFORM (B, A, OR F)
C
	DOUBLE PRECISION
     1   V,C,S,KF,KB,KEQ,SS,VK,VNH4
C
C       S-      SUBSTRATE CONCENTRATION (MMOL/ML)
C	KF-	FORWARD RATE CONSTANT = "Kon" (mol/l)^-1 * (s^-1)
C	KB-	BACKWARD RATE CONSTANT "Koff" (s^-1)
C	KEQ-	EQUILIBRIUM CONSTANT = KF/KB = Kon/Koff
C       C-	UNKNOWN CONCENTRATION (MMOL/ML)
C       V-	NET FORWARD REACTION VELOCITY
C	VK-	NET FLUX OF K (IN NA-K-2CL)  = KF(5)*C(5) - KB(5)*C(10)
C	VNH4-	NET FLUX OF NH4 (IN NA-NH4-2CL) = KF(13)*C(12) - KB(13)*C(15)
C       SS-     SUBSTRATE CONCENTRATION (MMOL/ML) -- IDENTICAL TO S
C
C Transported Solutes:
C	S1-	Na(1)
C	S2-	Na(2)
C	S3-	K(1)
C	S4-	K(2)
C	S5-	Cl(1)
C	S6-	Cl(2)
C	S7-	NH4(1)
C	S8-	NH4(2)
C
	DIMENSION
     1   S(10),KF(20),KB(20),KEQ(20),V(20),C(20),SS(10)
C
	DOUBLE PRECISION
     1    KFB(16),KBB(16),KEQB(16),
     1    KFA(16),KBA(16),KEQA(16),
     1    KFF(16),KBF(16),KEQF(16)
C
C PARAMETER SETS FOR EACH OF THE ISOFORMS
C
	DOUBLE PRECISION PHI,GAMMA,MAXPHI,PREMAX
	REAL*8 DELGAM,NDERIV,RPHI,EPSI
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C       PREMAX- MAXPHI FOR THE PREVIOUS NEWTON ITERATE
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C       EPSI-   TOLERANCE FOR THE ERROR VECTOR
C
	DIMENSION PHI(20),GAMMA(20),NDERIV(20,20),
     1    RPHI(20),DELGAM(20)
C
	REAL*8 SW,PTOL
	INTEGER PR,PFL
C
	DIMENSION PR(20),PFL(20)
C
        COMMON/NKCCKIN/ SOLS,SUBS,VELS,
     1   V,C,S,KF,KB,KEQ
C
C
	DATA KFB/
     1   .1000D+09 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+05 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+09 , .2517D+06 , .1000D+09 , .1000D+09 ,
     1   .2000D+04 , .1000D+09 , .1000D+09 , .1000D+09/
	DATA KBB/
     1   .2750D+08 , .8157D+04 , .5578D+09 , .8157D+04 ,
     1   .9695D+04 , .2750D+08 , .8157D+04 , .5578D+09 ,
     1   .8157D+04 , .2596D+06 , .5578D+09 , .8157D+04 ,
     1   .1939D+04 , .2750D+08 , .8157D+04 , .5578D+09/
	DATA KEQB/
     1   .3636D+01 , .1226D+05 , .1793D+00 , .1226D+05 ,
     1   .1031D+01 , .3636D+01 , .1226D+05 , .1793D+00 ,
     1   .1226D+05 , .9695D+00 , .1793D+00 , .1226D+05 ,
     1   .1031D+01 , .3636D+01 , .1226D+05 , .1793D+00/
	DATA KFA/
     1   .1000D+09 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+05 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+09 , .7535D+05 , .1000D+09 , .1000D+09 ,
     1   .2000D+04 , .1000D+09 , .1000D+09 , .1000D+09/
	DATA KBA/
     1   .1188D+08 , .8834D+04 , .1871D+10 , .8834D+04 ,
     1   .2904D+04 , .1188D+08 , .8834D+04 , .1871D+10 ,
     1   .8834D+04 , .2594D+06 , .1871D+10 , .8834D+04 ,
     1   .5808D+03 , .1188D+08 , .8834D+04 , .1871D+10/
	DATA KEQA/
     1   .8415D+01 , .1132D+05 , .5344D-01 , .1132D+05 ,
     1   .3443D+01 , .8415D+01 , .1132D+05 , .5344D-01 ,
     1   .1132D+05 , .2904D+00 , .5344D-01 , .1132D+05 ,
     1   .3443D+01 , .8415D+01 , .1132D+05 , .5344D-01/
	DATA KFF/
     1   .1000D+09 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+05 , .1000D+09 , .1000D+09 , .1000D+09 ,
     1   .1000D+09 , .3928D+05 , .1000D+09 , .1000D+09 ,
     1   .2000D+04 , .1000D+09 , .1000D+09 , .1000D+09/
	DATA KBF/
     1   .5892D+07 , .1312D+07 , .9149D+06 , .1312D+07 ,
     1   .1098D+04 , .5892D+07 , .1312D+07 , .9149D+06 ,
     1   .1312D+07 , .3578D+06 , .9149D+06 , .1312D+07 ,
     1   .2196D+03 , .5892D+07 , .1312D+07 , .9149D+06/
	DATA KEQF/
     1   .1697D+02 , .7623D+02 , .1093D+03 , .7623D+02 ,
     1   .9110D+01 , .1697D+02 , .7623D+02 , .1093D+03 ,
     1   .7623D+02 , .1098D+00 , .1093D+03 , .7623D+02 ,
     1   .9110D+01 , .1697D+02 , .7623D+02 , .1093D+03/
C
	DATA GAMMA/20*0.1D0/
	SOLS=15
	SUBS=8
	VELS=16
C
C  SOLUTE INDEX:
C
C Variable Concentrations:
C	C1-	E1
C	C2-	E1-Na
C	C3-	E1-NaCl
C	C4-	E1-NaKCl
C	C5-	E1-NaKCl2
C	C6-	E2
C	C7-	E2-Cl
C	C8-	E2-KCl
C	C9-	E2-KCl2
C	C10-	E2-NaKCl2
C
C	C11-	E1-NaNH4Cl
C	C12-	E1-NaNH4Cl2
C	C13-	E2-NH4Cl
C	C14-	E2-NH4Cl2
C	C15-	E2-NaNH4Cl2
C
C Transported Solutes:
C	S1-	Na(1)
C	S2-	Na(2)
C	S3-	K(1)
C	S4-	K(2)
C	S5-	Cl(1)
C	S6-	Cl(2)
C	S7-	NH4(1)
C	S8-	NH4(2)
C
C Reactions:
C	V1-	C1 + S1 -> C2
C	V2-	C2 + S5 -> C3
C	V3-	C3 + S3 -> C4
C	V4-	C4 + S5 -> C5
C	V5-	C5 -> C10
C
C	V6-	C9 + S2 -> C10
C	V7-	C8 + S6 -> C9
C	V8-	C7 + S4 -> C8
C	V9-	C6 + S6 -> C7
C	V10-	C6 -> C1
C
C	V11-	C3 + S7 -> C11
C	V12-	C11 + S5 -> C12
C	V13-	C12 -> C15
C	V14-	C14 + S2 -> C15
C	V15-	C13 + S6 -> C14
C	V16-	C7 + S8 -> C13
C
C
	NN=SOLS
	EPSI=1.D-10
	PTOL=1.D-15
	M1=20
	IERR=0
C
C SELECT PARAMETERS FOR THE ISOFORM
C
	IF (ISOFM.EQ.'B') THEN
	DO 12 I=1,VELS
	KF(I) = KFB(I)
	KB(I) = KBB(I)
	KEQ(I) = KEQB(I)
   12   CONTINUE
	ELSE IF (ISOFM.EQ.'A') THEN
	DO 14 I=1,VELS
	KF(I) = KFA(I)
	KB(I) = KBA(I)
	KEQ(I) = KEQA(I)
   14   CONTINUE
	ELSE IF (ISOFM.EQ.'F') THEN
	DO 16 I=1,VELS
	KF(I) = KFF(I)
	KB(I) = KBF(I)
	KEQ(I) = KEQF(I)
   16   CONTINUE
	ELSE
	PRINT 18
   18   FORMAT('ERROR IN READING NKCC ISOFORM')
	STOP
	ENDIF
C
	DO 20 I=1,SUBS
   20   S(I) = SS(I)
C
C
   70   ITER=0
   80   ITER=ITER+1
	IF (ITER.GT.1) PREMAX=MAXPHI
	DO 36 I=1,NN
   36   C(I)=GAMMA(I)
	CALL NKCCERR(PHI)
	MAXPHI=DABS(PHI(1))
	DO 230 J=2,NN
  230   IF (MAXPHI.LT.DABS(PHI(J))) MAXPHI=DABS(PHI(J))
C	PRINT 232, ITER, MAXPHI
  232   FORMAT ('ITER=',I3,5X,'MAXPHI=',D12.4)
	IF(MAXPHI.LT.EPSI) GO TO 558
C
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
	CALL NKCCJAC(NDERIV,M1)
	GO TO 923
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C
  923   DO 924 J=1,NN
  924   RPHI(J)=PHI(J)
C
	CALL LES8(NDERIV,RPHI,PR,PFL,NN,SW,PTOL,DELGAM,M1)
	IF (SW) 930,982,930
C
  930   DO 940 J=1,NN
  940   GAMMA(J)=GAMMA(J)-DELGAM(J)
C
	IF (ITER.LT.25) GO TO 80
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE WITHIN 15 PASSES
C THE PROGRAM STOPS.
	PRINT 980
  980   FORMAT (' TOO MANY ITERATIONS IN NKCC')
  982   STOP
C
  558   VK = KF(5)*C(5) - KB(5)*C(10)
	VNH4 = KF(13)*C(12) - KB(13)*C(15)
	RETURN
	END
	SUBROUTINE NKCCERR(PHI)
C
	INTEGER SOLS,SUBS,VELS
C
	DOUBLE PRECISION
     1   V,C,S,KF,KB,KEQ
C
	DIMENSION
     1   S(10),KF(20),KB(20),KEQ(20),V(20),C(20)
C
	DOUBLE PRECISION PHI
	DIMENSION PHI(20)
C
        COMMON/NKCCKIN/ SOLS,SUBS,VELS,
     1   V,C,S,KF,KB,KEQ
C
C  SOLUTE INDEX:
C
C Variable Concentrations:
C	C1-	E1
C	C2-	E1-Na
C	C3-	E1-NaCl
C	C4-	E1-NaKCl
C	C5-	E1-NaKCl2
C	C6-	E2
C	C7-	E2-Cl
C	C8-	E2-KCl
C	C9-	E2-KCl2
C	C10-	E2-NaKCl2
C
C	C11-	E1-NaNH4Cl
C	C12-	E1-NaNH4Cl2
C	C13-	E2-NH4Cl
C	C14-	E2-NH4Cl2
C	C15-	E2-NaNH4Cl2
C
C Transported Solutes:
C	S1-	Na(1)
C	S2-	Na(2)
C	S3-	K(1)
C	S4-	K(2)
C	S5-	Cl(1)
C	S6-	Cl(2)
C	S7-	NH4(1)
C	S8-	NH4(2)
C
C Reactions:
C	V1-	C1 + S1 -> C2
C	V2-	C2 + S5 -> C3
C	V3-	C3 + S3 -> C4
C	V4-	C4 + S5 -> C5
C	V5-	C5 -> C10
C
C	V6-	C9 + S2 -> C10
C	V7-	C8 + S6 -> C9
C	V8-	C7 + S4 -> C8
C	V9-	C6 + S6 -> C7
C	V10-	C6 -> C1
C
C	V11-	C3 + S7 -> C11
C	V12-	C11 + S5 -> C12
C	V13-	C12 -> C15
C	V14-	C14 + S2 -> C15
C	V15-	C13 + S6 -> C14
C	V16-	C7 + S8 -> C13
C
C FLUXES
C
	V(1) = KEQ(1)*C(1)*S(1) - C(2)
	V(2) = KEQ(2)*C(2)*S(5) - C(3)
	V(3) = KEQ(3)*C(3)*S(3) - C(4)
	V(4) = KEQ(4)*C(4)*S(5) - C(5)
	V(5) = KF(5)*C(5) - KB(5)*C(10)
C
	V(6) = KEQ(6)*C(9)*S(2) - C(10)
	V(7) = KEQ(7)*C(8)*S(6) - C(9)
	V(8) = KEQ(8)*C(7)*S(4) - C(8)
	V(9) = KEQ(9)*C(6)*S(6) - C(7)
	V(10) = KF(10)*C(6) - KB(10)*C(1)
C
	V(11) = KEQ(11)*C(3)*S(7) - C(11)
	V(12) = KEQ(12)*C(11)*S(5) - C(12)
	V(13) = KF(13)*C(12) - KB(13)*C(15)
	V(14) = KEQ(14)*C(14)*S(2) - C(15)
	V(15) = KEQ(15)*C(13)*S(6) - C(14)
	V(16) = KEQ(16)*C(7)*S(8) - C(13)
C
C
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C  FIRST, CONSERVATION OF TOTAL CARRIER (SUM OF CI = 1.0)
C
	PHI(1)=-1.D0
	DO 110 I=1,SOLS
  110   PHI(1)=PHI(1) + C(I)
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
C
	PHI(2) = V(1)
	PHI(3) = V(2)
	PHI(4) = V(3)
	PHI(5) = V(4)
	PHI(6) = V(5) + V(13) - V(10)
C
	PHI(7) = V(6)
	PHI(8) = V(7)
	PHI(9) = V(8)
	PHI(10) = V(9)
C
	PHI(11) = V(11)
	PHI(12) = V(12)
	PHI(13) = V(14)
	PHI(14) = V(15)
	PHI(15) = V(16)
C
  170   CONTINUE
	RETURN
	END
	SUBROUTINE NKCCJAC(DPHIDC,M1)
C
	INTEGER SOLS,SUBS,VELS
C
	DOUBLE PRECISION
     1   V,C,S,KF,KB,KEQ
C
	DIMENSION
     1   S(10),KF(20),KB(20),KEQ(20),V(20),C(20)
C
	DOUBLE PRECISION DVDC,DPHIDC,DPHIDV
	DIMENSION DVDC(M1,M1),DPHIDC(M1,M1),DPHIDV(M1,M1)
C
        COMMON/NKCCKIN/ SOLS,SUBS,VELS,
     1   V,C,S,KF,KB,KEQ
C
	DO 103 J=1,M1
	DO 103 I=1,M1
  103   DVDC(I,J)=0.D0
C
	DO 105 J=1,M1
	DO 105 I=1,M1
  105   DPHIDC(I,J)=0.D0
C
	DO 107 J=1,M1
	DO 107 I=1,M1
  107   DPHIDV(I,J)=0.D0
C
C FLUX DERIVATIVES
C
C	V(1) = KEQ(1)*C(1)*S(1) - C(2)
C	V(2) = KEQ(2)*C(2)*S(5) - C(3)
C	V(3) = KEQ(3)*C(3)*S(3) - C(4)
C	V(4) = KEQ(4)*C(4)*S(5) - C(5)
C	V(5) = KF(5)*C(5) - KB(5)*C(10)
C
	DVDC(1,1) = KEQ(1)*S(1)
	DVDC(1,2) = - 1.0
C
	DVDC(2,2) = KEQ(2)*S(5)
	DVDC(2,3) = - 1.0
C
	DVDC(3,3) = KEQ(3)*S(3)
	DVDC(3,4) = - 1.0
C
	DVDC(4,4) = KEQ(4)*S(5)
	DVDC(4,5) = - 1.0
C
	DVDC(5,5) = KF(5)
	DVDC(5,10) = -KB(5) 
C
C	V(6) = KEQ(6)*C(9)*S(2) - C(10)
C	V(7) = KEQ(7)*C(8)*S(6) - C(9)
C	V(8) = KEQ(8)*C(7)*S(4) - C(8)
C	V(9) = KEQ(9)*C(6)*S(6) - C(7)
C	V(10) = KF(10)*C(6) - KB(10)*C(1)
C
	DVDC(6,9) = KEQ(6)*S(2)
	DVDC(6,10) = - 1.0
C
	DVDC(7,8) = KEQ(7)*S(6)
	DVDC(7,9) = - 1.0
C
	DVDC(8,7) = KEQ(8)*S(4)
	DVDC(8,8) = - 1.0
C
	DVDC(9,6) = KEQ(9)*S(6)
	DVDC(9,7) = - 1.0
C
	DVDC(10,6) = KF(10)
	DVDC(10,1) = -KB(10) 
C
C	V(11) = KEQ(11)*C(3)*S(7) - C(11)
C	V(12) = KEQ(12)*C(11)*S(5) - C(12)
C	V(13) = KF(13)*C(12) - KB(13)*C(15)
C	V(14) = KEQ(14)*C(14)*S(2) - C(15)
C	V(15) = KEQ(15)*C(13)*S(6) - C(14)
C	V(16) = KEQ(16)*C(7)*S(8) - C(13)
C
	DVDC(11,3) = KEQ(11)*S(7)
	DVDC(11,11) = - 1.0
C
	DVDC(12,11) = KEQ(12)*S(5)
	DVDC(12,12) = - 1.0
C
	DVDC(13,12) = KF(13)
	DVDC(13,15) = -KB(13) 
C
	DVDC(14,14) = KEQ(14)*S(2)
	DVDC(14,15) = - 1.0
C
	DVDC(15,13) = KEQ(15)*S(6)
	DVDC(15,14) = - 1.0
C
	DVDC(16,7) = KEQ(16)*S(8)
	DVDC(16,13) = - 1.0
C
C MODEL EQUATIONS ESTABLISH DPHIDV
C
C	PHI(1)=-1.D0
C	DO 110 I=1,SOLS
C  110   PHI(1)=PHI(1) + C(I)
C
C	PHI(2) = V(1)
C	PHI(3) = V(2)
C	PHI(4) = V(3)
C	PHI(5) = V(4)
C	PHI(6) = V(5) + V(13) - V(10)
C
	DPHIDV(2,1) = 1.D0
	DPHIDV(3,2) = 1.D0
	DPHIDV(4,3) = 1.D0
	DPHIDV(5,4) = 1.D0
	DPHIDV(6,5) = 1.D0
	DPHIDV(6,10) = -1.D0
	DPHIDV(6,13) = 1.D0
C
C
C	PHI(7) = V(6)
C	PHI(8) = V(7)
C	PHI(9) = V(8)
C	PHI(10) = V(9)
C
	DPHIDV(7,6) = 1.D0
	DPHIDV(8,7) = 1.D0
	DPHIDV(9,8) = 1.D0
	DPHIDV(10,9) = 1.D0
C
C
C	PHI(11) = V(11)
C	PHI(12) = V(12)
C	PHI(13) = V(14)
C	PHI(14) = V(15)
C	PHI(15) = V(16)
C
	DPHIDV(11,11) = 1.D0
	DPHIDV(12,12) = 1.D0
	DPHIDV(13,14) = 1.D0
	DPHIDV(14,15) = 1.D0
	DPHIDV(15,16) = 1.D0
C
C ESTABLISH THE ERROR VECTORS, THE "PHI" ARRAY.
C  FIRST, CONSERVATION OF TOTAL CARRIER (SUM OF CI = 1.0)
C
	DO 110 I=1,SOLS
  110   DPHIDC(1,I)=1.D0
C
C   MASS CONSERVATION IN THE STEADY-STATE CASE
C       DPHIDC = DPHIDV*DVDC
C
	DO 150 I=2,SOLS
	DO 150 K=1,SOLS
	DO 149 J=1,VELS
  149   DPHIDC(I,K)=DPHIDC(I,K) + DPHIDV(I,J)*DVDC(J,K)
  150	CONTINUE
C
C
	RETURN
	END
