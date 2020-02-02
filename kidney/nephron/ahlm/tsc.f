	SUBROUTINE TSC(V,S)
C
	DOUBLE PRECISION
     1   V,S,KF,KB,KEQ
	DIMENSION
     1   S(4),KF(10),KB(10),KEQ(10)
C
C
	DOUBLE PRECISION
     1   STSC(2,2),PTSC(2,2),KTSC(4,2),CTSC(4,2)
C
C       S-      SUBSTRATE CONCENTRATION (MMOL/ML)
C	P-	TRANSLOCATION RATE CONSTANT
C	KEQ-	EQUILIBRIUM CONSTANT (MMOL/ML)
C       C-	UNKNOWN CONCENTRATION 
C       V-	NET FORWARD REACTION VELOCITY
C
C THE FIRST INDEX OF (I,J) IS THE SOLUTE 
C  THE SECOND IS THE SIDE (J=1, EXTERNAL; J=2, INTERNAL)
C
	DOUBLE PRECISION SK(3,2),DTSC(2),ETSC(2),SIG
C
	DATA KEQ/
     1     .3416D+04,  .8870D+01,  .1771D+04,  .6819D+06,
     1     .3416D+04,  .8870D+01,  .1771D+04,  .6819D+06/
	KF(9)= 0.4295D+07
	KB(9)= 0.1000D+06
	KF(10)= 0.7692D+04
	KB(10)= 0.1791D+03
C
C Original TSC coefficients
C	DATA KEQ/
C     1     .9263D+03,  .2585D+04,  .1999D+03,  .7165D+02,
C     1     .9263D+03,  .2585D+04,  .1999D+03,  .7165D+02/
C	KF(9)= 0.9597D+06
C	KB(9)= 0.1000D+06
C	KF(10)= 0.1139D+04
C	KB(10)= 0.1187D+03
C
C Variable Concentrations
C	CTSC(1,1)-	E1
C	CTSC(2,1)-	E1Na
C	CTSC(3,1)-	E1Cl
C	CTSC(4,1)-	E1NaCl
C	CTSC(1,2)-	E2
C	CTSC(2,2)-	E2Na
C	CTSC(3,2)-	E2Cl
C	CTSC(4,2)-	E2NaCl
C
C Specified Solutes - S
C	S1-	STSC(1,1)-	Na(1)
C	S2-	STSC(2,1)-	Cl(1)
C	S3-	STSC(1,2)-	Na(2)
C	S4-	STSC(2,2)-	Cl(2)
C
C
C Translate into the TSC variables.
C
	KTSC(1,1)=1.D0/KEQ(1)
	KTSC(2,1)=1.D0/KEQ(2)
	KTSC(3,1)=1.D0/KEQ(3)
	KTSC(4,1)=1.D0/KEQ(4)
	KTSC(1,2)=1.D0/KEQ(5)
	KTSC(2,2)=1.D0/KEQ(6)
	KTSC(3,2)=1.D0/KEQ(7)
	KTSC(4,2)=1.D0/KEQ(8)
C
	PTSC(1,1)=KF(9)
	PTSC(2,1)=KF(10)
	PTSC(1,2)=KB(9)
	PTSC(2,2)=KB(10)
C
	STSC(1,1)=S(1)
	STSC(2,1)=S(2)
	STSC(1,2)=S(3)
	STSC(2,2)=S(4)
C
C Solve the TSC equations
C
	SK(1,1)=STSC(1,1)/KTSC(1,1)
	SK(2,1)=STSC(2,1)/KTSC(2,1)
	SK(3,1)=STSC(1,1)*STSC(2,1)/(KTSC(3,1)*KTSC(1,1))
	SK(1,2)=STSC(1,2)/KTSC(1,2)
	SK(2,2)=STSC(2,2)/KTSC(2,2)
	SK(3,2)=STSC(1,2)*STSC(2,2)/(KTSC(3,2)*KTSC(1,2))
C	
	DTSC(1)=1.D0 + SK(1,1) + SK(2,1) + SK(3,1)
	DTSC(2)=1.D0 + SK(1,2) + SK(2,2) + SK(3,2)
	ETSC(1)=PTSC(1,1) + SK(3,1)*PTSC(2,1)
	ETSC(2)=PTSC(1,2) + SK(3,2)*PTSC(2,2)
	SIG=DTSC(1)*ETSC(2) + DTSC(2)*ETSC(1)
C
	CTSC(1,1)=ETSC(2)/SIG
	CTSC(2,1)=SK(1,1)*CTSC(1,1)
	CTSC(3,1)=SK(2,1)*CTSC(1,1)
	CTSC(4,1)=SK(3,1)*CTSC(1,1)
	CTSC(1,2)=ETSC(1)/SIG
	CTSC(2,2)=SK(1,2)*CTSC(1,2)
	CTSC(3,2)=SK(2,2)*CTSC(1,2)
	CTSC(4,2)=SK(3,2)*CTSC(1,2)
C
	V = PTSC(2,1)*CTSC(4,1) - PTSC(2,2)*CTSC(4,2)
C
	RETURN
	END
