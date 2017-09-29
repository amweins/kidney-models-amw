	PROGRAM MIX
C
C TO TAKE SIX INFLOWS WITH BUFFERED SOLUTES AND COMBINE TO A SINGLE OUTFLOW
C THE MIXING IS ASSUMED TO BE INSTANTANEOUS, WITH NO CO2 HYDRATION/DEHYDRATION
C OTHER BUFFERS ARE AT EQUILIBRIUM; CO2 EQUILIBRATION OCCURS ALONG THE TUBULE
C
	INTEGER SOLS,NTUB
C
C       SOLS-   NUMBER OF SOLUTES
C	NTUB-	NUMBER OF INPUT TUBULES
C       TUBULE INDICES: 1- SF, 2-6- JM, 7- OUTFLOW
C
C LUMINAL VARIABLES
	DOUBLE PRECISION PKC,
     1   CM(15,10),IMPM(10),LCHM(10),
     1   FVM(10),FKM(16,10)
C
C       CM-     LUMINAL SOLUTE CONCENTRATION
C       IMPM-   LUMINAL IMPERMEANT (INCL GLUCOSE)
C       LCHM-   LUMINAL PH
C	FVM-	LUMINAL VOLUME FLOW RATE
C	FKM-	LUMINAL SOLUTE FLUX RATE
C
	CHARACTER*5 SOL(12)
C
        COMMON SOLS,NTUB,
     1   CM,IMPM,LCHM,FVM,FKM
C
	OPEN (71,FILE='result.dat')
	OPEN (72,FILE='bound.dat')
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
C
C  SOLUTE INDEX:
C
	SOL(1)='  NA '
	SOL(2)='  K  '
	SOL(3)='  CL '
	SOL(4)=' HCO3'
	SOL(5)='H2CO3'
	SOL(6)=' CO2 '
	SOL(7)=' HPO4'
	SOL(8)='H2PO4'
	SOL(9)=' UREA'
	SOL(10)=' NH3 '
	SOL(11)=' NH4 '
	SOL(12)='  H  '
C
C
	PKC=3.57
	SOLS=12
	NTUB=6
C
C INPUT FLOWS AND CONCENTRATIONS:
	READ(72,50) 
     1   (FVM(J),J=1,NTUB),
     1   (IMPM(J),J=1,NTUB),
     1   ((CM(I,J),J=1,NTUB),I=1,SOLS-1)
   50   FORMAT (6D14.5,/,(6F14.9))
C
        DO 61 J=1,NTUB
	LCHM(J)=PKC + DLOG10(CM(4,J)/CM(5,J))
	CM(SOLS,J)=10.**(-LCHM(J))
        DO 60 I=1,SOLS
   60   FKM(I,J)=FVM(J)*CM(I,J)
        FKM(SOLS+1,J)=FVM(J)*IMPM(J)
   61   CONTINUE
C
	CALL MIXNEWT
C
        WRITE(71,90) (J, J=1,7)
        WRITE(71,91)
        DO 70 I=1,SOLS-1
   70   WRITE(71,92) SOL(I),(CM(I,J),J=1,7)
        WRITE(71,93) (LCHM(J),J=1,7)
        WRITE(71,94) (IMPM(J),J=1,7)
        WRITE(71,95) (J, J=1,7)
        WRITE(71,96) (FVM(J),J=1,7)
        WRITE(71,97)
        DO 75 I=1,SOLS
   75   WRITE(71,98) SOL(I),(FKM(I,J),J=1,7)
        WRITE(71,99) (FKM(SOLS+1,J),J=1,7)
C
   90   FORMAT(4X,'TUBE:',5X,7(6X,I2,6X))
   91   FORMAT(4X,' CM: ')
   92   FORMAT(4X,A5,5X,7F14.9)
   93   FORMAT(4X,' PH ',7X,7(4X,F6.4,4X))
   94   FORMAT(4X,'IMPM:',5X,7F14.9)
   95   FORMAT(//,4X,'TUBE:',5X,7(6X,I2,6X))
   96   FORMAT(4X,' FVM:',5X,7D14.5)
   97   FORMAT(4X,'FKM: ')
   98   FORMAT(4X,A5,5X,7D14.5)
   99   FORMAT(4X,'IMPM:',5X,7D14.5)
  100   CONTINUE
C
	STOP
	END