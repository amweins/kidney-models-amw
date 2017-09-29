C	TRNSIT
C  PROGRAM TO READ AN OUTPUT FILE FROM PROX.FOR
C  AND TO COMPUTE A TRANSIT TIME FROM THE FLOW PROFILE
C
C
	INTEGER IWR,JWR,CHOP
	REAL RES(63,12),TIME,TAU,X(21),VM(21),AM
        CHARACTER*12 IPTFIL
        CHARACTER*4 VNAME1,VNAME2
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A12)
	PRINT 8
    8   FORMAT(' TYPE THE SPATIAL CHOP:',$)
	READ *, CHOP
	CHOP=CHOP+1
C
	OPEN(UNIT=1,FILE=IPTFIL)
	OPEN(UNIT=2,FILE='trnsit.dat')
C
	AM=4.91
  506   IWR=1+2*CHOP
	JWR=3
C
C
	DO 810 KX=1,20
C
  300   READ (1,310,END=810) TIME
  310   FORMAT (100(/),100(/),100(/),48(/),11X,F15.4)
	READ (1,325) ((RES(I,J),J=1,12),I=1,CHOP)
  325   FORMAT (///,(3X,F8.6,11F11.6))
	READ (1,335) ((RES(I,J),J=1,12),I=1+CHOP,2*CHOP)
  335   FORMAT (///,(3X,F8.6,11F11.6))
	READ (1,345) (RES(I,1),(RES(I,J),J=3,12),I=1+2*CHOP,3*CHOP)
  345   FORMAT (///,(3X,F8.6,11X,10F11.6))
C
C
	DO 399 I=1,CHOP
	X(I)=RES(I,1)
  399   VM(I)=RES(IWR-1+I,JWR)/AM
	TAU=-.5/VM(1) - .5/VM(CHOP)
	DO 400 I=1,CHOP
  400   TAU=TAU + 1./VM(I)
	TAU=TAU*(X(2)-X(1))
C
	WRITE (2,415) TAU,VM(1),VM(CHOP)
  415   FORMAT (3E16.8)
C
  810   CONTINUE
C
	STOP
	END