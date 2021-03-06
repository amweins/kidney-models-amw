	PROGRAM MMULT
C TO READ TWO MATRICES AND STORE THEIR PRODUCT
	DOUBLE PRECISION M1(100,100),M2(100,100),M3(100,100)
	CHARACTER*13 IPTFL1,IPTFL2,OPTFIL
C
	PRINT 2
    2   FORMAT(/,' PRINT THE INPUT FILENAME FOR MATRIX #1:',$)
	READ 4, IPTFL1
    4   FORMAT (A13)
	OPEN(1,FILE=IPTFL1)
C
	PRINT 5, IPTFL1
    5   FORMAT(/,' PRINT THE NUMBER OF COLUMNS IN ',A13,$)
	READ *, NCOLM1
C
	PRINT 12
   12   FORMAT(/,' PRINT THE INPUT FILENAME FOR MATRIX #2:',$)
	READ 14, IPTFL2
   14   FORMAT (A13)
	OPEN(2,FILE=IPTFL2)
C
	PRINT 15, IPTFL2
   15   FORMAT(/,' PRINT THE NUMBER OF COLUMNS IN ',A13,$)
	READ *, NCOLM2
C
	PRINT 6
    6   FORMAT(/,' PRINT THE OUTPUT FILENAME:',$)
	READ 8, OPTFIL
    8   FORMAT (A13)
	OPEN(3,FILE=OPTFIL)
C
	DO 20 I=1,100
	READ(1,30,END=22) (M1(I,J),J=1,NCOLM1)
   20   CONTINUE
   22   NROWM1=I-1
	DO 25 I=1,100
	READ(2,30,END=26) (M2(I,J),J=1,NCOLM2)
   25   CONTINUE
   26   NROWM2=I-1
   30   FORMAT(8D16.8)
C
C CHECK THAT THE MATRICES ARE COMPATIBLE
	IF (NCOLM1.EQ.NROWM2) GO TO 40
	PRINT 35
   35   FORMAT('INCOMPATIBLE MATRICES FOR MULTIPLICATION')
	STOP
C
   40   DO 45 I=1,NROWM1
	DO 45 J=1,NCOLM2
	M3(I,J)=0.D0
	DO 45 K=1,NCOLM1
	M3(I,J)=M3(I,J) + M1(I,K)*M2(K,J)
   45   CONTINUE
C
	DO 47 I=1,NROWM1
   47   WRITE (3,48) (M3(I,J),J=1,NCOLM2)
   48   FORMAT(8D16.8)
C
	STOP
	END
