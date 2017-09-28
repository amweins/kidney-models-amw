	PROGRAM REFLUX
C READS FLUXES.DAT AND OUTPUTS THE REQUESTED SOLUTE IN GRAPHABLE FORMAT
C
	INTEGER CHOP,DCHOP,ISOL
	DOUBLE PRECISION MNA(405,9),MK(405,9),MCL(405,8),MNH4(405,12)
        CHARACTER*15 IPTFIL,OPTFIL,SOL,SWI
C
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A15)
	OPEN(1,FILE=IPTFIL)
C
	PRINT 3
    3   FORMAT(' MESH SPACING - DATA FILE OR INTERACTIVE? (d/i):',$)
	READ *, SWI
	IF (SWI.EQ.'i') GO TO 7
	OPEN(4,FILE='ahlmesh.dat')
	READ (4,5) CHOP,DCHOP
    5   FORMAT(2I5)
	GO TO 12
C
    7   PRINT 8  
    8   FORMAT(' SPATIAL CHOP, AND PRINT SPACING:',$) 
	READ *, CHOP,DCHOP
   10   CONTINUE
   12   CHOP=CHOP+1
C
	DO 910 KX=1,5
C
	I1=(KX-1)*CHOP+1
	I2=(KX)*CHOP
C
	READ (1,320,END=440) ((MNA(I,J),J=1,9),I=I1,I2)
	READ (1,350) ((MK(I,J),J=1,9),I=I1,I2)
	READ (1,380) ((MCL(I,J),J=1,8),I=I1,I2)
	READ (1,410) ((MNH4(I,J),J=1,12),I=I1,I2)
  320   FORMAT (///,(F7.4,2X,8F10.4))
  350   FORMAT (///,(F7.4,2X,8F10.4))
  380   FORMAT (///,(F7.4,2X,7F10.4))
  410   FORMAT (///,(F7.4,2X,11F10.4))
C
  910   CONTINUE
C
  440   OPEN(2,FILE='xyf.na')
  510   DO 515 I=1,(KX-1)*CHOP,DCHOP
  515   WRITE(2,388) (MNA(I,J),J=1,9)
	CLOSE (2)
C
	OPEN(2,FILE='xyf.k')
  520   DO 525 I=1,(KX-1)*CHOP,DCHOP
  525   WRITE(2,388) (MK(I,J),J=1,9)
	CLOSE (2)
C	
	OPEN(2,FILE='xyf.cl')
  530   DO 535 I=1,(KX-1)*CHOP,DCHOP
  535   WRITE(2,388) (MCL(I,J),J=1,8)
	CLOSE (2)
C
	OPEN(2,FILE='xyf.nh4')
  540   DO 545 I=1,(KX-1)*CHOP,DCHOP
  545   WRITE(2,388) (MNH4(I,J),J=1,12)
	CLOSE (2)
C
	OPEN(2,FILE='xyf.cell')
  550   DO 555 I=1,(KX-1)*CHOP,DCHOP
  555   WRITE(2,388) MNA(I,1),
     1   (MNA(I,J),J=4,9),
     1   (MCL(I,J),J=4,8),
     1   (MK(I,J),J=4,9),
     1   (MNH4(I,J),J=4,12)
	CLOSE (2)
C
  388   FORMAT(8E16.8)
C
	END
C
