	PROGRAM DCTWR
C  PROGRAM TO READ AN OUTPUT FILE FROM DCT.F
C
C
	INTEGER NWR,IWR(20),JWR(20),CHOP,DCHOP
	REAL WRRES(3601,200),RES(3601,20),TIME(201)
        CHARACTER*12 IPTFIL,OPTFIL,VARNAM
        CHARACTER*5 VNAME1(20),VNAME2(20),SWI,SWT
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A12)
	PRINT 6
    6   FORMAT(' TYPE THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
C
	PRINT 3
    3   FORMAT(' MESH SPACING - DATA FILE OR INTERACTIVE? (d/i):',$)
	READ *, SWI
	IF (SWI.EQ.'i') GO TO 7
	OPEN(4,FILE='dctmesh.dat')
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
	OPEN(1,FILE=IPTFIL)
	OPEN(2,FILE=OPTFIL)
	IWR(1)=1
	JWR(1)=1
C
C
	PRINT 500
  500   FORMAT(' IS INPUT FROM A DATA FILE OR INTERACTIVE? (d/i): ',$)
	READ 30, SWI
  510   IF (SWI.EQ.'i') GO TO 15
	PRINT 516
  516   FORMAT(' TYPE THE VARIABLE FILENAME:',$)
	READ 4, VARNAM
	OPEN(3,FILE=VARNAM)
	READ (3,515) NWR
  515   FORMAT(3X,I2)
	DO 525 KWR=2,NWR+1
  525   READ (3,530) VNAME1(KWR),VNAME2(KWR)
  530   FORMAT(2A5)
	GO TO 801
C
C
   15   PRINT 16
   16   FORMAT(' TYPE THE NUMBER OF VARIABLES TO BE REWRITTEN')
	READ *, NWR
	DO 50 KWR=2,NWR+1
	PRINT 20, KWR-1
   20   FORMAT(' TYPE THE NAME OF VARIABLE NUMBER-',I2)
	READ 30, VNAME1(KWR)
   30   FORMAT(A5)
	IF(VNAME1(KWR).EQ.'cm') GO TO 40
	IF(VNAME1(KWR).EQ.'xm') GO TO 40
	IF(VNAME1(KWR).EQ.'fkm') GO TO 40
	GO TO 50
C
   40   PRINT 45
   45   FORMAT(' TYPE THE SPECIES')
	READ 30, VNAME2(KWR)
   50   CONTINUE
C
  801   DO 810 KWR=1,NWR+1
  810   CALL SELEC 
     1   (VNAME1(KWR),VNAME2(KWR),IWR(KWR),JWR(KWR),CHOP)
C
	PRINT 802
  802   FORMAT(' STEADY STATE (ss) OR TRANSIENT (tr) EXPT.? ',$)
	READ 30, SWT
	IF(SWT.EQ.'tr') GO TO 900
C
	DO 405 KX=1,150
C
  300   READ (1,310,END=420) TIME(KX)
  310   FORMAT (2(/),11X,F15.4)
	READ (1,325) ((RES(I,J),J=1,15),I=1,CHOP)
  325   FORMAT (///,(F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,D9.3,2F9.6))
	READ (1,335) ((RES(I,J),J=1,15),I=1+CHOP,2*CHOP)
  335   FORMAT (///,(F7.4,14F9.4))
	READ (1,345) ((RES(I,J),J=1,15),I=1+2*CHOP,3*CHOP)
  345   FORMAT (///,(F7.4,6F9.6,D9.3,4F9.6,D9.3,2F9.6))
C
	DO 399 I=1,CHOP
  399   WRRES(I,1)=RES(I,1)
	DO 400 J=2,NWR+1
	DO 400 I=1,CHOP
  400   WRRES(I,J+(KX-1)*NWR)=RES(IWR(J)+I-1,JWR(J))
C
  405   CONTINUE
C
  420   DO 410 I=1,CHOP,DCHOP
  410   WRITE (2,415) (WRRES(I,J),J=1,(KX-1)*NWR+1)
  415   FORMAT (8E16.8)
	STOP
C
C
  900   CONTINUE
C
	KX=1
	READ (1,310,END=440) TIME(KX)
 	READ (1,325) ((RES(I,J),J=1,16),I=1,CHOP)
	READ (1,335) ((RES(I,J),J=1,16),I=1+CHOP,2*CHOP)
	READ (1,345) (RES(I,1),(RES(I,J),J=3,16),I=1+2*CHOP,3*CHOP)
C
	DO 358 I=1,CHOP
  358   WRRES(I,1)=RES(I,1)
	DO 359 J=2,NWR+1
	DO 359 I=1,CHOP
  359   WRRES(I,J+(KX-1)*NWR)=RES(IWR(J)+I-1,JWR(J))
C
C
	DO 910 KX=2,150
C
	READ (1,360,END=440) TIME(KX)
  360   FORMAT (100(/),61(/),11X,F15.4)
	READ (1,325) ((RES(I,J),J=1,16),I=1,CHOP)
	READ (1,335) ((RES(I,J),J=1,16),I=1+CHOP,2*CHOP)
	READ (1,345) ((RES(I,J),J=1,16),I=1+2*CHOP,3*CHOP)
C
C
	DO 370 J=2,NWR+1
	DO 370 I=1,CHOP
  370   WRRES(I,J+(KX-1)*NWR)=RES(IWR(J)+I-1,JWR(J))
C
  910   CONTINUE
C
  440   WRITE(2,386) (WRRES(I,1),I=1,CHOP,DCHOP)
  386   FORMAT(16X,7E16.8,/,(8E16.8))
	DO 387 J=2,1+(KX-1)*NWR
  387   WRITE(2,388) TIME(J-1),(WRRES(I,J),I=1,CHOP,DCHOP)
  388   FORMAT(8E16.8)
C
C
	STOP
	END
C
C
	SUBROUTINE SELEC(VNAME1,VNAME2,IWR,JWR,CHOP)
C PROGRAM TO ASSIGN THE ARRAY COORDINATES FOR THE OUTPUT VARIABLES
C
	INTEGER IWR,JWR,CHOP
        CHARACTER*5 VNAME1,VNAME2
C
C
	IF(VNAME1.EQ.'vm') GO TO 501
	IF(VNAME1.EQ.'pm') GO TO 502
	IF(VNAME1.EQ.'cm') GO TO 503
	IF(VNAME1.EQ.'phm') GO TO 504
	IF(VNAME1.EQ.'xm') GO TO 505
	IF(VNAME1.EQ.'fvm') GO TO 506
	IF(VNAME1.EQ.'fkm') GO TO 507
C	IF(VNAME1.EQ.'pheq') GO TO 508
	IF(VNAME1.EQ.'ta') GO TO 508
	IF(VNAME1.EQ.'osmm') GO TO 509
C
C
  501   IWR=1
	JWR=2
	GO TO 700
  502   IWR=1
	JWR=3
	GO TO 700
  503   IWR=1
	JWR=4
	GO TO 600
  508   IWR=1+CHOP
	JWR=2
	GO TO 700
  504   IWR=1+CHOP
	JWR=3
	GO TO 700
  505   IWR=1+CHOP
	JWR=4
	GO TO 600
  506   IWR=1+2*CHOP
	JWR=3
	GO TO 700
  507   IWR=1+2*CHOP
	JWR=4
	GO TO 600
  509   IWR=1+2*CHOP
	JWR=2
	GO TO 700
C
C
  600   CONTINUE
C
	IF(VNAME2.EQ.'na') GO TO 621
	IF(VNAME2.EQ.'k') GO TO 622
	IF(VNAME2.EQ.'cl') GO TO 623
	IF(VNAME2.EQ.'hco3') GO TO 624
	IF(VNAME2.EQ.'h2co3') GO TO 625
	IF(VNAME2.EQ.'co2') GO TO 626
	IF(VNAME2.EQ.'hpo4') GO TO 627
	IF(VNAME2.EQ.'h2po4') GO TO 628
	IF(VNAME2.EQ.'urea') GO TO 629
	IF(VNAME2.EQ.'nh3') GO TO 630
	IF(VNAME2.EQ.'nh4') GO TO 631
	IF(VNAME2.EQ.'h') GO TO 632
	IF(VNAME2.EQ.'impm') GO TO 632
	IF(VNAME2.EQ.'ta') GO TO 632
C
C
  621   JWR=JWR+0
	GO TO 700
  622   JWR=JWR+1
	GO TO 700
  623   JWR=JWR+2
	GO TO 700
  624   JWR=JWR+3
	GO TO 700
  625   JWR=JWR+4
	GO TO 700
  626   JWR=JWR+5
	GO TO 700
  627   JWR=JWR+6
	GO TO 700
  628   JWR=JWR+7
	GO TO 700
  629   JWR=JWR+8
	GO TO 700
  630   JWR=JWR+9
	GO TO 700
  631   JWR=JWR+10
	GO TO 700
  632   JWR=JWR+11
	GO TO 700
C
  700   CONTINUE
	RETURN
	END
