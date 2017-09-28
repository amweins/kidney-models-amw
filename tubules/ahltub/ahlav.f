C	DCTAV
C  PROGRAM TO READ AN OUTPUT FILE FROM PROX.FOR
C  AND AVERAGE THE INDICATED COLUMNS
C
C
	INTEGER NWR,IWR(20),JWR(20),CHOP
	REAL RES(243,16),TIME,AVG(40,10)
        CHARACTER*12 IPTFIL,OPTFIL
        CHARACTER*5 VNAME1,VNAME2,SWI
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A12)
	PRINT 6
    6   FORMAT(' TYPE THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
	PRINT 8
    8   FORMAT(' TYPE THE SPATIAL CHOP:',$)
	READ *, CHOP
	CHOP=CHOP+1
C
	OPEN(UNIT=1,FILE=IPTFIL)
	OPEN(UNIT=2,FILE=OPTFIL)
C
	PRINT 500
  500   FORMAT(' IS INPUT FROM A DATA FILE OR INTERACTIVE? (d/i): ',$)
	READ 30, SWI
  510   IF (SWI.EQ.'i') GO TO 9
	OPEN(UNIT=3,FILE='ptwr.dat')
	READ (3,515) NWR
  515   FORMAT(3X,I2)
	DO 700 KWR=1,NWR
	READ (3,30) VNAME1
C
	IF(VNAME1.EQ.'vm') GO TO 501
	IF(VNAME1.EQ.'pm') GO TO 502
	IF(VNAME1.EQ.'cm') GO TO 503
	IF(VNAME1.EQ.'phm') GO TO 504
	IF(VNAME1.EQ.'xm') GO TO 505
	IF(VNAME1.EQ.'fvm') GO TO 506
	IF(VNAME1.EQ.'fkm') GO TO 507
	IF(VNAME1.EQ.'osmm') GO TO 508
C
  501   IWR(KWR)=1
	JWR(KWR)=2
	GO TO 700
  502   IWR(KWR)=1
	JWR(KWR)=3
	GO TO 700
  503   IWR(KWR)=1
	JWR(KWR)=4
	GO TO 600
  508   IWR(KWR)=1+CHOP
	JWR(KWR)=2
	GO TO 700
  504   IWR(KWR)=1+CHOP
	JWR(KWR)=3
	GO TO 700
  505   IWR(KWR)=1+CHOP
	JWR(KWR)=4
	GO TO 600
  506   IWR(KWR)=1+2*CHOP
	JWR(KWR)=3
	GO TO 700
  507   IWR(KWR)=1+2*CHOP
	JWR(KWR)=4
	GO TO 600
C
C
  600   READ (3,30) VNAME2
C
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
C
C
  621   JWR(KWR)=JWR(KWR)+0
	GO TO 700
  622   JWR(KWR)=JWR(KWR)+1
	GO TO 700
  623   JWR(KWR)=JWR(KWR)+2
	GO TO 700
  624   JWR(KWR)=JWR(KWR)+3
	GO TO 700
  625   JWR(KWR)=JWR(KWR)+4
	GO TO 700
  626   JWR(KWR)=JWR(KWR)+5
	GO TO 700
  627   JWR(KWR)=JWR(KWR)+6
	GO TO 700
  628   JWR(KWR)=JWR(KWR)+7
	GO TO 700
  629   JWR(KWR)=JWR(KWR)+8
	GO TO 700
  630   JWR(KWR)=JWR(KWR)+9
	GO TO 700
  631   JWR(KWR)=JWR(KWR)+10
	GO TO 700
  632   JWR(KWR)=JWR(KWR)+11
	GO TO 700
C
C
  700   CONTINUE
	GO TO 801
C
C
    9   PRINT 10
   10   FORMAT(' TYPE THE NUMBER OF VARIABLES TO BE REWRITTEN')
	READ *, NWR
	DO 800 KWR=1,NWR
	PRINT 20, KWR
   20   FORMAT(' TYPE THE NAME OF VARIABLE NUMBER-',I2)
	READ 30, VNAME1
   30   FORMAT(A4)
C
C
	IF(VNAME1.EQ.'vm') GO TO 101
	IF(VNAME1.EQ.'pm') GO TO 102
	IF(VNAME1.EQ.'cm') GO TO 103
	IF(VNAME1.EQ.'phm') GO TO 104
	IF(VNAME1.EQ.'xm') GO TO 105
	IF(VNAME1.EQ.'fvm') GO TO 106
	IF(VNAME1.EQ.'fkm') GO TO 107
	IF(VNAME1.EQ.'osmm') GO TO 108
C
C
  101   IWR(KWR)=1
	JWR(KWR)=2
	GO TO 800
  102   IWR(KWR)=1
	JWR(KWR)=3
	GO TO 800
  103   IWR(KWR)=1
	JWR(KWR)=4
	GO TO 200
  108   IWR(KWR)=1+CHOP
	JWR(KWR)=2
	GO TO 800
  104   IWR(KWR)=1+CHOP
	JWR(KWR)=3
	GO TO 800
  105   IWR(KWR)=1+CHOP
	JWR(KWR)=4
	GO TO 200
  106   IWR(KWR)=1+2*CHOP
	JWR(KWR)=3
	GO TO 800
  107   IWR(KWR)=1+2*CHOP
	JWR(KWR)=4
	GO TO 200
C
C
  200   PRINT 201
  201   FORMAT(' TYPE THE SPECIES')
	READ 30, VNAME2
C
C
	IF(VNAME2.EQ.'na') GO TO 221
	IF(VNAME2.EQ.'k') GO TO 222
	IF(VNAME2.EQ.'cl') GO TO 223
	IF(VNAME2.EQ.'hco3') GO TO 224
	IF(VNAME2.EQ.'h2co3') GO TO 225
	IF(VNAME2.EQ.'co2') GO TO 226
	IF(VNAME2.EQ.'hpo4') GO TO 227
	IF(VNAME2.EQ.'h2po4') GO TO 228
	IF(VNAME2.EQ.'urea') GO TO 229
	IF(VNAME2.EQ.'nh3') GO TO 230
	IF(VNAME2.EQ.'nh4') GO TO 231
	IF(VNAME2.EQ.'h') GO TO 232
	IF(VNAME2.EQ.'impm') GO TO 232
C
C
  221   JWR(KWR)=JWR(KWR)+0
	GO TO 800
  222   JWR(KWR)=JWR(KWR)+1
	GO TO 800
  223   JWR(KWR)=JWR(KWR)+2
	GO TO 800
  224   JWR(KWR)=JWR(KWR)+3
	GO TO 800
  225   JWR(KWR)=JWR(KWR)+4
	GO TO 800
  226   JWR(KWR)=JWR(KWR)+5
	GO TO 800
  227   JWR(KWR)=JWR(KWR)+6
	GO TO 800
  228   JWR(KWR)=JWR(KWR)+7
	GO TO 800
  229   JWR(KWR)=JWR(KWR)+8
	GO TO 800
  230   JWR(KWR)=JWR(KWR)+9
	GO TO 800
  231   JWR(KWR)=JWR(KWR)+10
	GO TO 800
  232   JWR(KWR)=JWR(KWR)+11
	GO TO 800
  233   JWR(KWR)=JWR(KWR)+12
	GO TO 800
C
C
  800   CONTINUE
  801   CONTINUE
C
	DO 810 KX=1,100
C
  300   READ (1,310,END=420) TIME
  310   FORMAT (421(/),11X,F15.4)
	READ (1,325) ((RES(I,J),J=1,15),I=1,CHOP)
  325   FORMAT (///,(F6.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,D9.3,2F9.6))
	READ (1,335) ((RES(I,J),J=1,15),I=1+CHOP,2*CHOP)
  335   FORMAT (///,(F6.4,14F9.4))
	READ (1,345) (RES(I,1),(RES(I,J),J=3,15),I=1+2*CHOP,3*CHOP)
  345   FORMAT (///,(F6.4,9X,5F9.6,D9.3,4F9.6,D9.3,2F9.6))
C
	DO 401 J=1,NWR
	AVG(KX,J)=.5*RES(IWR(J),JWR(J))
	DO 400 I=2,CHOP-1
  400   AVG(KX,J)=AVG(KX,J)+RES(IWR(J)+I-1,JWR(J))
	AVG(KX,J)=AVG(KX,J)+.5*RES(IWR(J)+CHOP-1,JWR(J))
  401   AVG(KX,J)=AVG(KX,J)/FLOAT(CHOP-1)
C
  810   CONTINUE
C
  420   DO 410 I=1,KX-1
  410   WRITE (2,415) (AVG(I,J),J=1,NWR)
  415   FORMAT (8E16.8)
C
	STOP
	END
