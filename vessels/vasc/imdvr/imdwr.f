	PROGRAM IMDWR
C  PROGRAM TO READ AN OUTPUT FILE FROM CAP.FOR
C
C
	INTEGER NWR,IWR(20),JWR(20),IC
	REAL RES(901,15),WRRES(101,80),TIME(80)
        CHARACTER*16 IPTFIL,OPTFIL
        CHARACTER*8 VNAME1(20),VNAME2(20),SWI,SWSOU
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A16)
	PRINT 6
    6   FORMAT(' TYPE THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
	PRINT 8
    8   FORMAT(' TYPE THE SPATIAL CHOP:',$)
	READ *, IC
	IC=IC+1
C
	OPEN(UNIT=1,FILE=IPTFIL,STATUS='old')
	OPEN(UNIT=2,FILE=OPTFIL)
	IWR(1)=1
	JWR(1)=1
C
	PRINT 490
  490   FORMAT(' IS INPUT FROM A DATA FILE OR INTERACTIVE? (d/i): ',$)
	READ *, SWI
  495   IF (SWI.EQ.'i') GO TO 15
	OPEN(UNIT=3,FILE='imdwr.dat',STATUS='old')
	READ (3,498) NWR
  498   FORMAT(3X,I2)
	DO 500 KWR=2,NWR+1
        READ (3,30) VNAME1(KWR),VNAME2(KWR)
   30   FORMAT(2A8)
        CALL SELEC 
     1   (VNAME1(KWR),VNAME2(KWR),IWR(KWR),JWR(KWR),KWR,IC,'fil')
  500   CONTINUE
	GO TO 801
C
C
   15   PRINT 16
   16   FORMAT(' TYPE THE NUMBER OF VARIABLES TO BE REWRITTEN')
	READ *, NWR
        SWSOU='int'
	DO 50 KWR=2,NWR+1
        CALL SELEC 
     1   (VNAME1(KWR),VNAME2(KWR),IWR(KWR),JWR(KWR),KWR,IC,SWSOU)
   50   CONTINUE
C
C THE WRITTEN VARIABLES HAVE BEEN SPECIFIED
C
  801   DO 810 KX=1,150
C
  300   READ (1,310,END=400) TIME(KX)
  310   FORMAT (38(/),11X,F15.4)
	READ (1,320) ((RES(I,J),J=1,13),I=1,IC)
	READ (1,325) ((RES(I,J),J=1,13),I=1+IC,2*IC)
	READ (1,325) ((RES(I,J),J=1,13),I=1+2*IC,3*IC)
C
	READ (1,330) ((RES(I,J),J=1,14),I=1+3*IC,4*IC)
	READ (1,335) (RES(1+4*IC,J),J=2,14)
	READ (1,340) ((RES(I,J),J=1,13),I=2+4*IC,1+5*IC)
	READ (1,345) 
     1   ((RES(I,J),J=1,9),RES(I,11),RES(I,12),I=2+5*IC,1+6*IC)
	READ (1,350) 
     1   ((RES(I,J),J=1,9),RES(I,11),RES(I,12),I=2+6*IC,1+7*IC)
  320   FORMAT (///,(3X,F8.6,12F11.6))
  325   FORMAT (///,(3X,F8.6,12F11.6))
  330   FORMAT (///,(3X,F8.6,13F11.6))
  335   FORMAT (//,(11X,13F11.6))
  340   FORMAT (///,(3X,F8.6,12D11.3))
  345   FORMAT (6(/),(3X,F8.6,F10.6,D12.4,6F11.6,11X,2F11.6))
  350   FORMAT (///,(3X,F8.6,8F11.6,11X,2F11.6))
C
C
  460   DO 462 I=1,IC
  462   WRRES(I,1)=RES(I,1)
	DO 464 J=2,NWR+1
	DO 464 I=1,IC
  464   WRRES(I,J+(KX-1)*NWR)=RES(IWR(J)+I-1,JWR(J))
C
  810   CONTINUE
C
  400   DO 412 I=1,IC
  412   WRITE (2,415) (WRRES(I,J),J=1,(KX-1)*NWR+1)
  415   FORMAT (8E16.8)
C
	STOP
	END
C
        SUBROUTINE SELEC (VNAME1,VNAME2,IWR,JWR,KWR,IC,SWSOU)
C
	INTEGER NWR,IWR,JWR,KWR,IC
        CHARACTER*8 VNAME1,VNAME2,SWSOU
C
   18   IF(SWSOU.EQ.'int') THEN
        PRINT 20, KWR-1
   20   FORMAT(' TYPE THE NAME OF VARIABLE NUMBER-',I2)
	READ 30, VNAME1
   30   FORMAT(A8)
        ELSE
        ENDIF
C
	IF(VNAME1.EQ.'vc') THEN 
        IWR=1+2*IC
        JWR=2
        GO TO 800
        ELSE IF(VNAME1.EQ.'pc') THEN 
        IWR=1+IC
        JWR=2
        GO TO 800
        ELSE IF(VNAME1.EQ.'cc') THEN 
        IWR=1+IC
        JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'phc') THEN 
        IWR=1+2*IC
        JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'xc') THEN 
        IWR=1+2*IC
        JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'fvc') THEN 
        IWR=1+3*IC
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'fkc') THEN 
        IWR=1+3*IC
	JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'fbc') THEN 
        IWR=1+3*IC
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'hct') THEN 
        IWR=1+3*IC
	JWR=14
        GO TO 800
C
	ELSE IF(VNAME1.EQ.'ps') THEN 
        IWR=1
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'cs') THEN 
        IWR=1
	JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'phs') THEN 
        IWR=1
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'fvcs') THEN 
        IWR=2+4*IC
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'fkcs') THEN 
        IWR=2+4*IC
	JWR=3
        GO TO 200
C
	ELSE IF(VNAME1.EQ.'apr') THEN 
        IWR=1+4*IC
	JWR=2
        GO TO 200
	ELSE IF(VNAME1.EQ.'cur') THEN 
        IWR=2+4*IC
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'h') THEN 
        IWR=2+5*IC
	JWR=3
        GO TO 800
	ELSE IF(VNAME1.EQ.'hgb') THEN 
        GO TO 250
        ELSE
        PRINT 40
   40   FORMAT(' NO SUCH VARIABLE NAME, TRY AGAIN')
        IF(SWSOU.EQ.'int') THEN
        GO TO 18
        ELSE
        STOP
        ENDIF
        ENDIF
C
C
  200   IF(SWSOU.EQ.'int') THEN
        PRINT 201
  201   FORMAT(' TYPE THE SPECIES')
	READ 30, VNAME2
        ELSE
        ENDIF
C
	IF(VNAME2.EQ.'vol') THEN 
        JWR=JWR-1
	GO TO 800
	ELSE IF(VNAME2.EQ.'na') THEN 
        JWR=JWR+0
	GO TO 800
	ELSE IF(VNAME2.EQ.'k') THEN 
        JWR=JWR+1
	GO TO 800
	ELSE IF(VNAME2.EQ.'cl') THEN 
        JWR=JWR+2
	GO TO 800
	ELSE IF(VNAME2.EQ.'hco3') THEN 
        JWR=JWR+3
	GO TO 800
	ELSE IF(VNAME2.EQ.'hpo') THEN 
        JWR=JWR+4
	GO TO 800
	ELSE IF(VNAME2.EQ.'h2po') THEN 
        JWR=JWR+5
	GO TO 800
	ELSE IF(VNAME2.EQ.'urea') THEN 
        JWR=JWR+6
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh4') THEN 
        JWR=JWR+7
	GO TO 800
	ELSE IF(VNAME2.EQ.'hco2') THEN 
        JWR=JWR+8
	GO TO 800
	ELSE IF(VNAME2.EQ.'gluc') THEN 
        JWR=JWR+9
	GO TO 800
	ELSE IF(VNAME2.EQ.'prot') THEN 
        JWR=13
	GO TO 800
	ELSE IF(VNAME2.EQ.'osm') THEN 
        IWR=1+2*IC
        JWR=13
	GO TO 800
        ELSE IF(VNAME2.EQ.'ta') THEN
        JWR=11
          IF(VNAME1.EQ.'cc') THEN
          IWR=2+5*IC
	  GO TO 800
          ELSE IF(VNAME1.EQ.'fkc') THEN
          IWR=2+6*IC
	  GO TO 800
          ELSE IF(VNAME1.EQ.'apr') THEN
          JWR=13
	  GO TO 800
          ELSE
          PRINT 235
  235     FORMAT(' CANT FIND VNAME1 MATCHING TA, TRY AGAIN')
            IF(SWSOU.EQ.'int') THEN
            GO TO 200
            ELSE
            STOP
            ENDIF
          ENDIF
        ELSE IF(VNAME2.EQ.'nae') THEN
        JWR=12
          IF(VNAME1.EQ.'cc') THEN
          IWR=2+5*IC
	  GO TO 800
          ELSE IF(VNAME1.EQ.'fkc') THEN
          IWR=2+6*IC
	  GO TO 800
          ELSE IF(VNAME1.EQ.'apr') THEN
          JWR=14
	  GO TO 800
          ELSE
          PRINT 236
  236     FORMAT(' CANT FIND VNAME1 MATCHING NEA, TRY AGAIN')
            IF(SWSOU.EQ.'int') THEN
            GO TO 200
            ELSE
            STOP
            ENDIF
          ENDIF
        ELSE
        PRINT 240
  240   FORMAT(' NO SUCH SOLUTE NAME, TRY AGAIN')
        IF(SWSOU.EQ.'int') THEN
        GO TO 200
        ELSE
        STOP
        ENDIF
        ENDIF
C
C
C A HGB SPECIES HAS BEEN DESIGNATED
C
  250   IF(SWSOU.EQ.'int') THEN
        PRINT 201
        READ 30, VNAME2
        ELSE
        ENDIF
C
	IF(VNAME2.EQ.'onh2') THEN 
        IWR=2+5*IC
        JWR=4
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh2') THEN 
        IWR=2+5*IC
        JWR=5
	GO TO 800
	ELSE IF(VNAME2.EQ.'onh3') THEN 
        IWR=2+5*IC
        JWR=6
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh3') THEN 
        IWR=2+5*IC
        JWR=7
	GO TO 800
	ELSE IF(VNAME2.EQ.'onhco2') THEN 
        IWR=2+5*IC
        JWR=8
	GO TO 800
	ELSE IF(VNAME2.EQ.'nhco2') THEN 
        IWR=2+5*IC
        JWR=9
	GO TO 800
	ELSE IF(VNAME2.EQ.'oxh') THEN 
        IWR=2+5*IC
        JWR=2
	GO TO 800
	ELSE IF(VNAME2.EQ.'xh') THEN 
        IWR=2+6*IC
        JWR=3
	GO TO 800
	ELSE IF(VNAME2.EQ.'ox') THEN 
        IWR=2+6*IC
        JWR=4
	GO TO 800
	ELSE IF(VNAME2.EQ.'x') THEN 
        IWR=2+6*IC
        JWR=5
	GO TO 800
	ELSE IF(VNAME2.EQ.'oyh') THEN 
        IWR=2+6*IC
        JWR=6
	GO TO 800
	ELSE IF(VNAME2.EQ.'yh') THEN 
        IWR=2+6*IC
        JWR=7
	GO TO 800
	ELSE IF(VNAME2.EQ.'oy') THEN 
        IWR=2+6*IC
        JWR=8
	GO TO 800
	ELSE IF(VNAME2.EQ.'y') THEN 
        IWR=2+6*IC
        JWR=9
	GO TO 800
        ELSE
        PRINT 260
  260   FORMAT(' NO SUCH HGB NAME, TRY AGAIN')
        IF(SWSOU.EQ.'int') THEN
        GO TO 250
        ELSE
        STOP
        ENDIF
        ENDIF
C
  800   RETURN
        END