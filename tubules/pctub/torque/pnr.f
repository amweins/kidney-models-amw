	PROGRAM PNR
C  PROGRAM TO READ AN OUTPUT FILE FROM PCT.F
C   AND COMPUTE THE REABSORPTIVE FLUXES IN EACH SEGMENT
C
C
	INTEGER NWR,IWR(20),JWR(20),CHOP(3),DCHOP(3)
	REAL WRRES(6,20,1500),RDRES(20),WRDL(6,20,1500),DL(20),
     1    WREX(6,20,1500),EX(20),PRSC
        CHARACTER*12 IPTFIL,OPTFIL,VARNAM
        CHARACTER*5 VNAME1(20),VNAME2(20),SWI
	CHARACTER*4 SEGNAM(5)
C
	COMMON/TWR/NWR,IWR,JWR,RDRES,DL,EX
	DATA SEGNAM/' PCT','    ','    ','    ','    '/
C
	PRINT 2
    2   FORMAT(' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL
    4   FORMAT(A12)
	OPEN(1,FILE=IPTFIL)
	OPEN(2,FILE='tab.pnr')
	OPEN(3,FILE='xy.pnr')
C
	PRINT 3
    3   FORMAT(' MESH SPACING - DATA FILE OR INTERACTIVE? (d/i):',$)
	READ *, SWI
	IF (SWI.EQ.'i') GO TO 7
	OPEN(4,FILE='pctmesh.dat')
	READ (4,5) (CHOP(L),DCHOP(L),L=1,1)
    5   FORMAT(2I5)
	GO TO 12
C
    7   DO 10 K=1,1
	PRINT 8, SEGNAM(K)
    8   FORMAT(A4,' SPATIAL CHOP, AND PRINT SPACING:',$) 
	READ *, CHOP(K),DCHOP(K)
   10   CONTINUE
C
   12   DO 13 L=1,1
   13   CHOP(L)=CHOP(L)+1
C
	PRINT 500
  500   FORMAT(' IS INPUT FROM A DATA FILE OR INTERACTIVE? (d/i): ',$)
	READ 30, SWI
  510   IF (SWI.EQ.'i') GO TO 15
	PRINT 516
  516   FORMAT(' TYPE THE VARIABLE FILENAME:',$)
	READ 514, VARNAM
  514   FORMAT(A12)
	OPEN(5,FILE=VARNAM)
	READ (5,515) NWR
  515   FORMAT(3X,I2)
	DO 525 KWR=1,NWR
  525   READ (5,530) VNAME1(KWR),VNAME2(KWR)
  530   FORMAT(2A5)
	GO TO 801
C
C
   15   PRINT 16
   16   FORMAT(' TYPE THE NUMBER OF VARIABLES TO BE REWRITTEN:',$)
	READ *, NWR
	DO 50 KWR=1,NWR
	PRINT 17, KWR
   17   FORMAT(' TYPE THE NAME OF VARIABLE NUMBER',I2,":",2X,$)
	READ 30, VNAME1(KWR)
	IF(VNAME1(KWR).EQ.'cm') GO TO 19
	IF(VNAME1(KWR).EQ.'xm') GO TO 19
	IF(VNAME1(KWR).EQ.'fkm') GO TO 19
	GO TO 50	
   19   PRINT 20
   20   FORMAT(' TYPE THE SPECIES NUMBER:',2X,$)
   30   FORMAT(A5)
	READ 30, VNAME2(KWR)
   50   CONTINUE
C
C THE WRITTEN VARIABLES HAVE BEEN SPECIFIED, AS HAS TUBULE SEGMENTATION
C
  801   DO 900 KXP=1,1500
C
C DATA FOR EACH SEGMENT IS READ
C
	DO 820 ISEG=1,1
	DO 810 KWR=1,NWR
  810   CALL SELEC 
     1   (VNAME1(KWR),VNAME2(KWR),IWR(KWR),JWR(KWR),CHOP(ISEG))
	CALL RDSEG(CHOP(ISEG),*811,*901)
  811   DO 817 J=1,NWR
	PRSC=1.0
	IF(VNAME1(J).EQ.'fvm') PRSC=60.
	IF(VNAME1(J).EQ.'fkm') PRSC=60.
	WRDL(ISEG,J,KXP)=PRSC*DL(J)
	WREX(ISEG,J,KXP)=PRSC*EX(J)
  817   WRRES(ISEG,J,KXP)=PRSC*RDRES(J)
  820   CONTINUE
C	
  900   CONTINUE
C
C WRRES HAS BEEN FILLED. 
C
  901   KXP=KXP-1
	WRITE (2,902) (SEGNAM(I),I=1,1)
  902   FORMAT(///,10X,A4,8X)
	DO 910 K=1,KXP
	WRITE (2,903) 
  903   FORMAT(/,2X,'absolute delivery',/)
	DO 923 J=1,NWR
  923   WRITE (2,915) VNAME2(J),
     1    (WRDL(ISEG,J,K),ISEG=1,1),WREX(1,J,K)
	WRITE (2,904) 
  904   FORMAT(/,2X,'absolute reabsorption',/)
	DO 924 J=1,NWR
  924   WRITE (2,915) VNAME2(J),
     1    (WRRES(ISEG,J,K),ISEG=1,1)
	WRITE (2,905) 
  905   FORMAT(/,2X,'reabsorption relative to segmental delivery',/)
	DO 925 J=1,NWR
  925   WRITE (2,916) VNAME2(J),
     1     (WRRES(ISEG,J,K)/WRDL(ISEG,J,K),ISEG=1,1)
C
	DO 930 J=1,NWR
  930   WRITE (3,917) WRDL(1,J,K),
     1     (WRRES(ISEG,J,K),
     1     WRRES(ISEG,J,K)/WRDL(ISEG,J,K),ISEG=1,1)
C
  910   CONTINUE
  915   FORMAT (1X,A5,8E12.4)
  916   FORMAT (1X,A5,8F12.4)
  917   FORMAT (8E16.8)
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
	IF(VNAME2.EQ.'hco2') GO TO 632
	IF(VNAME2.EQ.'h2co2') GO TO 633
	IF(VNAME2.EQ.'gluc') GO TO 634
	IF(VNAME2.EQ.'h') GO TO 635
	IF(VNAME2.EQ.'impm') GO TO 635
	IF(VNAME2.EQ.'ta') GO TO 635
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
  633   JWR=JWR+12
	GO TO 700
  634   JWR=JWR+13
	GO TO 700
  635   JWR=JWR+14
	GO TO 700
C
  700   CONTINUE
	RETURN
	END
C
	SUBROUTINE RDSEG(CHOP,*,*)
C PROGRAM TO READ THE TUBULE SEGMENT RECORDS AND OUTPUT THE SELECTED VARIABLES
C
	INTEGER NWR,IWR(20),JWR(20),CHOP
	REAL RES(4803,18),TIME,RDRES(20),DL(20),EX(20)
C
	COMMON/TWR/NWR,IWR,JWR,RDRES,DL,EX
C
C
  300   READ (1,310,END=360) TIME
  310   FORMAT (2(/),11X,F15.4)
	READ (1,325) ((RES(I,J),J=1,18),I=1,CHOP)
  325   FORMAT (///,(F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,
     1   D9.3,2F9.6,D9.3,2F9.6))
	READ (1,335) ((RES(I,J),J=1,18),I=1+CHOP,2*CHOP)
  335   FORMAT (///,(F7.4,17F9.4))
	READ (1,345) ((RES(I,J),J=1,18),I=1+2*CHOP,3*CHOP)
  345   FORMAT (///,(F7.4,6F9.6,D9.3,4F9.6,
     1   D9.3,2F9.6,D9.3,2F9.6))
C
	DO 350 J=1,NWR
	DL(J) = RES(IWR(J),JWR(J))
	EX(J) = RES(IWR(J)+CHOP-1,JWR(J))
  350   RDRES(J)=RES(IWR(J),JWR(J)) - RES(IWR(J)+CHOP-1,JWR(J))
C
	RETURN 1
  360   RETURN 2
	END
