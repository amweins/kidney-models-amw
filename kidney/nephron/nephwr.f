	PROGRAM NEPHWR
C  PROGRAM TO READ LUMEN.DAT FROM NEPHRON
C
C  Calculation plan: 
C  SFPCT -> SFPST -> SDHL -> AHLm -> AHLc -> DCT
C  JMPCT -> JMPST -> LDHLu -> LDHLl -> tAHL -> AHLm -> AHLc -> DCT  (2-6)
C  Merged CNT -> CCD -> OMCD -> IMCD
C 
C  In all, there are 6 SF segments, 40 JM segments, and 4 ASDN segments,
C   so 50 segments and 150 pages of data for a single run.
C  In will be read as named data in anticipation for processing.
C
C Segment numbers (ISEG)
C      1- PCT
C      2- PST
C      3- sDHL, lDHLu
C      4- lDHLl
C      5- tAHL
C      6- AHLm
C      7- AHLc
C      8- DCT
C      9- CNT
C      10-CCD
C      11-OMCD
C      12-IMCD
C
C Plan for 24000 SF nephrons (@30 nl/min -> 720 mul/min) 
C    12400 JM nephrons (@60 nl/min -> 720 mul/min) 
C     Of the 12000 JM, turns can parallel CD coning:
C    1 mm 6400
C    2 mm 3200
C    3 mm 1600
C    4 mm  800
C    5 mm  400
C
C Index to solutes in FKM:
C	FKM(1)='  NA '
C	FKM(2)='  K  '
C	FKM(3)='  CL '
C	FKM(4)=' HCO3'
C	FKM(5)='H2CO3'
C	FKM(6)=' CO2 '
C	FKM(7)=' HPO4'
C	FKM(8)='H2PO4'
C	FKM(9)=' UREA'
C	FKM(10)=' NH3 '
C	FKM(11)=' NH4 '
C	FKM(12)='  H  '
C	FKM(13)=' HCO2'
C	FKM(14)='H2CO2'
C	FKM(15)=' GLUC'
C	FKM(16)=' IMPM'
C	FKM(17)='  TA '
C	FKM(18)=' NAE '
C
C For this program, the file lumen.dat willl be read into archive
C and then output according to (nephron,segment,distance)
C
C LUMINAL AND PERITUBULAR VARIABLES
        INTEGER SOLS,CHOP,DCHOP,JWR(20),
     1    RCHOP(6,12),RDCHOP(6,12),RMUM(6)
	DOUBLE PRECISION RDIST(6,12,901),
     1   RVM(6,12,901), RPM(6,12,901), RCM(6,12,15,901),
     1   RIMPM(6,12,901), RLCHM(6,12,901), RXM(6,12,15,901),
     1   ROSMM(7,12,901),RTA(6,12,901),RFVM(7,12,901), 
     1   RFKM(7,12,18,901),RNAE(6,12,901),
     1   RDPAGE(60,901),WRPAGE(10,901),SCAL(20)
	DOUBLE PRECISION DIST(901),
     1   VM(901),PM(901),CM(15,901),IMPM(901),
     1   LCHM(901),XM(15,901),OSMM(901),TA(901),
     1   FVM(901),FKM(18,901),NAE(901)
	DOUBLE PRECISION 
     1   WRRES(7,12,10),WRDL(7,12,10),WREX(7,12,10)
C
        CHARACTER*12 IPTFIL,OPTFIL,VARNAM
        CHARACTER*5 VNAME1(20),VNAME2(20),SWI,SWO,
     1    SEGSF(13),SEGJM(13),NEPNAM(7),BARNAM(10)
C
	COMMON/TWR/CHOP,DIST,VM,PM,CM,IMPM,
     1    LCHM,XM,OSMM,TA,FVM,FKM,NAE
	DATA 
     1    SEGSF /'  pct','  pst',' sdhl','     ',
     1    '     ',' ahlm',' ahlc','  dct',
     1    '  cnt','  ccd',' omcd',' imcd','total'/,
     1    SEGJM /'  pct','  pst','ldhlu','ldhll',
     1    ' tahl',' ahlm',' ahlc','  dct',
     1    '  cnt','  ccd',' omcd',' imcd','total'/,
     1    NEPNAM  /'   sf',
     1    '  jm1','  jm2','  jm3','  jm4','  jm5','alljm'/,
     1    RMUM /24000,6400,3200,1600,800,400/
C
        NNEPH=6
        NSEG=12
        IPTFIL = 'lumen.dat'
	OPEN(1,FILE=IPTFIL)
        OPEN(20,FILE='xy.na')
        OPEN(21,FILE='xy.k')
        OPEN(22,FILE='xy.hco3')
        OPEN(23,FILE='xy.nh4')
        OPEN(24,FILE='xy.nae')
        OPEN(25,FILE='xy.urea')
        OPEN(26,FILE='xy.h2o')
        OPEN(30,FILE='bars.na')
        OPEN(31,FILE='bars.k')
        OPEN(32,FILE='bars.hco3')
        OPEN(33,FILE='bars.nh4')
        OPEN(34,FILE='bars.nae')
        OPEN(35,FILE='bars.urea')
        OPEN(36,FILE='bars.h2o')
        OPEN(40,FILE='tab.dat')
C
	OPEN(3,FILE='nephmesh.dat')
        DO 5 ISEG=1,NSEG
    5   READ (3,8) (RCHOP(L,ISEG),RDCHOP(L,ISEG),L=1,NNEPH)
    8   FORMAT(12I5)
C
C
C DATA FOR EACH SEGMENT IS READ: 
C  FIRST SF (SFPCT -> DCT), THEN JM2-6 (JMPCT -> DCT), AND FINALLY CNT -> IMCT
C
        INEPH=1
        DO 100 ISEG=1,3
        CHOP=RCHOP(INEPH,ISEG)
        IF (ISEG.LE.2) THEN
        SOLS=15
        CALL RDSEG(1,*92,*190)
        ELSE
        SOLS=12
        CALL RDSEG(2,*92,*190)
        ENDIF
   92   CONTINUE
C
C LUMINAL VARIABLES
	DO 96 IX = 1,CHOP+1
	RDIST(INEPH,ISEG,IX) = DIST(IX)
	RVM(INEPH,ISEG,IX) = VM(IX)
	RPM(INEPH,ISEG,IX) = PM(IX)
	RIMPM(INEPH,ISEG,IX) = IMPM(IX)
	RLCHM(INEPH,ISEG,IX) = LCHM(IX)
	ROSMM(INEPH,ISEG,IX) = OSMM(IX)
	RTA(INEPH,ISEG,IX) = TA(IX)/1.D3
	RFVM(INEPH,ISEG,IX) = FVM(IX)
	DO 95 ISOL = 1,SOLS
	RCM(INEPH,ISEG,ISOL,IX) = CM(ISOL,IX)
	RXM(INEPH,ISEG,ISOL,IX) = XM(ISOL,IX)
	RFKM(INEPH,ISEG,ISOL,IX) = FKM(ISOL,IX)
   95   CONTINUE
	RFKM(INEPH,ISEG,16,IX) = FKM(16,IX)
	RFKM(INEPH,ISEG,17,IX) = FKM(17,IX)
        NAE(IX) = CM(11,IX)-CM(4,IX)+TA(IX)/1.D3 
        FKM(18,IX) = FKM(11,IX)-FKM(4,IX)+FKM(17,IX)
        RNAE(INEPH,ISEG,IX)=NAE(IX)
        RFKM(INEPH,ISEG,18,IX)=FKM(18,IX)
   96   CONTINUE
  100   CONTINUE
C
        INEPH=1
        DO 110 ISEG=6,8
        CHOP=RCHOP(INEPH,ISEG)
        CALL RDSEG(2,*102,*190)
  102   SOLS=12
C
C LUMINAL VARIABLES
	DO 106 IX = 1,CHOP+1
	RDIST(INEPH,ISEG,IX) = DIST(IX) + 1.0
	RVM(INEPH,ISEG,IX) = VM(IX)
	RPM(INEPH,ISEG,IX) = PM(IX)
	RIMPM(INEPH,ISEG,IX) = IMPM(IX)
	RLCHM(INEPH,ISEG,IX) = LCHM(IX)
	ROSMM(INEPH,ISEG,IX) = OSMM(IX)
	RTA(INEPH,ISEG,IX) = TA(IX)/1.D3
	RFVM(INEPH,ISEG,IX) = FVM(IX)
	DO 105 ISOL = 1,SOLS
	RCM(INEPH,ISEG,ISOL,IX) = CM(ISOL,IX)
	RXM(INEPH,ISEG,ISOL,IX) = XM(ISOL,IX)
	RFKM(INEPH,ISEG,ISOL,IX) = FKM(ISOL,IX)
  105   CONTINUE
	RFKM(INEPH,ISEG,16,IX) = FKM(16,IX)
	RFKM(INEPH,ISEG,17,IX) = FKM(17,IX)
        NAE(IX) = CM(11,IX)-CM(4,IX)+TA(IX)/1.D3 
        FKM(18,IX) = FKM(11,IX)-FKM(4,IX)+FKM(17,IX)
        RNAE(INEPH,ISEG,IX)=NAE(IX)
        RFKM(INEPH,ISEG,18,IX)=FKM(18,IX)
  106   CONTINUE
  110   CONTINUE
C
        DO 120 INEPH=2,6
        DO 120 ISEG=1,8
        CHOP=RCHOP(INEPH,ISEG)
        IF (ISEG.LE.2) THEN
        SOLS=15
        CALL RDSEG(1,*112,*190)
        ELSE
        SOLS=12
        CALL RDSEG(2,*112,*190)
        ENDIF
  112   CONTINUE
C
C LUMINAL VARIABLES
	DO 116 IX = 1,CHOP+1
        IF (ISEG.LE.4) THEN
	RDIST(INEPH,ISEG,IX) = DIST(IX)
        ELSE
	RDIST(INEPH,ISEG,IX) = DIST(IX) + 0.2*FLOAT(6-INEPH)
        ENDIF
	RVM(INEPH,ISEG,IX) = VM(IX)
	RPM(INEPH,ISEG,IX) = PM(IX)
	RIMPM(INEPH,ISEG,IX) = IMPM(IX)
	RLCHM(INEPH,ISEG,IX) = LCHM(IX)
	ROSMM(INEPH,ISEG,IX) = OSMM(IX)
	RTA(INEPH,ISEG,IX) = TA(IX)/1.D3
	RFVM(INEPH,ISEG,IX) = FVM(IX)
	DO 115 ISOL = 1,SOLS
	RCM(INEPH,ISEG,ISOL,IX) = CM(ISOL,IX)
	RXM(INEPH,ISEG,ISOL,IX) = XM(ISOL,IX)
	RFKM(INEPH,ISEG,ISOL,IX) = FKM(ISOL,IX)
  115   CONTINUE
	RFKM(INEPH,ISEG,16,IX) = FKM(16,IX)
	RFKM(INEPH,ISEG,17,IX) = FKM(17,IX)
        NAE(IX) = CM(11,IX)-CM(4,IX)+TA(IX)/1.D3 
        FKM(18,IX) = FKM(11,IX)-FKM(4,IX)+FKM(17,IX)
        RNAE(INEPH,ISEG,IX)=NAE(IX)
        RFKM(INEPH,ISEG,18,IX)=FKM(18,IX)
  116   CONTINUE
  120   CONTINUE
C
        INEPH=1
        DO 130 ISEG=9,12
        CHOP=RCHOP(INEPH,ISEG)
        CALL RDSEG(2,*122,*190)
  122   SOLS=12
C
C LUMINAL VARIABLES
	DO 126 IX = 1,CHOP+1
	RDIST(INEPH,ISEG,IX) = DIST(IX)
	RVM(INEPH,ISEG,IX) = VM(IX)
	RPM(INEPH,ISEG,IX) = PM(IX)
	RIMPM(INEPH,ISEG,IX) = IMPM(IX)
	RLCHM(INEPH,ISEG,IX) = LCHM(IX)
	ROSMM(INEPH,ISEG,IX) = OSMM(IX)
	RTA(INEPH,ISEG,IX) = TA(IX)
	RFVM(INEPH,ISEG,IX) = FVM(IX)
	DO 125 ISOL = 1,SOLS
	RCM(INEPH,ISEG,ISOL,IX) = CM(ISOL,IX)
	RXM(INEPH,ISEG,ISOL,IX) = XM(ISOL,IX)
	RFKM(INEPH,ISEG,ISOL,IX) = FKM(ISOL,IX)
  125   CONTINUE
	RFKM(INEPH,ISEG,16,IX) = FKM(16,IX)
	RFKM(INEPH,ISEG,17,IX) = FKM(17,IX)
        NAE(IX) = CM(11,IX)-CM(4,IX)+TA(IX)/1.D3 
        FKM(18,IX) = FKM(11,IX)-FKM(4,IX)+FKM(17,IX)
        RNAE(INEPH,ISEG,IX)=NAE(IX)
        RFKM(INEPH,ISEG,18,IX)=FKM(18,IX)
  126   CONTINUE
  130   CONTINUE
C
C ALLOW FOR THE POSSIBILITY FOR OUTPUT OF TOTAL JM FLOWS (INEPH=7)
C  FOR THE LONG THIN LIMBS (LDHLL AND TAHL) THERE NEEDS TO BE CARE FOR THE NUMBERING
C
        DO 150 ISEG=1,3
        CHOP=RCHOP(6,ISEG)
        DO 150 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(6,ISEG,IX)*RFVM(6,ISEG,IX)
        RFVM(7,ISEG,IX)=RFVM(6,ISEG,IX)
        DO 148 INEPH=2,5
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)+
     1    ROSMM(INEPH,ISEG,IX)*RFVM(INEPH,ISEG,IX)
  148   RFVM(7,ISEG,IX)=RFVM(7,ISEG,IX)+RFVM(INEPH,ISEG,IX)
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)/RFVM(7,ISEG,IX)
        DO 149 ISOL=1,18
        RFKM(7,ISEG,ISOL,IX)=RFKM(6,ISEG,ISOL,IX)
        DO 149 INEPH=2,5
  149   RFKM(7,ISEG,ISOL,IX)=
     1   RFKM(7,ISEG,ISOL,IX)+RFKM(INEPH,ISEG,ISOL,IX)
  150   CONTINUE
C
        DO 160 ISEG=6,8
        CHOP=RCHOP(6,ISEG)
        DO 160 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(6,ISEG,IX)*RFVM(6,ISEG,IX)
        RFVM(7,ISEG,IX)=RFVM(6,ISEG,IX)
        DO 158 INEPH=2,5
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)+
     1    ROSMM(INEPH,ISEG,IX)*RFVM(INEPH,ISEG,IX)
  158   RFVM(7,ISEG,IX)=RFVM(7,ISEG,IX)+RFVM(INEPH,ISEG,IX)
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)/RFVM(7,ISEG,IX)
        DO 159 ISOL=1,18
        RFKM(7,ISEG,ISOL,IX)=RFKM(6,ISEG,ISOL,IX)
        DO 159 INEPH=2,5
  159   RFKM(7,ISEG,ISOL,IX)=
     1   RFKM(7,ISEG,ISOL,IX)+RFKM(INEPH,ISEG,ISOL,IX)
  160   CONTINUE
C
        ISEG=4
        CHOP=RCHOP(6,ISEG)
        DO 165 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(6,ISEG,IX)*RFVM(6,ISEG,IX)
        RFVM(7,ISEG,IX)=RFVM(6,ISEG,IX)
        DO 165 ISOL=1,18
        RFKM(7,ISEG,ISOL,IX)=RFKM(6,ISEG,ISOL,IX)
  165   CONTINUE
        DO 170 INEPH=5,2,-1
        CHOP=RCHOP(INEPH,ISEG)
        DO 170 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)+
     1    ROSMM(INEPH,ISEG,IX)*RFVM(INEPH,ISEG,IX)
        RFVM(7,ISEG,IX)=RFVM(7,ISEG,IX)+RFVM(INEPH,ISEG,IX)
        DO 169 ISOL=1,18
        RFKM(7,ISEG,ISOL,IX)=
     1   RFKM(7,ISEG,ISOL,IX)+RFKM(INEPH,ISEG,ISOL,IX)
  169   CONTINUE
  170   CONTINUE
        CHOP=RCHOP(6,ISEG)
        DO 171 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)/RFVM(7,ISEG,IX)
  171   CONTINUE
C
        ISEG=5
        CHOP=RCHOP(6,ISEG)
        ICHOP6=CHOP
        DO 175 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(6,ISEG,IX)*RFVM(6,ISEG,IX)
        RFVM(7,ISEG,IX)=RFVM(6,ISEG,IX)
        DO 175 ISOL=1,18
        RFKM(7,ISEG,ISOL,IX)=RFKM(6,ISEG,ISOL,IX)
  175   CONTINUE
        DO 180 INEPH=5,2,-1
        CHOP=RCHOP(INEPH,ISEG)
        DO 180 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)+
     1    ROSMM(INEPH,ISEG,IX)*RFVM(INEPH,ISEG,IX)
        RFVM(7,ISEG,ICHOP6-CHOP+IX)=
     1    RFVM(7,ISEG,ICHOP6-CHOP+IX)+RFVM(INEPH,ISEG,IX)
        DO 179 ISOL=1,18
        RFKM(7,ISEG,ISOL,ICHOP6-CHOP+IX)=
     1   RFKM(7,ISEG,ISOL,ICHOP6-CHOP+IX)+RFKM(INEPH,ISEG,ISOL,IX)
  179   CONTINUE
  180   CONTINUE
        CHOP=RCHOP(6,ISEG)
        DO 181 IX=1,CHOP+1
        ROSMM(7,ISEG,IX)=ROSMM(7,ISEG,IX)/RFVM(7,ISEG,IX)
  181   CONTINUE
C
        GO TO 199
C
  190   PRINT 192, INEPH,ISEG
  192   FORMAT ('PREMATURE EOF:',2X,'INEPH= ', I2, 10X,'ISEG= ',I2)
        STOP
C
C ALL VARIABLES HAVE BEEN READ INTO THE ARCHIVE
C THE SELECTED VARIABLES SHOULD BE PRINTED INTO A FILE READABLE BY GRX:
C   SEGMENT BY SEGMENT EACH VARIABLE IS PRINTED
C FIRST A SURROGATE LUMEN.DAT PAGE IS CREATED:
C        DIST   VM    PM     CM    IMPM   TA   LCHM    XM     OSMM    FVM     FKM     
C COL:    1     2     3     4-18    19    20    21    22-36    37     38     39-56
C
C One optional output file may be generated, along with the standard output files.
C
  199   DO 300 KK=1,7
        IUNIT=19+KK
        JUNIT=29+KK
        NWR=2
        SCAL(1)=1.D1
        SCAL(2)=1.D3
        SCAL(3)=1.D3
        IF (IUNIT.EQ.20) THEN
        VNAME1(2) = 'cm'
        VNAME2(2) = 'na'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'na'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.21) THEN
        VNAME1(2) = 'cm'
        VNAME2(2) = 'k'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'k'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.22) THEN
        VNAME1(2) = 'cm'
        VNAME2(2) = 'hco3'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'hco3'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.23) THEN
        VNAME1(2) = 'cm'
        VNAME2(2) = 'nh4'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'nh4'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.24) THEN
        VNAME1(2) = 'nae'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'nae'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.25) THEN
        VNAME1(2) = 'cm'
        VNAME2(2) = 'urea'
        VNAME1(3) = 'fkm'
        VNAME2(3) = 'urea'
        BARNAM(KK)=VNAME2(3)
        ELSE IF (IUNIT.EQ.26) THEN
        VNAME1(2) = 'osmm'
        VNAME1(3) = 'fvm'
        BARNAM(KK)=VNAME1(3)
        ENDIF
C
        JWR(1)=1
        DO 197 KWR=2,NWR+1
  197   CALL SELEC (VNAME1(KWR),VNAME2(KWR),JWR(KWR))
C
  200   CONTINUE
        INEPH=1
        DO 220 ISEG=1,3
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 204 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 202 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  202   CONTINUE
  204   CONTINUE
        DO 208 IX=1,CHOP+1
        DO 206 KWR=1,NWR+1
  206   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  208   CONTINUE
C
        WRDL(INEPH,ISEG,KK)=WRPAGE(3,1)
        WREX(INEPH,ISEG,KK)=WRPAGE(3,CHOP+1)
        WRRES(INEPH,ISEG,KK)=WRDL(INEPH,ISEG,KK)-WREX(INEPH,ISEG,KK)
C
        WRITE(IUNIT,298) CHOP/DCHOP
  298   FORMAT(I3)
        DO 210 IX=1,CHOP+1,DCHOP
  210   WRITE(IUNIT,299) (SCAL(KWR)*WRPAGE(KWR,IX),KWR=1,NWR+1)
  299   FORMAT (8D16.8)
  220   CONTINUE
C
        DO 240 ISEG=6,8
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 224 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 222 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  222   CONTINUE
  224   CONTINUE
        DO 228 IX=1,CHOP+1
        DO 226 KWR=1,NWR+1
  226   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  228   CONTINUE
C
        WRDL(INEPH,ISEG,KK)=WRPAGE(3,1)
        WREX(INEPH,ISEG,KK)=WRPAGE(3,CHOP+1)
        WRRES(INEPH,ISEG,KK)=WRDL(INEPH,ISEG,KK)-WREX(INEPH,ISEG,KK)
C
        WRITE(IUNIT,298) CHOP/DCHOP
        DO 230 IX=1,CHOP+1,DCHOP
  230   WRITE(IUNIT,299) (SCAL(KWR)*WRPAGE(KWR,IX),KWR=1,NWR+1)
  240   CONTINUE
C
C What is below if full output of all juxtamedullary segments.  
C We will try to graph an average JM concentration and total JM flows
C
        DO 241 ISEG=1,8
        WRDL(7,ISEG,KK)=0.D0
        WREX(7,ISEG,KK)=0.D0
        WRRES(7,ISEG,KK)=0.D0
  241   CONTINUE
C
        DO 260 INEPH=2,6
        DO 260 ISEG=1,8
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 244 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 242 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  242   CONTINUE
  244   CONTINUE
        DO 248 IX=1,CHOP+1
        DO 246 KWR=1,NWR+1
  246   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  248   CONTINUE
C
        WRDL(INEPH,ISEG,KK)=WRPAGE(3,1)
        WREX(INEPH,ISEG,KK)=WRPAGE(3,CHOP+1)
        WRRES(INEPH,ISEG,KK)=WRDL(INEPH,ISEG,KK)-WREX(INEPH,ISEG,KK)
C
        WRDL(7,ISEG,KK)=WRDL(7,ISEG,KK) + WRDL(INEPH,ISEG,KK)
        WREX(7,ISEG,KK)=WREX(7,ISEG,KK) + WREX(INEPH,ISEG,KK)
        WRRES(7,ISEG,KK)=WRRES(7,ISEG,KK) + WRRES(INEPH,ISEG,KK)
C
        GO TO 255 
        WRITE(IUNIT,298) CHOP/DCHOP
        DO 250 IX=1,CHOP+1,DCHOP
  250   WRITE(IUNIT,299) (SCAL(KWR)*WRPAGE(KWR,IX),KWR=1,NWR+1)
  255   CONTINUE
  260   CONTINUE
C
C
  275   CONTINUE
        INEPH=6
        DO 290 ISEG=1,8
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 284 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(7,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(7,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(7,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(7,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 282 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(7,ISEG,ISOL,IX)
  282   CONTINUE
  284   CONTINUE
        DO 288 IX=1,CHOP+1
        DO 286 KWR=1,NWR+1
  286   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  288   CONTINUE
C
        WRITE(IUNIT,298) CHOP/DCHOP
        DO 289 IX=1,CHOP+1,DCHOP
  289   WRITE(IUNIT,299) (SCAL(KWR)*WRPAGE(KWR,IX),KWR=1,NWR+1)
  290   CONTINUE
C
C
        INEPH=1
        DO 280 ISEG=9,12
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 264 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 262 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  262   CONTINUE
  264   CONTINUE
        DO 268 IX=1,CHOP+1
        DO 266 KWR=1,NWR+1
  266   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  268   CONTINUE
C
        WRDL(INEPH,ISEG,KK)=WRPAGE(3,1)
        WREX(INEPH,ISEG,KK)=WRPAGE(3,CHOP+1)
        WRRES(INEPH,ISEG,KK)=WRDL(INEPH,ISEG,KK)-WREX(INEPH,ISEG,KK)
C
        WRITE(IUNIT,298) CHOP/DCHOP
        DO 270 IX=1,CHOP+1,DCHOP
  270   WRITE(IUNIT,299) (SCAL(KWR)*WRPAGE(KWR,IX),KWR=1,NWR+1)
  280   CONTINUE
C
C Bargraph data for SF and composite JM and then distal nephron
C
        IF ((KK.EQ.2).OR.(KK.EQ.4).OR.(KK.EQ.6)) THEN
	NBARS=18
        ELSE 
        NBARS=4
        ENDIF
	WRITE (JUNIT,965) NBARS
	ISEG=0
	WRITE (JUNIT,965) ISEG, 1.D3*WRDL(1,1,KK),1
	WRITE (JUNIT,965) ISEG, 1.D3*WRDL(7,1,KK),2
	ISEG=1
        WRITE (JUNIT,965) ISEG, -1.D3*WRRES(1,ISEG,KK),1
        WRITE (JUNIT,965) ISEG, -1.D3*WRRES(7,ISEG,KK),2
C
        IF ((KK.NE.2).AND.(KK.NE.4).AND.(KK.NE.6)) THEN
	NBARS=14
	WRITE (JUNIT,965) NBARS
        ELSE 
        ENDIF
	DO 962 ISEG=2,3
        WRITE (JUNIT,965) ISEG, -1.D3*WRRES(1,ISEG,KK),1
  962   WRITE (JUNIT,965) ISEG, -1.D3*WRRES(7,ISEG,KK),2
        WRITE (JUNIT,965) 4, -1.D3*WRRES(7,4,KK),2
        WRITE (JUNIT,965) 5, -1.D3*WRRES(7,5,KK),2
	DO 964 ISEG=6,8
        WRITE (JUNIT,965) ISEG, -1.D3*WRRES(1,ISEG,KK),1
  964   WRITE (JUNIT,965) ISEG, -1.D3*WRRES(7,ISEG,KK),2
C
        WRITE (JUNIT,965) 9, 1.D3*WREX(1,8,KK),1
        WRITE (JUNIT,965) 9, 1.D3*WREX(7,8,KK),2
C
	NBARS=6
	WRITE (JUNIT,965) NBARS
	ISEG=0
	WRITE (JUNIT,965) ISEG, 1.D3*WRDL(1,9,KK),3
	DO 968 ISEG=9,12
  968   WRITE (JUNIT,965) ISEG, -1.D3*WRRES(1,ISEG,KK),3
        ISEG=13
	WRITE (JUNIT,965) ISEG, 1.D3*WREX(1,12,KK),3
  965   FORMAT (I2,F9.3,2X,I2)
C
  300   CONTINUE
C
C Output tabular data
C
        INEPH=1
	WRITE (40,902) NEPNAM(INEPH),
     1   (SEGSF(ISEG),ISEG=1,3),(SEGSF(ISEG),ISEG=6,8)
  902   FORMAT(///,2X,A5,3X,13(A5,7X))
	WRITE (40,903) 
  903   FORMAT(/,2X,'absolute delivery',/)
	DO 923 KK=1,7
  923   WRITE (40,915) BARNAM(KK),
     1    (WRDL(INEPH,ISEG,KK),ISEG=1,3),
     1    (WRDL(INEPH,ISEG,KK),ISEG=6,8)
	WRITE (40,904) 
  904   FORMAT(/,2X,'absolute reabsorption',/)
	DO 924 KK=1,7
  924   WRITE (40,915) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK),ISEG=1,3),
     1    (WRRES(INEPH,ISEG,KK),ISEG=6,8)
	WRITE (40,905) 
  905   FORMAT(/,2X,'reabsorption relative to segmental delivery',/)
	DO 925 KK=1,7
  925   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,ISEG,KK),ISEG=1,3),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,ISEG,KK),ISEG=6,8)
	WRITE (40,906) 
  906   FORMAT(/,2X,'reabsorption relative to pct delivery',/)
	DO 926 KK=1,7
  926   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,1,KK),ISEG=1,3),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,1,KK),ISEG=6,8)
  915   FORMAT (1X,A5,13E12.4)
  916   FORMAT (1X,A5,13F12.4)
C
        DO 940 INEPH=2,7
	WRITE (40,902) NEPNAM(INEPH),
     1   (SEGJM(ISEG),ISEG=1,8)
	WRITE (40,903) 
	DO 933 KK=1,7
  933   WRITE (40,915) BARNAM(KK),
     1    (WRDL(INEPH,ISEG,KK),ISEG=1,8)
	WRITE (40,904) 
	DO 934 KK=1,7
  934   WRITE (40,915) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK),ISEG=1,8)
	WRITE (40,905) 
	DO 935 KK=1,7
  935   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,ISEG,KK),ISEG=1,8)
	WRITE (40,906) 
	DO 936 KK=1,7
  936   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,1,KK),ISEG=1,8)
  940   CONTINUE
C
        INEPH=1
	WRITE (40,952) 
     1   (SEGJM(ISEG),ISEG=9,13)
  952   FORMAT(///,10X,13(A5,7X))
	WRITE (40,903) 
	DO 953 KK=1,7
  953   WRITE (40,915) BARNAM(KK),
     1    (WRDL(INEPH,ISEG,KK),ISEG=9,12),
     1    WREX(INEPH,12,KK)
	WRITE (40,904) 
	DO 954 KK=1,7
  954   WRITE (40,915) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK),ISEG=9,12),
     1    (WRDL(1,1,KK)+WRDL(7,1,KK))-WREX(INEPH,12,KK)
	WRITE (40,905) 
	DO 955 KK=1,7
  955   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,ISEG,KK),ISEG=9,12),
     1    1.D0-WREX(INEPH,12,KK)/(WRDL(1,1,KK)+WRDL(7,1,KK))
	WRITE (40,906) 
	DO 956 KK=1,7
  956   WRITE (40,916) BARNAM(KK),
     1    (WRRES(INEPH,ISEG,KK)/WRDL(INEPH,1,KK),ISEG=9,12),
     1    1.D0-WREX(INEPH,12,KK)/(WRDL(1,1,KK)+WRDL(7,1,KK))
C
C
        PRINT 2
    2   FORMAT (' DO YOU NEED A CUSTOM OUTPUT FILE? (y,n): ',$)
        READ 30, SWO
        IF (SWO.EQ.'n') GO TO 400
        IUNIT=2
    4   FORMAT(A12)
	PRINT 6
    6   FORMAT(' TYPE THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
	OPEN(IUNIT,FILE=OPTFIL)
C
	PRINT 10
   10   FORMAT(' IS INPUT FROM A DATA FILE OR INTERACTIVE? (d/i): ',$)
	READ 30, SWI
        IF (SWI.EQ.'i') GO TO 15
	PRINT 11
   11   FORMAT(' TYPE THE VARIABLE FILENAME:',$)
	READ 4, VARNAM
	OPEN(4,FILE=VARNAM)
	READ (4,12) NWR
   12   FORMAT(3X,I2)
	DO 13 KWR=2,NWR+1
   13   READ (4,14) VNAME1(KWR),VNAME2(KWR)
   14   FORMAT(2A5)
	GO TO 55
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
   55   JWR(1)=1
        DO 57 KWR=2,NWR+1
   57   CALL SELEC (VNAME1(KWR),VNAME2(KWR),JWR(KWR))
C   
        INEPH=1
        DO 320 ISEG=1,3
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 304 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 302 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  302   CONTINUE
  304   CONTINUE
        DO 308 IX=1,CHOP+1
        DO 306 KWR=1,NWR+1
  306   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  308   CONTINUE
C
        WRITE(IUNIT,398) CHOP/DCHOP
  398   FORMAT(I3)
        DO 310 IX=1,CHOP+1,DCHOP
  310   WRITE(IUNIT,399) (WRPAGE(KWR,IX),KWR=1,NWR+1)
  399   FORMAT (8D16.8)
  320   CONTINUE
C
        DO 340 ISEG=6,8
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 324 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 322 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  322   CONTINUE
  324   CONTINUE
        DO 328 IX=1,CHOP+1
        DO 326 KWR=1,NWR+1
  326   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  328   CONTINUE
C
        WRITE(IUNIT,398) CHOP/DCHOP
        DO 330 IX=1,CHOP+1,DCHOP
  330   WRITE(IUNIT,399) (WRPAGE(KWR,IX),KWR=1,NWR+1)
  340   CONTINUE
C
        DO 360 INEPH=3,6,3
        DO 360 ISEG=1,8
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 344 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 342 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  342   CONTINUE
  344   CONTINUE
        DO 348 IX=1,CHOP+1
        DO 346 KWR=1,NWR+1
  346   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  348   CONTINUE
C
        WRITE(IUNIT,398) CHOP/DCHOP
        DO 350 IX=1,CHOP+1,DCHOP
  350   WRITE(IUNIT,399) (WRPAGE(KWR,IX),KWR=1,NWR+1)
  360   CONTINUE
C
        INEPH=1
        DO 380 ISEG=9,12
        CHOP=RCHOP(INEPH,ISEG)
        DCHOP=RDCHOP(INEPH,ISEG)
        DO 364 IX=1,CHOP+1
        RDPAGE(1,IX)=RDIST(INEPH,ISEG,IX)
        RDPAGE(2,IX)=RVM(INEPH,ISEG,IX)
        RDPAGE(3,IX)=RPM(INEPH,ISEG,IX)
        RDPAGE(19,IX)=RIMPM(INEPH,ISEG,IX)
        RDPAGE(20,IX)=RTA(INEPH,ISEG,IX)
        RDPAGE(21,IX)=RLCHM(INEPH,ISEG,IX)
        RDPAGE(37,IX)=ROSMM(INEPH,ISEG,IX)
        RDPAGE(38,IX)=RFVM(INEPH,ISEG,IX)/SCAL(3)
        RDPAGE(54,IX)=RFKM(INEPH,ISEG,16,IX)
        RDPAGE(55,IX)=RFKM(INEPH,ISEG,17,IX)
        RDPAGE(56,IX)=RFKM(INEPH,ISEG,18,IX)
        RDPAGE(57,IX)=RNAE(INEPH,ISEG,IX)
        DO 362 ISOL=1,15
        RDPAGE(3+ISOL,IX)=RCM(INEPH,ISEG,ISOL,IX)
        RDPAGE(21+ISOL,IX)=RXM(INEPH,ISEG,ISOL,IX)
        RDPAGE(38+ISOL,IX)=RFKM(INEPH,ISEG,ISOL,IX)
  362   CONTINUE
  364   CONTINUE
        DO 368 IX=1,CHOP+1
        DO 366 KWR=1,NWR+1
  366   WRPAGE(KWR,IX)=RDPAGE(JWR(KWR),IX)
  368   CONTINUE
C
        WRITE(IUNIT,398) CHOP/DCHOP
        DO 370 IX=1,CHOP+1,DCHOP
  370   WRITE(IUNIT,399) (WRPAGE(KWR,IX),KWR=1,NWR+1)
  380   CONTINUE
  400   CONTINUE
C
	STOP
	END
C
C
	SUBROUTINE SELEC(VNAME1,VNAME2,JWR)
C PROGRAM TO ASSIGN THE ARRAY COORDINATES FOR THE OUTPUT VARIABLES
C
C FIRST A SURROGATE LUMEN.DAT PAGE IS CREATED:
C        DIST   VM    PM     CM    IMPM   TA   LCHM    XM     OSMM    FVM     FKM    NAE 
C COL:    1     2     3     4-18    19    20    21    22-36    37     38     39-56    57
C
	INTEGER JWR
        CHARACTER*5 VNAME1,VNAME2
C
        IF (VNAME1.EQ.'vm') THEN
        JWR=2
        RETURN
        ELSE IF (VNAME1.EQ.'pm') THEN
        JWR=3
        RETURN
        ELSE IF (VNAME1.EQ.'cm') THEN
        JWR=3
        GO TO 20
        ELSE IF (VNAME1.EQ.'impm') THEN
        JWR=19
        RETURN
        ELSE IF (VNAME1.EQ.'ta') THEN
        JWR=20
        RETURN
        ELSE IF (VNAME1.EQ.'lchm') THEN
        JWR=21
        RETURN
        ELSE IF (VNAME1.EQ.'nae') THEN
        JWR=57
        RETURN
        ELSE IF (VNAME1.EQ.'xm') THEN
        JWR=21
        GO TO 20
        ELSE IF (VNAME1.EQ.'osmm') THEN
        JWR=37
        RETURN
        ELSE IF (VNAME1.EQ.'fvm') THEN
        JWR=38
        RETURN
        ELSE IF (VNAME1.EQ.'fkm') THEN
        JWR=38
        GO TO 20
        ENDIF
C
   20   CONTINUE
	IF(VNAME2.EQ.'na') THEN
        JWR=JWR+1
	ELSE IF(VNAME2.EQ.'k') THEN
        JWR=JWR+2
	ELSE IF(VNAME2.EQ.'cl') THEN
        JWR=JWR+3
	ELSE IF(VNAME2.EQ.'hco3') THEN
        JWR=JWR+4
	ELSE IF(VNAME2.EQ.'h2co3') THEN
        JWR=JWR+5
	ELSE IF(VNAME2.EQ.'co2') THEN
        JWR=JWR+6
	ELSE IF(VNAME2.EQ.'hpo4') THEN
        JWR=JWR+7
	ELSE IF(VNAME2.EQ.'h2po4') THEN
        JWR=JWR+8
	ELSE IF(VNAME2.EQ.'urea') THEN
        JWR=JWR+9
	ELSE IF(VNAME2.EQ.'nh3') THEN
        JWR=JWR+10
	ELSE IF(VNAME2.EQ.'nh4') THEN
        JWR=JWR+11
	ELSE IF(VNAME2.EQ.'h') THEN
        JWR=JWR+12
	ELSE IF(VNAME2.EQ.'gluc') THEN
        JWR=JWR+15
	ELSE IF(VNAME2.EQ.'impm') THEN
        JWR=JWR+16
	ELSE IF(VNAME2.EQ.'ta') THEN
        JWR=JWR+17
	ELSE IF(VNAME2.EQ.'nae') THEN
        JWR=JWR+18
        ENDIF
C
	RETURN
	END
C
C
	SUBROUTINE RDSEG(IRD,*,*)
C PROGRAM TO READ THE TUBULE SEGMENT RECORDS AND OUTPUT THE SELECTED VARIABLES
C
        INTEGER CHOP,IRD
	DOUBLE PRECISION DIST(901),
     1   VM(901),PM(901),CM(15,901),IMPM(901),
     1   LCHM(901),XM(15,901),OSMM(901),TA(901),
     1   FVM(901),FKM(18,901),NAE(901)
C
	COMMON/TWR/CHOP,DIST,VM,PM,CM,IMPM,
     1    LCHM,XM,OSMM,TA,FVM,FKM,NAE
C
C
        GO TO (100,200) IRD
C
  100   CONTINUE
        READ (1,104)
        DO 110 KX=1,CHOP+1
  110   READ (1,105,END=300) DIST(KX),
     1  VM(KX),PM(KX),(CM(I,KX),I=1,11),
     1  (CM(I,KX),I=13,15),IMPM(KX)
  104   FORMAT (5/)
  105   FORMAT (F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,
     1   D9.3,2F9.6,D9.3,2F9.6)
C
        READ (1,114)
        DO 120 KX=1,CHOP+1
  120   READ (1,115) DIST(KX),TA(KX),LCHM(KX),
     1    (XM(I,KX),I=1,11),(XM(I,KX),I=13,15),XM(12,KX)
  114   FORMAT (2/)
  115   FORMAT (F7.4,17F9.4)
C
        READ (1,124)
        DO 130 KX=1,CHOP+1
  130   READ (1,125) DIST(KX),OSMM(KX),
     1    FVM(KX),(FKM(I,KX),I=1,11),(FKM(I,KX),I=13,17)
  124   FORMAT (2/)
  125   FORMAT (F7.4,6F9.6,D9.3,4F9.6,D9.3,2F9.6,D9.3,3F9.6)
C
        RETURN 1
C
  200   CONTINUE
        READ (1,204)
        DO 210 KX=1,CHOP+1
  210   READ (1,205,END=300) DIST(KX),
     1  VM(KX),PM(KX),(CM(I,KX),I=1,11),IMPM(KX)
  204   FORMAT (5/)
  205   FORMAT (F7.4,F9.4,F9.2,4F9.6,D9.3,4F9.6,D9.3,2F9.6)
C
        READ (1,214)
        DO 220 KX=1,CHOP+1
  220   READ (1,215) DIST(KX),TA(KX),LCHM(KX),
     1    (XM(I,KX),I=1,11),XM(12,KX)
  214   FORMAT (2/)
  215   FORMAT (F7.4,F9.6,13F9.4)
C
        READ (1,224)
        DO 230 KX=1,CHOP+1
  230   READ (1,225) DIST(KX),OSMM(KX),
     1    FVM(KX),(FKM(I,KX),I=1,11),FKM(16,KX),FKM(17,KX)
  224   FORMAT (2/)
  225   FORMAT (F7.4,6F9.6,D9.3,4F9.6,D9.3,3F9.6)
C
	RETURN 1
C
  300   RETURN 2
	END
C
