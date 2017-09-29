C CPBCR-PROGRAM TO CREATE THE FILE CBOUND.DAT TO BE READ BY CAP.FOR
C   IT IS ASSUMED THAT THE SOLUTION IS MADE UP VIA HCL TITRATION:
C   I.E. FIRST THE HCO3 CONCENTRATION IS SPECIFIED ALONG WITH THE TOTAL
C   CONCENTRATION OF THE OTHER BUFFERS; NEXT THE BUFFER PAIRS ARE 
C   ADJUSTED; AND FINALLY THE CL CONCENTRATION IS CALCULATED FROM 
C   ELECTRONEUTRALITY.
C
	INTEGER SOLS,TLIM,SLICE
C
	DOUBLE PRECISION
     1   Z(15),ZIMPC,ZIMPS,
     1   PKC,PKF,PKP,PKN,KCO2,
     1   TFC,TPC,TNC,QFC,QPC,QNC,
     1   TFS(2),TPS(2),TNS(2),QFS(2),QPS(2),QNS(2),
     1   VC,PC,CC(16),IMPC,LCHC,
     1   VS(2),PS(2),CS(16,2),IMPS(2),LCHS(2),
     1   DT,FVC,POX(20),HCTC
C
C
	INTEGER NWR,IWR(9),JWR(9)
	CHARACTER*5 VNAME1,VNAME2
	DOUBLE PRECISION VMAT(20,4),VINCR(9)
C
C
	OPEN (UNIT=19,FILE='oidbound.tem',STATUS='old')
	OPEN (UNIT=22,FILE='oidbound.dat')
C
C
C  SOLUTE INDEX:
C	1-	NA+
C	2-	K+
C	3-	CL-
C	4-	HCO3-
C	5-	H2CO3
C	6-	CO2
C	7-	HPO4--
C	8-	H2PO4-
C	9-	UREA
C	10-	NH3
C	11-	NH4+
C	12-	H+
C	13-	HCO2-
C	14-	H2CO2
C	15-	GLUC
C
	Z(1)=1.D0
	Z(2)=1.D0
	Z(3)=-1.D0
	Z(4)=-1.D0
	Z(5)=0.D0
	Z(6)=0.D0
	Z(7)=-2.D0
	Z(8)=-1.D0
	Z(9)=0.D0
	Z(10)=0.D0
	Z(11)=1.D0
	Z(12)=1.D0
	Z(13)=-1.D0
	Z(14)=0.D0
	Z(15)=0.D0
	ZIMPC=-1.D0
	ZIMPS=0.D0
C
	PKC=3.57
	PKF=3.76
	PKP=6.8
	PKN=9.15
	KCO2=340.
C
	SOLS=15
	EPSI=1.D-12
        SLICE=2
C
C
	PRINT 10
   10   FORMAT(' TYPE THE NUMBER OF VARIABLES TO BE INCREMENTED')
	READ *, NWR
	IF(NWR.EQ.0) GO TO 801
	DO 800 KWR=1,NWR
   15  PRINT 20, KWR
   20   FORMAT(' TYPE THE NAME OF VARIABLE NUMBER-',I2)
	READ 30, VNAME1
   30   FORMAT(A4)
C
	IF(VNAME1.EQ.'fvc') THEN
	IWR(KWR)=1
	JWR(KWR)=1
	GO TO 700
	ELSE IF(VNAME1.EQ.'hct') THEN
        IWR(KWR)=1
	JWR(KWR)=2
	GO TO 700
	ELSE IF(VNAME1.EQ.'pox') THEN
        IWR(KWR)=1
	JWR(KWR)=4
        PRINT 31
   31   FORMAT(' PO2 AT SLICE#(0-2): ',$)
        READ *, ISLICE
        IWR(KWR)=1+ISLICE
	GO TO 700
	ELSE IF(VNAME1.EQ.'vc') THEN
	IWR(KWR)=2
	JWR(KWR)=1
	GO TO 700
	ELSE IF(VNAME1.EQ.'pc') THEN
	IWR(KWR)=3
	JWR(KWR)=1
	GO TO 700
	ELSE IF(VNAME1.EQ.'impc') THEN
	IWR(KWR)=4
	JWR(KWR)=1
	GO TO 700
	ELSE IF(VNAME1.EQ.'cc') THEN
	IWR(KWR)=5
	JWR(KWR)=1
	GO TO 200
	ELSE IF(VNAME1.EQ.'vs0') THEN
	IWR(KWR)=2
	JWR(KWR)=2
	GO TO 700
	ELSE IF(VNAME1.EQ.'ps0') THEN
	IWR(KWR)=3
	JWR(KWR)=2
	GO TO 700
	ELSE IF(VNAME1.EQ.'imps0') THEN
	IWR(KWR)=4
	JWR(KWR)=2
	GO TO 700
	ELSE IF(VNAME1.EQ.'cs0') THEN
	IWR(KWR)=5
	JWR(KWR)=2
	GO TO 200
	ELSE IF(VNAME1.EQ.'vs1') THEN
	IWR(KWR)=2
	JWR(KWR)=3
	GO TO 700
	ELSE IF(VNAME1.EQ.'ps1') THEN
	IWR(KWR)=3
	JWR(KWR)=3
	GO TO 700
	ELSE IF(VNAME1.EQ.'imps1') THEN
	IWR(KWR)=4
	JWR(KWR)=3
	GO TO 700
	ELSE IF(VNAME1.EQ.'cs1') THEN
	IWR(KWR)=5
	JWR(KWR)=3
	GO TO 200
	ELSE
	PRINT 35
   35   FORMAT (' NO SUCH VARIABLE, TRY AGAIN')
	GO TO 15
	ENDIF
C
C
  200   PRINT 201
  201   FORMAT(' PRINT THE SPECIES')
	READ 30, VNAME2
C
	IF(VNAME2.EQ.'na') THEN
	 IWR(KWR)=IWR(KWR)+0
	ELSE IF(VNAME2.EQ.'k') THEN
	 IWR(KWR)=IWR(KWR)+1
	ELSE IF(VNAME2.EQ.'cl') THEN
	 IWR(KWR)=IWR(KWR)+2
	ELSE IF(VNAME2.EQ.'hco3') THEN
	 IWR(KWR)=IWR(KWR)+3
	ELSE IF(VNAME2.EQ.'h2co3') THEN
	 IWR(KWR)=IWR(KWR)+4
	ELSE IF(VNAME2.EQ.'co2') THEN
	 IWR(KWR)=IWR(KWR)+5
	ELSE IF(VNAME2.EQ.'hpo4') THEN
	 IWR(KWR)=IWR(KWR)+6
	ELSE IF(VNAME2.EQ.'h2po4') THEN
	 IWR(KWR)=IWR(KWR)+7
	ELSE IF(VNAME2.EQ.'urea') THEN
	 IWR(KWR)=IWR(KWR)+8
	ELSE IF(VNAME2.EQ.'nh3') THEN
	 IWR(KWR)=IWR(KWR)+9
	ELSE IF(VNAME2.EQ.'nh4') THEN
	 IWR(KWR)=IWR(KWR)+10
	ELSE IF(VNAME2.EQ.'hco2') THEN
	 IWR(KWR)=IWR(KWR)+12
	ELSE IF(VNAME2.EQ.'h2co2') THEN
	 IWR(KWR)=IWR(KWR)+13
	ELSE IF(VNAME2.EQ.'gluc') THEN
	 IWR(KWR)=IWR(KWR)+14
	ELSE
	PRINT 235
  235   FORMAT (' NO SUCH SPECIES, TRY AGAIN')
	GO TO 200
	ENDIF
C
C
  700   PRINT 710
  710   FORMAT(' VARIABLE INCREMENT=',$)
	READ *, VINCR(KWR)
C
  800   CONTINUE
  801   CONTINUE
C
C
C BOUNDARY VALUES,ACTIVE TRANSPORT,AND THE TIME INCREMENT ARE
C READ. THE TIME STEP WILL CYCLE FROM THIS POINT.
C
	READ(19,50) DT,FVC,HCTC,(POX(I),I=1,SLICE+1)
        READ(19,51) VC,VS(1),VS(2),
     1  PC,PS(1),PS(2),
     1  CC(16),CS(16,1),CS(16,2),
     1  (CC(I),CS(I,1),CS(I,2),I=1,11),
     1  (CC(I),CS(I,1),CS(I,2),I=13,15)
   50   FORMAT (2D12.4,F8.4,/,9F8.4)
   51   FORMAT (3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
	VMAT(1,1)=FVC
	VMAT(1,2)=HCTC
        DO 39 I=1,SLICE+1
   39   VMAT(I,4)=POX(I)
	VMAT(2,1)=VC
	VMAT(3,1)=PC
	VMAT(4,1)=IMPC
	VMAT(2,2)=VS(1)
	VMAT(3,2)=PS(1)
	VMAT(4,2)=IMPS(1)
	VMAT(2,3)=VS(2)
	VMAT(3,3)=PS(2)
	VMAT(4,3)=IMPS(2)
	DO 40 I=1,SOLS
	VMAT(4+I,1)=CC(I)
	VMAT(4+I,2)=CS(I,1)
   40   VMAT(4+I,3)=CS(I,2)
C
	PRINT 60
   60   FORMAT(' TLIM=',$)
	READ *, TLIM
C
	DO 900 K=1,TLIM
C
C ASSUMING EQUILIBRIUM BETWEEN DISSOLVED CO2 AND H2CO3
C
	CC(5)=CC(6)/KCO2
	CS(5,1)=CS(6,1)/KCO2
	CS(5,2)=CS(6,2)/KCO2
	LCHC=PKC + DLOG10(CC(4)/CC(5))
	LCHS(1)=PKC + DLOG10(CS(4,1)/CS(5,1))
	LCHS(2)=PKC + DLOG10(CS(4,2)/CS(5,2))
	CC(12)=10**(-LCHC)
	CS(12,1)=10**(-LCHS(1))
	CS(12,2)=10**(-LCHS(2))
C
	TPC=CC(7)+CC(8)
	TNC=CC(10)+CC(11)
	TFC=CC(13)+CC(14)
	TPS(1)=CS(7,1)+CS(8,1)
	TPS(2)=CS(7,2)+CS(8,2)
	TNS(1)=CS(10,1)+CS(11,1)
	TNS(2)=CS(10,2)+CS(11,2)
	TFS(1)=CS(13,1)+CS(14,1)
	TFS(2)=CS(13,2)+CS(14,2)
C
	QPC=10.**(LCHC-PKP)
	QNC=10.**(LCHC-PKN)
	QFC=10.**(LCHC-PKF)
	QPS(1)=10.**(LCHS(1)-PKP)
	QPS(2)=10.**(LCHS(2)-PKP)
	QNS(1)=10.**(LCHS(1)-PKN)
	QNS(2)=10.**(LCHS(2)-PKN)
	QFS(1)=10.**(LCHS(1)-PKF)
	QFS(2)=10.**(LCHS(2)-PKF)
C
	CC(8)=TPC/(1.+QPC)
	CC(7)=TPC-CC(8)
	CC(11)=TNC/(1.+QNC)
	CC(10)=TNC-CC(11)
	CC(14)=TFC/(1.+QFC)
	CC(13)=TFC-CC(14)
C
	CS(8,1)=TPS(1)/(1.+QPS(1))
	CS(7,1)=TPS(1)-CS(8,1)
	CS(11,1)=TNS(1)/(1.+QNS(1))
	CS(10,1)=TNS(1)-CS(11,1)
	CS(14,1)=TFS(1)/(1.+QFS(1))
	CS(13,1)=TFS(1)-CS(14,1)
C
	CS(8,2)=TPS(2)/(1.+QPS(2))
	CS(7,2)=TPS(2)-CS(8,2)
	CS(11,2)=TNS(2)/(1.+QNS(2))
	CS(10,2)=TNS(2)-CS(11,2)
	CS(14,2)=TFS(2)/(1.+QFS(2))
	CS(13,2)=TFS(2)-CS(14,2)
C
	CC(3)=IMPC*ZIMPC + CC(1)*Z(1) + CC(2)*Z(2)
	CS(3,1)=IMPS(1)*ZIMPS + CS(1,1)*Z(1) + CS(2,1)*Z(2)
	CS(3,2)=IMPS(2)*ZIMPS + CS(1,2)*Z(1) + CS(2,2)*Z(2)
	DO 130 I=4,SOLS
	CC(3)=CC(3)+CC(I)*Z(I)
	CS(3,1)=CS(3,1)+CS(I,1)*Z(I)
	CS(3,2)=CS(3,2)+CS(I,2)*Z(I)
  130   CONTINUE
C
C
	WRITE (22,150) DT,FVC,HCTC,(POX(I),I=1,SLICE+1)
        WRITE (22,151) VC,VS(1),VS(2),
     1  PC,PS(1),PS(2),
     1  CC(16),CS(16,1),CS(16,2),
     1  (CC(I),CS(I,1),CS(I,2),I=1,11),
     1  (CC(I),CS(I,1),CS(I,2),I=13,15)
  150   FORMAT (2D12.4,F8.4,/,9F8.4)
  151   FORMAT (3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
C
	IF(NWR.EQ.0) GO TO 145
	DO 140 KWR=1,NWR
	VMAT(IWR(KWR),JWR(KWR))=VMAT(IWR(KWR),JWR(KWR))+VINCR(KWR)
  140   CONTINUE
  145   CONTINUE
C
C
	FVC=VMAT(1,1)
	HCTC=VMAT(1,2)
        DO 159 I=1,SLICE+1
  159   POX(I)=VMAT(I,4)
	VC=VMAT(2,1)
	PC=VMAT(3,1)
	IMPC=VMAT(4,1)
	VS(1)=VMAT(2,2)
	PS(1)=VMAT(3,2)
	IMPS(1)=VMAT(4,2)
	VS(2)=VMAT(2,3)
	PS(2)=VMAT(3,3)
	IMPS(2)=VMAT(4,3)
	DO 160 I=1,SOLS
	CC(I)=VMAT(4+I,1)
	CS(I,1)=VMAT(4+I,2)
  160   CS(I,2)=VMAT(4+I,3)
C
  900   CONTINUE
	STOP
	END
