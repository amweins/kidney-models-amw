	PROGRAM BCR
C
C   IT IS ASSUMED THAT THE SOLUTION IS MADE UP VIA HCL TITRATION:
C   I.E. FIRST THE HCO3 CONCENTRATION IS SPECIFIED ALONG WITH THE TOTAL
C   CONCENTRATION OF THE OTHER BUFFERS; NEXT THE BUFFER PAIRS ARE 
C   ADJUSTED; AND FINALLY THE CL CONCENTRATION IS CALCULATED FROM 
C   ELECTRONEUTRALITY.
C
	INTEGER SOLS,TLIM
C
	DOUBLE PRECISION
     1   Z(15),ZIMPM,ZIMPS,
     1   PKC,PKF,PKP,PKN,KCO2,
     1   TFM,TPM,TNM,QFM,QPM,QNM,
     1   TFS(2),TPS(2),TNS(2),QFS(2),QPS(2),QNS(2),
     1   VM,PM,CM(15),IMPM,LCHM,
     1   VS(2),PS(2),CS(15,2),IMPS(2),LCHS(2),
     1   DT,FVM
C
	INTEGER NWR,IWR(9),JWR(9)
	CHARACTER*5 VNAME1,VNAME2
	DOUBLE PRECISION VMAT(20,3),VINCR(9)
C
C
	OPEN (19,FILE='jmpstbound.tem')
	OPEN (22,FILE='jmpstbound.dat')
C
	SOLS=15
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
	ZIMPM=1.D0
	ZIMPS=0.D0
C
	PKC=3.57
	PKF=3.76
	PKP=6.8
	PKN=9.15
	KCO2=340.
C
C
	PRINT 10
   10   FORMAT(' PRINT THE NUMBER OF VARIABLES TO BE INCREMENTED')
	READ *, NWR
	IF(NWR.EQ.0) GO TO 801
	DO 800 KWR=1,NWR
   15   PRINT 20, KWR
   20   FORMAT(' PRINT THE NAME OF VARIABLE NUMBER-',I2)
	READ 30, VNAME1
   30   FORMAT(A4)
C
C
	IF(VNAME1.EQ.'fvm') THEN
	IWR(KWR)=1
	JWR(KWR)=1
	ELSE IF(VNAME1.EQ.'vm') THEN
	IWR(KWR)=2
	JWR(KWR)=1
	ELSE IF(VNAME1.EQ.'pm') THEN
	IWR(KWR)=3
	JWR(KWR)=1
	ELSE IF(VNAME1.EQ.'impm') THEN
	IWR(KWR)=4
	JWR(KWR)=1
	ELSE IF(VNAME1.EQ.'cm') THEN
	IWR(KWR)=5
	JWR(KWR)=1
	GO TO 200
	ELSE IF(VNAME1.EQ.'vs0') THEN
	IWR(KWR)=2
	JWR(KWR)=2
	ELSE IF(VNAME1.EQ.'ps0') THEN
	IWR(KWR)=3
	JWR(KWR)=2
	ELSE IF(VNAME1.EQ.'imps0') THEN
	IWR(KWR)=4
	JWR(KWR)=2
	ELSE IF(VNAME1.EQ.'cs0') THEN
	IWR(KWR)=5
	JWR(KWR)=2
	GO TO 200
	ELSE IF(VNAME1.EQ.'vs1') THEN
	IWR(KWR)=2
	JWR(KWR)=3
	ELSE IF(VNAME1.EQ.'ps1') THEN
	IWR(KWR)=3
	JWR(KWR)=3
	ELSE IF(VNAME1.EQ.'imps1') THEN
	IWR(KWR)=4
	JWR(KWR)=3
	ELSE IF(VNAME1.EQ.'cs1') THEN
	IWR(KWR)=5
	JWR(KWR)=3
	GO TO 200
	ELSE
	PRINT 35
   35   FORMAT (' NO SUCH VARIABLE, TRY AGAIN')
	GO TO 15
	ENDIF
	GO TO 700
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
	READ(19,50) DT,FVM,
     1  VM,VS(1),VS(2),PM,PS(1),PS(2),IMPM,IMPS(1),IMPS(2),
     1  (CM(I),CS(I,1),CS(I,2),I=1,11),
     1  (CM(I),CS(I,1),CS(I,2),I=13,15)
   50   FORMAT (2D12.4,/,3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
	VMAT(1,1)=FVM
	VMAT(2,1)=VM
	VMAT(3,1)=PM
	VMAT(4,1)=IMPM
	VMAT(2,2)=VS(1)
	VMAT(3,2)=PS(1)
	VMAT(4,2)=IMPS(1)
	VMAT(2,3)=VS(2)
	VMAT(3,3)=PS(2)
	VMAT(4,3)=IMPS(2)
	DO 40 I=1,SOLS
	VMAT(4+I,1)=CM(I)
	VMAT(4+I,2)=CS(I,1)
   40   VMAT(4+I,3)=CS(I,2)
C
	PRINT 60
   60   FORMAT(' TLIM=',$)
	READ *, TLIM
C
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
	DO 900 K=1,TLIM
C
C ASSUMING EQUILIBRIUM BETWEEN DISSOLVED CO2 AND H2CO3
C
	CM(5)=CM(6)/KCO2
	CS(5,1)=CS(6,1)/KCO2
	CS(5,2)=CS(6,2)/KCO2
	LCHM=PKC + DLOG10(CM(4)/CM(5))
	LCHS(1)=PKC + DLOG10(CS(4,1)/CS(5,1))
	LCHS(2)=PKC + DLOG10(CS(4,2)/CS(5,2))
	CM(12)=10**(-LCHM)
	CS(12,1)=10**(-LCHS(1))
	CS(12,2)=10**(-LCHS(2))
C
	TPM=CM(7)+CM(8)
	TNM=CM(10)+CM(11)
	TFM=CM(13)+CM(14)
	TPS(1)=CS(7,1)+CS(8,1)
	TPS(2)=CS(7,2)+CS(8,2)
	TNS(1)=CS(10,1)+CS(11,1)
	TNS(2)=CS(10,2)+CS(11,2)
	TFS(1)=CS(13,1)+CS(14,1)
	TFS(2)=CS(13,2)+CS(14,2)
C
	QPM=10.**(LCHM-PKP)
	QNM=10.**(LCHM-PKN)
	QFM=10.**(LCHM-PKF)
	QPS(1)=10.**(LCHS(1)-PKP)
	QPS(2)=10.**(LCHS(2)-PKP)
	QNS(1)=10.**(LCHS(1)-PKN)
	QNS(2)=10.**(LCHS(2)-PKN)
	QFS(1)=10.**(LCHS(1)-PKF)
	QFS(2)=10.**(LCHS(2)-PKF)
C
	CM(8)=TPM/(1.+QPM)
	CM(7)=TPM-CM(8)
	CM(11)=TNM/(1.+QNM)
	CM(10)=TNM-CM(11)
	CM(14)=TFM/(1.+QFM)
	CM(13)=TFM-CM(14)
C
	CS(8,1)=TPS(1)/(1.+QPS(1))
	CS(8,2)=TPS(2)/(1.+QPS(2))
	CS(7,1)=TPS(1)-CS(8,1)
	CS(7,2)=TPS(2)-CS(8,2)
	CS(11,1)=TNS(1)/(1.+QNS(1))
	CS(11,2)=TNS(2)/(1.+QNS(2))
	CS(10,1)=TNS(1)-CS(11,1)
	CS(10,2)=TNS(2)-CS(11,2)
	CS(14,1)=TFS(1)/(1.+QFS(1))
	CS(14,2)=TFS(2)/(1.+QFS(2))
	CS(13,1)=TFS(1)-CS(14,1)
	CS(13,2)=TFS(2)-CS(14,2)
C
	CM(3)=IMPM*ZIMPM + CM(1)*Z(1) + CM(2)*Z(2)
	CS(3,1)=IMPS(1)*ZIMPS + CS(1,1)*Z(1) + CS(2,1)*Z(2)
	CS(3,2)=IMPS(2)*ZIMPS + CS(1,2)*Z(1) + CS(2,2)*Z(2)
	DO 130 I=4,SOLS
	CM(3)=CM(3)+CM(I)*Z(I)
	CS(3,1)=CS(3,1)+CS(I,1)*Z(I)
	CS(3,2)=CS(3,2)+CS(I,2)*Z(I)
  130   CONTINUE
C
	WRITE (22,150) DT,FVM,
     1  VM,VS(1),VS(2),PM,PS(1),PS(2),IMPM,IMPS(1),IMPS(2),
     1  (CM(I),CS(I,1),CS(I,2),I=1,11),
     1  (CM(I),CS(I,1),CS(I,2),I=13,15)
  150   FORMAT (2D12.4,/,3F8.4,/,3F8.4,/,3F8.4,/,(3F14.9))
C
C
	IF(NWR.EQ.0) GO TO 145
	DO 140 KWR=1,NWR
	VMAT(IWR(KWR),JWR(KWR))=VMAT(IWR(KWR),JWR(KWR))+VINCR(KWR)
  140   CONTINUE
  145   CONTINUE
C
	FVM=VMAT(1,1)
	VM=VMAT(2,1)
	PM=VMAT(3,1)
	IMPM=VMAT(4,1)
	VS(1)=VMAT(2,2)
	PS(1)=VMAT(3,2)
	IMPS(1)=VMAT(4,2)
	VS(2)=VMAT(2,3)
	PS(2)=VMAT(3,3)
	IMPS(2)=VMAT(4,3)
	DO 160 I=1,SOLS
	CM(I)=VMAT(4+I,1)
	CS(I,1)=VMAT(4+I,2)
  160   CS(I,2)=VMAT(4+I,3)
C
C
  900   CONTINUE
	STOP
	END
