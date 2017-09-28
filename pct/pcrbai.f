	PROGRAM PCR
C PROGRAM TO CREATE THE PARAM.DAT FILE FOR EXECUTION OF COMP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE PARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,NMI(3),NIS(3),II(3,50),JJ(3,50)
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   EPSI,KHY(4),KDHY(4),L0,
     1   AME,AE0,MUA,CHVL0,MUV,
     1   LPME,LPES,SME(15),SES(15),HME(15),HES(15),
     1   AIE(3),AMI(3),AIS(3),CLVL0(3),IMP0(3),
     1   ZIMP(3),TBUF(3),PKB(3),PMI(3,50),PIS(3,50),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),HMI(3,15),HIS(3,15),
     1   NP(3),KNH4(3),NPHK(3),LHP(3),XIHP(3),XHP(3),
     1   NAE1(3),NTSC(3),NNHE3(3),QIAMM
C
	double precision dlpme
C
	OPEN (19,FILE='param.tem')
	OPEN (20,FILE='param.dat')
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
C	15-	GLUCOSE
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
C LIE WILL CONTAIN THE CONDUCTANCE COMPONENT OF LIS
C 7/29/91- IN THE OLDER VERSIONS THE DIAGONAL LMI AND LIS WERE FIRST READ
C AND THEN THE COUPLING COEFFICIENTS.  NOW EACH COTRANSPORTER IS READ IN
C COEFFICIENT BY COEFFICIENT.  THIS ALLOWS MULTIPLE TRANSPORTERS TO SHARE
C A SINGLE SPECIES WITHOUT HAVING TO BACK-CALCULATE THE COUPLING COEFFICIENTS.
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,EPSI,
     1  AME,AE0,MUA,KHY(4),KDHY(4),
     1  L0,CHVL0,MUV,LPME,LPES,
     1  (SME(I),SES(I),HME(I),HES(I),I=1,SOLS)
   30   FORMAT (2I5,D12.4,/,5D12.4,/,5D12.4,/,(2F6.3, 2D12.4))
C
	DO 50 K=1,1
	READ (19,25)
     1  AIE(K),AMI(K),AIS(K),
     1  CLVL0(K),IMP0(K),ZIMP(K),
     1  KHY(K),KDHY(K),TBUF(K),PKB(K),QIAMM,
     1  NP(K),KNH4(K),NPHK(K),
     1  LHP(K),XIHP(K),XHP(K),
     1  NAE1(K),NTSC(K),NNHE3(K)
	READ (19,26)
     1  LPMI(K),LPIS(K),
     1  (SMI(K,I),SIS(K,I),HMI(K,I),HIS(K,I),I=1,SOLS)
   25   FORMAT (3D12.4,/,3D12.4,/,5D12.4,/,(3D12.4))
   26   FORMAT (2D12.4,/,(2F6.3, 2D12.4))
C
	READ(19,31) NMI(K),NIS(K)
	IF(NMI(K).EQ.0) GO TO 36
	READ(19,31) (II(K,LL),JJ(K,LL),PMI(K,LL),LL=1,NMI(K))
   36   IF(NIS(K).EQ.0) GO TO 50
	READ(19,31) (II(K,LL),JJ(K,LL),PIS(K,LL),
     1   LL=NMI(K)+1,NMI(K)+NIS(K))
   31   FORMAT (2I5,D12.4)
C
   50   CONTINUE
C
	PRINT 60
   60   FORMAT(' NUMBER OF PARAMETER SETS NEEDED=',$)
	READ *, NSETS
C
	dlpme=lpme/float(nsets)
C
	DO 200 KSETS=1,NSETS+1
C
	lpme=float(nsets+1-ksets)*dlpme
C
	do 200 kimp=1,27
C
	WRITE (20,130) TAU,TLIM,EPSI,
     1  AME,AE0,MUA,KHY(4),KDHY(4),
     1  L0,CHVL0,MUV,LPME,LPES,
     1  (SME(I),SES(I),HME(I),HES(I),I=1,SOLS)
  130   FORMAT (2I5,D12.4,/,5D12.4,/,5D12.4,/,(2F6.3, 2D12.4))
C
	DO 150 K=1,1
	WRITE (20,125)
     1  AIE(K),AMI(K),AIS(K),
     1  CLVL0(K),IMP0(K),ZIMP(K),
     1  KHY(K),KDHY(K),TBUF(K),PKB(K),QIAMM,
     1  NP(K),KNH4(K),NPHK(K),
     1  LHP(K),XIHP(K),XHP(K),
     1  NAE1(K),NTSC(K),NNHE3(K)
	WRITE (20,126)
     1  LPMI(K),LPIS(K),
     1  (SMI(K,I),SIS(K,I),HMI(K,I),HIS(K,I),I=1,SOLS)
  125   FORMAT (3D12.4,/,3D12.4,/,5D12.4,/,(3D12.4))
  126   FORMAT (2D12.4,/,(2F6.3, 2D12.4))
C
	WRITE (20,131) NMI(K),NIS(K)
	IF(NMI(K).EQ.0) GO TO 136
	WRITE (20,131) (II(K,LL),JJ(K,LL),PMI(K,LL),LL=1,NMI(K))
  136   IF(NIS(K).EQ.0) GO TO 150
	WRITE (20,131) (II(K,LL),JJ(K,LL),PIS(K,LL),
     1   LL=NMI(K)+1,NMI(K)+NIS(K))
  131   FORMAT (2I5,D12.4)
C
  150   CONTINUE
C
  200   CONTINUE
	STOP
	END
