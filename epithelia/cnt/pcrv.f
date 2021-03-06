	PROGRAM PCRV
C PROGRAM TO CREATE THE PARAM.DAT FILE FOR EXECUTION OF COMP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE PARAM.TEM 
C A SELECTED PARAMETER IS INCREMENTED, AND OUTPUT TO AN ABCISSA FILE
C
	INTEGER SOLS,TAU,TLIM,NMI(3),NIS(3),II(3,50),JJ(3,50)
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   EPSI,KHY(4),KDHY(4),L0,
     1   AME,AE0,MUA,CHVL0,MUV,
     1   LPME,LPES,SME(12),SES(12),HME(12),HES(12),
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),
     1   ZIMP(3),TBUF(3),PKB(3),PMI(3,50),PIS(3,50),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),HMI(3,12),HIS(3,12),
     1   NP(3),KNH4(3),NPHK(3),LHP(3),XIHP(3),XHP(3),NAE1(3)
C
 	DOUBLE PRECISION DPAR,DPAR2,DPAR3,DPAR4,DPAR5,
     1    MPAR,PAR,PAR2,PAR3,PAR4,PAR5
C
C	EQUIVALENCE (PAR,hme(2))
C	EQUIVALENCE (PAR,lhp(2))
C	EQUIVALENCE (PAR,pmi(1,1)),(PAR2,pmi(1,2)),(PAR3,pmi(1,3))
C     1    (PAR4,his(1,3)),(PAR5,his(1,4))
C	EQUIVALENCE (PAR,AME)
C 	EQUIVALENCE (PAR,hmi(1,2))
C	EQUIVALENCE (PAR,hmi(1,3))
C	EQUIVALENCE (PAR,hmi(1,2)),(PAR2,his(1,2))
C	EQUIVALENCE (PAR,hmi(1,1))
C	EQUIVALENCE (PAR,hmi(1,1)),(PAR2,np(1))
C 	EQUIVALENCE (PAR,hmi(1,1)),(PAR2,np(1)),(PAR3,hmi(1,2))
C 	EQUIVALENCE (PAR,hmi(1,1)),(PAR2,np(1)),(PAR3,hmi(1,2)),(PAR4,AME)
C 	EQUIVALENCE (PAR,hmi(1,1)),(PAR2,hmi(1,2))
 	EQUIVALENCE (PAR,hmi(1,1)),(PAR2,hmi(1,2)),(PAR3,AME)
C	EQUIVALENCE (PAR,lpmi(1)),(PAR2,lpis(1))
C
C	EQUIVALENCE (PAR,his(2,3)),(PAR2,his(2,4))
C	EQUIVALENCE (PAR,his(2,4))
C	EQUIVALENCE (PAR,his(1,2))
C	EQUIVALENCE (PAR,his(1,3))
C	EQUIVALENCE (PAR,khy(2)),(PAR2,kdhy(2))
C	EQUIVALENCE (PAR,nphk(2)),(PAR2,his(2,2))
C
	OPEN (19,FILE='param.tem')
	OPEN (20,FILE='param.dat')
 	OPEN (21,FILE='abcis.dat')
C
	SOLS=12
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
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
C LIE WILL CONTAIN THE CONDUCTANCE COMPONENT OF LIS
C 7/29/91- IN THE OLDER VERSIONS THE DIAGONAL LMI AND LIS WERE FIRST READ
C AND THEN THE COUPLING COEFFICIENTS.  NOW EACH COTRANSPORTER IS READ IN
C COEFFICIENT BY COEFFICIENT.  THIS ALLOWS MULTIPLE TRANSPORTERS TO SHARE
C A SINGLE SPECIES WITHOUT HAVING TO BACK-CALCULATE THE COUPLING COEFFICIENTS.
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,EPSI,
     1  AME,AE0,MUA,KHY(4),KDHY(4),
     1  L0,CHVL0,MUV,LPME,LPES,
     1  (SME(I),SES(I),HME(I),HES(I),I=1,SOLS)
   30   FORMAT (2I5,D12.4,/,5D12.4,/,5D12.4,/,(2F6.3, 2D12.4))
C
	DO 50 K=1,3
	READ (19,25)
     1  AIE(K),AI0(K),CLVL0(K),IMP0(K),ZIMP(K),
     1  KHY(K),KDHY(K),TBUF(K),PKB(K),
     1  NP(K),KNH4(K),NPHK(K),
     1  LHP(K),XIHP(K),XHP(K),NAE1(K)
	READ (19,26)
     1  LPMI(K),LPIS(K),
     1  (SMI(K,I),SIS(K,I),HMI(K,I),HIS(K,I),I=1,SOLS)
   25   FORMAT (5D12.4,/,4D12.4,/,3D12.4,/,4D12.4)
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
 	NSETS=25
C
	GO TO 100
C This section for the CNT manuscript
C    Tight junction area
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	MPAR=5
	PAR=5.0*PAR
	GO TO 110
C    ROMK
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	MPAR=5
	PAR=5.0*PAR
	GO TO 110
C    Luminal Na permeability plus peritubular Na,K-ATPase 
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	DPAR2=10**(-1.0/FLOAT(NSETS))
	MPAR=5.0
	PAR=5.0*PAR
	PAR2=SQRT(10.)*PAR2
	GO TO 110
C    Luminal Na permeability, peritubular Na,K-ATPase, plus ROMK
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	DPAR2=10**(-1.0/FLOAT(NSETS))
	DPAR3=DPAR
	MPAR=5.0
	PAR=5.0*PAR
	PAR2=SQRT(10.)*PAR2
	PAR3=5.0*PAR3
	GO TO 110
C    Luminal Na permeability, peritubular Na,K-ATPase, ROMK, AME
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	DPAR2=10**(-1.0/FLOAT(NSETS))
	DPAR3=DPAR
	DPAR4=1./DPAR
	MPAR=5.0
	PAR=5.0*PAR
	PAR2=SQRT(10.)*PAR2
	PAR3=5.0*PAR3
	PAR4=PAR4/5.0
	GO TO 110
C    Luminal Na permeability plus ROMK
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	DPAR2=DPAR
	MPAR=5.0
	PAR=5.0*PAR
	PAR2=5.0*PAR2
	GO TO 110
C    Luminal Na permeability, ROMK, AME
  100   CONTINUE
	DPAR=10**(-ALOG10(25.)/FLOAT(NSETS))
	DPAR2=DPAR
	DPAR3=1./DPAR
	MPAR=5.0
	PAR=5.0*PAR
	PAR2=5.0*PAR2
	PAR3=PAR3/5.0
	GO TO 110
C
C This section for the K-perm and H,K-ATpase variation
C   also for the Na-perm and Na,K-ATpase variation
C 	DPAR=10**(-2.0/FLOAT(NSETS))
C 	DPAR2=10**(-1.2/FLOAT(NSETS))
C 	MPAR=0.3D+1
C 	PAR=PAR*0.3D+1
C 	PAR2=PAR2*0.2D+1
C
C This section for simultaneous luminal and peritubular K-perm variation
C	DPAR=10**(+1.0/FLOAT(NSETS))
C	DPAR2=10**(-1.0/FLOAT(NSETS))
C	MPAR=0.1D+1
C
C This section for water permeabilities
	DPAR=10**(ALOG10(25.)/FLOAT(NSETS))
	DPAR2=DPAR
	MPAR=1.0
	GO TO 110
C This section for luminal Na and K permeabilities plus
C   peritubular Na,K-ATPase and K permeability
	DPAR=10**(-2.0/FLOAT(NSETS))
	DPAR2=10**(-0.95/FLOAT(NSETS))
	DPAR3=10**(-2.0/FLOAT(NSETS))
	MPAR=0.1D+2
	PAR=PAR*0.1D+2
	PAR2=PAR2*0.3D+1
	PAR3=PAR3*0.1D+2
C This section for all other parameter variations
	GO TO 110
	DPAR=10**(-2.0/FLOAT(NSETS))
	DPAR2=10**(-2.0/FLOAT(NSETS))
	DPAR3=10**(-2.0/FLOAT(NSETS))
	DPAR4=10**(-2.0/FLOAT(NSETS))
	MPAR=0.1D+2
	PAR=PAR*0.1D+2
	PAR2=PAR2*0.1D+2
	PAR3=PAR3*0.1D+2
	PAR4=PAR4*0.1D+2
  110   CONTINUE
C
 	DO 200 KSETS=1,NSETS+1
C
	WRITE (20,130) TAU,TLIM,EPSI,
     1  AME,AE0,MUA,KHY(4),KDHY(4),
     1  L0,CHVL0,MUV,LPME,LPES,
     1  (SME(I),SES(I),HME(I),HES(I),I=1,SOLS)
  130   FORMAT (2I5,D12.4,/,5D12.4,/,5D12.4,/,(2F6.3, 2D12.4))
C
	DO 150 K=1,3
	WRITE (20,125)
     1  AIE(K),AI0(K),CLVL0(K),IMP0(K),ZIMP(K),
     1  KHY(K),KDHY(K),TBUF(K),PKB(K),
     1  NP(K),KNH4(K),NPHK(K),
     1  LHP(K),XIHP(K),XHP(K),NAE1(K)
	WRITE (20,126)
     1  LPMI(K),LPIS(K),
     1  (SMI(K,I),SIS(K,I),HMI(K,I),HIS(K,I),I=1,SOLS)
  125   FORMAT (5D12.4,/,4D12.4,/,3D12.4,/,4D12.4)
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
C	WRITE (21,155) DLOG10(MPAR)
 	WRITE (21,155) MPAR
  155   FORMAT (E16.8)
 	MPAR=MPAR*DPAR
 	PAR=PAR*DPAR
 	PAR2=PAR2*DPAR2
 	PAR3=PAR3*DPAR3
 	PAR4=PAR4*DPAR4
 	PAR5=PAR5*DPAR5
C
  200   CONTINUE
	STOP
	END
