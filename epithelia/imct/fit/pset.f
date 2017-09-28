	SUBROUTINE PSET(TAU,TLIM)
C PROGRAM TO READ AND INITIALIZE THE PARAMETERS FOR COMP
C
	INTEGER SOLS,TAU,T,TLIM
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY(4),KDHY(4),L0,L(2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   DUMVM,VM,PM,CM(12),IMPM,LCHM,XM(12),
     1   VS,PS,CS(12),IMPS,LCHS,XS(12)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(2),MUA,CHVL0,CHVL(2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(2),PE(2),CE(12,2),LCHE,XE(12),
     1   FEVM(2),FEKM(12,2),FEVS(2),FEKS(12,2),CURE
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,2),HCBUF(3,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12),ATIS(3,12),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),
     1   VI(3,2),PI(3,2),CI(3,12,2),IMP(3),LCHI(3),XI(3,12),
     1   FIVM(3,2),FIKM(3,12,2),FIVS(3,2),FIKS(3,12,2),CURI(3),
     1   JV(3,2),JK(3,12,2)
C
        COMMON SOLS,T,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   DUMVM,VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AI0,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
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
C
	RT=0.193D+05
	RTE=2.57
	F=0.965000D+05
	PKC=3.57
	PKP=6.8
	PKN=9.15
C
C 7/29/91- IN THE OLDER VERSIONS THE DIAGONAL LMI AND LIS WERE FIRST READ
C AND THEN THE COUPLING COEFFICIENTS.  NOW EACH COTRANSPORTER IS READ IN
C COEFFICIENT BY COEFFICIENT.  THIS ALLOWS MULTIPLE TRANSPORTERS TO SHARE
C A SINGLE SPECIES WITHOUT HAVING TO BACK-CALCULATE THE COUPLING COEFFICIENTS.
C
	DO 18 K=1,1
	DO 18 I=1,SOLS
	ATMI(K,I)=0.
	ATIS(K,I)=0.
	DO 18 J=1,SOLS
	LMI(K,I,J)=0.D0
   18   LIS(K,I,J)=0.D0
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (20,30) TAU,TLIM,EPSI,
     1  AME,AE0,MUA,KHY(4),KDHY(4),
     1  L0,CHVL0,MUV,LPME,LPES,
     1  (SME(I),SES(I),HME(I),HES(I),I=1,SOLS)
   30   FORMAT (2I5,D12.4,/,5D12.4,/,5D12.4,/,(2F6.3, 2D12.4))
C
	DO 50 K=1,1
	READ (20,25)
     1  AIE(K),AI0(K),CLVL0(K),IMP0(K),ZIMP(K),
     1  KHY(K),KDHY(K),TBUF(K),PKB(K),
     1  NP(K),KNH4(K),NPHK(K),
     1  LHP(K),XIHP(K),XHP(K)
	READ (20,26)
     1  LPMI(K),LPIS(K),
     1  (SMI(K,I),SIS(K,I),HMI(K,I),HIS(K,I),I=1,SOLS)
   25   FORMAT (5D12.4,/,4D12.4,/,3D12.4,/,3D12.4)
   26   FORMAT (2D12.4,/,(2F6.3, 2D12.4))
C
	READ(20,31) NMI,NIS
	IF(NMI.EQ.0) GO TO 36
	DO 44 LL=1,NMI
	READ(20,31) I,J,QMI
	LMI(K,I,J)=LMI(K,I,J)+QMI
   44   LMI(K,J,I)=LMI(K,I,J)
C
   36   IF(NIS.EQ.0) GO TO 50
	DO 46 LL=1,NIS
	READ(20,31) I,J,QIS
	LIS(K,I,J)=LIS(K,I,J)+QIS
   46   LIS(K,J,I)=LIS(K,I,J)
   31   FORMAT (2I5,D12.4)
C
   50   CONTINUE
C
	RETURN
	END
