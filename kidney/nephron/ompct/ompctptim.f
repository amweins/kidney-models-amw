	SUBROUTINE OMPCTPTIM(NPAR)
C PROGRAM TO READ AND INITIALIZE THE PARAMETERS FOR CCT
C
	INTEGER SOLS,T,X,NPAR,CHOP
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(15),RT,RTE,F,DIST,
     1   EPSI,DT,RTAU,TIME,L0,L(1601,2),
     1   PKC,PKF,PKN,PKP,KHY(5),KDHY(5),BCO2(3)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(1601,2),PM(1601,2),CM(15,1601,2),
     1   IMPM(1601,2),LCHM(1601),XM(15,1601),
     1   VS(1601,2),PS(1601,2),CS(15,1601,2),
     1   IMPS(1601,2),LCHS(1601),XS(15,1601),
     1   TL,DX,RM0,MUM,ETA,ZIMPS,
     1   SM(1601,2),AM(1601,2),FVM(1601,2),FKM(16,1601,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(1601,2),MUA,CHVL0,CHVL(1601,2),MUV,
     1   LPME,LPES,SME(15),SES(15),
     1   HME(15),HES(15),CME(15),CES(15),
     1   VE(1601,2),PE(1601,2),CE(15,1601,2),LCHE(1601),XE(15,1601),
     1   FEVM(1601,2),FEKM(15,1601,2),FEVS(1601,2),FEKS(15,1601,2),
     1   CURE(1601)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AMI(3),AIS(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,1601,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,1601,2),HCBUF(3,1601,2),
     1   LPMI(3),LPIS(3),SMI(3,15),SIS(3,15),
     1   HMI(3,15),HIS(3,15),CMI(3,15),CIE(3,15),CIS(3,15),
     1   LMI(3,15,15),LIS(3,15,15),
     1   ATMI(3,15,1601),ATIS(3,15,1601),ATIE(3,15,1601),
     1   VI(3,1601,2),PI(3,1601,2),CI(3,15,1601,2),IMP(3,1601),
     1   LCHI(3,1601),XI(3,15,1601),FIVM(3,1601,2),FIKM(3,15,1601,2),
     1   FIVS(3,1601,2),FIKS(3,15,1601,2),CURI(3,1601),JV(3,1601,2),
     1   JK(3,15,1601,2)
C SPECIAL TRANSPORTERS
	CHARACTER*1 ISOFM
	DOUBLE PRECISION
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),NTSC(3),NNHE3(3),
     1   JNAK(3,3,1601,2),JHK(3,1601,2),JHP(3,1601,2),
     1   JAE1(3,1601,2),JTSC(3,1601,2),JNHE3(3,3,1601,2),QIAMM,
     1   NNKCC(3),NKCL(3),JNKCC(3,4,1601,2),JKCC(3,3,1601,2)
C
C
        COMMON SOLS,T,X,CHOP,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKF,PKN,PKP,KHY,KDHY,BCO2,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   TL,DX,RM0,MUM,ETA,ZIMPS,
     1   SM,AM,FVM,FKM,
     1   AME,AE0,AE,MUA,CHVL0,CHVL,MUV,
     1   LPME,LPES,SME,SES,
     1   HME,HES,CME,CES,
     1   VE,PE,CE,LCHE,XE,
     1   FEVM,FEKM,FEVS,FEKS,CURE,
     1   AIE,AMI,AIS,AI0,CLVL0,IMP0,CLVL,
     1   ZIMP,TBUF,PKB,CBUF,HCBUF,
     1   LPMI,LPIS,SMI,SIS,
     1   HMI,HIS,CMI,CIE,CIS,
     1   LMI,LIS,ATMI,ATIS,ATIE,
     1   VI,PI,CI,IMP,LCHI,XI,
     1   FIVM,FIKM,FIVS,FIKS,CURI,
     1   JV,JK
	COMMON/KINET/ ISOFM,
     1   NP,KNPN,KNPK,KNH4,NPHK,
     1   LHP,XIHP,XHP,NAE1,NTSC,NNHE3,
     1   JNAK,JHK,JHP,JAE1,JTSC,JNHE3,QIAMM,
     1   NNKCC,NKCL,JNKCC,JKCC
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
	CHARACTER*5 PNAM
C
C
	GO TO (100,200) T
C
  100   READ (36,110) NPAR
  110   FORMAT(I3)
	RETURN
C
  200   IF (NPAR.EQ.0) RETURN
	DO 300 M=1,NPAR
	READ (36,210) PNAM,KPAR,IPAR,JPAR,PVAL
  210   FORMAT(A5,3I3,D12.4)
C
	IF(PNAM.EQ.'lpme') THEN
	LPME=PVAL
	ELSE IF(PNAM.EQ.'sme') THEN
	SME(IPAR)=PVAL
	ELSE IF(PNAM.EQ.'hme') THEN
	HME(IPAR)=PVAL
	ELSE IF(PNAM.EQ.'lpmi') THEN
	LPMI(KPAR)=PVAL
	ELSE IF(PNAM.EQ.'smi') THEN
	SMI(KPAR,IPAR)=PVAL
	ELSE IF(PNAM.EQ.'hmi') THEN
	HMI(KPAR,IPAR)=PVAL
	ELSE IF(PNAM.EQ.'lmi') THEN
	LMI(KPAR,IPAR,JPAR)=PVAL
	ELSE IF(PNAM.EQ.'lpis') THEN
	LPIS(KPAR)=PVAL
	ELSE IF(PNAM.EQ.'sis') THEN
	SIS(KPAR,IPAR)=PVAL
	ELSE IF(PNAM.EQ.'his') THEN
	HIS(KPAR,IPAR)=PVAL
	ELSE IF(PNAM.EQ.'lis') THEN
	LIS(KPAR,IPAR,JPAR)=PVAL
	ELSE IF(PNAM.EQ.'np') THEN
	NP(KPAR)=PVAL
	ELSE IF(PNAM.EQ.'lhp') THEN
	LHP(KPAR)=PVAL
	ELSE IF(PNAM.EQ.'khy') THEN
	KHY(KPAR)=PVAL
	ELSE IF(PNAM.EQ.'kdhy') THEN
	KDHY(KPAR)=PVAL
	ELSE
	PRINT 295
  295   FORMAT (' NO SUCH PNAM')
	STOP
	ENDIF
C
  300   CONTINUE
	RETURN
	END
