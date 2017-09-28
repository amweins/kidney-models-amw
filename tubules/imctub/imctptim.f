	SUBROUTINE IMCTPTIM(NPAR)
C PROGRAM TO READ AND INITIALIZE THE PARAMETERS FOR IMCT
C
	INTEGER SOLS,T,X,NPAR
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   Z(12),RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY(5),KDHY(5),L0,L(801,2)
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   VM(801,2),PM(801,2),CM(12,801,2),
     1   IMPM(801,2),LCHM(801),XM(12,801),
     1   VS(801,2),PS(801,2),CS(12,801,2),
     1   IMPS(801,2),LCHS(801),XS(12,801),
     1   TL,DX,RM0,MUM,ETA,
     1   SM(801,2),AM(801,2),FVM(801,2),FKM(13,801,2)
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   AME,AE0,AE(801,2),MUA,CHVL0,CHVL(801,2),MUV,
     1   LPME,LPES,SME(12),SES(12),
     1   HME(12),HES(12),CME(12),CES(12),
     1   VE(801,2),PE(801,2),CE(12,801,2),LCHE(801),XE(12,801),
     1   FEVM(801,2),FEKM(12,801,2),FEVS(801,2),FEKS(12,801,2),CURE(801)
C CELL PARAMETERS
	DOUBLE PRECISION
     1   AIE(3),AI0(3),CLVL0(3),IMP0(3),CLVL(3,801,2),
     1   ZIMP(3),TBUF(3),PKB(3),CBUF(3,801,2),HCBUF(3,801,2),
     1   LPMI(3),LPIS(3),SMI(3,12),SIS(3,12),
     1   HMI(3,12),HIS(3,12),CMI(3,12),CIE(3,12),CIS(3,12),
     1   LMI(3,12,12),LIS(3,12,12),ATMI(3,12,801),ATIS(3,12,801),
     1   NP(3),KNPN(3),KNPK(3),KNH4(3),NPHK(3),
     1   LHP(3),XIHP(3),XHP(3),NAE1(3),
     1   VI(3,801,2),PI(3,801,2),CI(3,12,801,2),IMP(3,801),LCHI(3,801),
     1   XI(3,12,801),FIVM(3,801,2),FIKM(3,12,801,2),FIVS(3,801,2),
     1   FIKS(3,12,801,2),CURI(3,801),JV(3,801,2),JK(3,12,801,2)
C
        COMMON SOLS,T,X,
     1   Z,RT,RTE,F,
     1   EPSI,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,KHY,KDHY,L0,L,
     1   VM,PM,CM,IMPM,LCHM,XM,
     1   VS,PS,CS,IMPS,LCHS,XS,
     1   TL,DX,RM0,MUM,ETA,
     1   SM,AM,FVM,FKM,
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
     1   LHP,XIHP,XHP,NAE1,
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
	CHARACTER*5 PNAM
C
C
	GO TO (100,200) T
C
  100   READ (46,110) NPAR
  110   FORMAT(I3)
	RETURN
C
  200   IF (NPAR.EQ.0) RETURN
	DO 300 M=1,NPAR
	READ (46,210) PNAM,KPAR,IPAR,JPAR,PVAL
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
