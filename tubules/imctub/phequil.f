C  PHEQUIL READS BICARBONATE VALUES AND OUTPUTS EQUILIBRIUM PH
C
	REAL M,PCO2
	DIMENSION M(3000,50)
	CHARACTER*15 IPTFIL,OPTFIL
	INTEGER NCOL(8),JCOL(50)
C
	PRINT 5
    5   FORMAT(/,' PRINT THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
    4   FORMAT (A15)
	OPEN(2,FILE=OPTFIL)
C
	NFIL=1
	NCOLS=0
	KFIL=1
C
	PRINT 2
    2   FORMAT(/,' PRINT THE INPUT FILENAME:',$)
	READ 4, IPTFIL
	OPEN(1,FILE=IPTFIL)
C
	PRINT 10, IPTFIL
   10   FORMAT(' PRINT THE NUMBER OF COLUMNS IN ',A15,$)
	READ *, NCOL(KFIL)
	DO 20 I=1,3000
   20   READ(1,30,END=40) (M(I,J),J=NCOLS+1,NCOLS+NCOL(KFIL))
   30   FORMAT(8E16.8)
   40   NROW=I-1
	CLOSE(1)
   90   NCOLS=NCOLS+NCOL(KFIL)
C
C
  290   PRINT 291
  291   FORMAT(' PRINT THE NUMBER OF HCO3 COLUMNS:',$) 
	READ *, JD
	PRINT 292
  292   FORMAT(' PRINT THE HCO3 COLUMNS: J1,...,JK')
	READ *, (JCOL(I),I=1,JD)
	PRINT 294
  294   FORMAT(' PRINT THE PCO2: ',$)
	READ *, PCO2
	DO 295 K=1,JD
	DO 295 I=1,NROW
  295   M(I,JCOL(K))=6.10 + LOG10(M(I,JCOL(K))/(.03*PCO2))
C
C
	DO 340 I=1,NROW
  340   WRITE(2,350) (M(I,J),J=1,NCOLS)
  350   FORMAT(8E16.8)
C
	STOP
	END
