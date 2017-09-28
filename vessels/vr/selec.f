        SUBROUTINE SELEC (VNAME1,VNAME2,IWR,JWR,IC)
C
	INTEGER NWR,IWR,JWR,IC
        CHARACTER*8 VNAME1,VNAME2
C
	IF(VNAME1.EQ.'vc') THEN 
        IWR=1+2*IC
        JWR=2
        GO TO 800
        ELSE IF(VNAME1.EQ.'pc') THEN 
        IWR=1+IC
        JWR=2
        GO TO 800
        ELSE IF(VNAME1.EQ.'cc') THEN 
        IWR=1+IC
        JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'phc') THEN 
        IWR=1+2*IC
        JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'xc') THEN 
        IWR=1+2*IC
        JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'fvc') THEN 
        IWR=1+3*IC
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'fkc') THEN 
        IWR=1+3*IC
	JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'fbc') THEN 
        IWR=1+3*IC
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'hct') THEN 
        IWR=1+3*IC
	JWR=14
        GO TO 800
C
	ELSE IF(VNAME1.EQ.'ps') THEN 
        IWR=1
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'cs') THEN 
        IWR=1
	JWR=3
        GO TO 200
	ELSE IF(VNAME1.EQ.'phs') THEN 
        IWR=1
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'fvcs') THEN 
        IWR=2+4*IC
	JWR=2
        GO TO 800
	ELSE IF(VNAME1.EQ.'fkcs') THEN 
        IWR=2+4*IC
	JWR=3
        GO TO 200
C
	ELSE IF(VNAME1.EQ.'apr') THEN 
        IWR=1+4*IC
	JWR=2
        GO TO 200
	ELSE IF(VNAME1.EQ.'cur') THEN 
        IWR=2+4*IC
	JWR=13
        GO TO 800
	ELSE IF(VNAME1.EQ.'h') THEN 
        IWR=2+5*IC
	JWR=3
        GO TO 800
	ELSE IF(VNAME1.EQ.'hgb') THEN 
        GO TO 250
	ELSE IF(VNAME1.EQ.'fbase') THEN 
        IWR=2+6*IC
	JWR=10
        GO TO 800
        ELSE
        PRINT 40
   40   FORMAT(' NO SUCH VARIABLE NAME, TRY AGAIN')
        ENDIF
C
C
  200   CONTINUE
	IF(VNAME2.EQ.'vol') THEN 
        JWR=JWR-1
	GO TO 800
	ELSE IF(VNAME2.EQ.'na') THEN 
        JWR=JWR+0
	GO TO 800
	ELSE IF(VNAME2.EQ.'k') THEN 
        JWR=JWR+1
	GO TO 800
	ELSE IF(VNAME2.EQ.'cl') THEN 
        JWR=JWR+2
	GO TO 800
	ELSE IF(VNAME2.EQ.'hco3') THEN 
        JWR=JWR+3
	GO TO 800
	ELSE IF(VNAME2.EQ.'hpo') THEN 
        JWR=JWR+4
	GO TO 800
	ELSE IF(VNAME2.EQ.'h2po') THEN 
        JWR=JWR+5
	GO TO 800
	ELSE IF(VNAME2.EQ.'urea') THEN 
        JWR=JWR+6
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh4') THEN 
        JWR=JWR+7
	GO TO 800
	ELSE IF(VNAME2.EQ.'hco2') THEN 
        JWR=JWR+8
	GO TO 800
	ELSE IF(VNAME2.EQ.'gluc') THEN 
        JWR=JWR+9
	GO TO 800
	ELSE IF(VNAME2.EQ.'prot') THEN 
        JWR=10
	GO TO 800
	ELSE IF(VNAME2.EQ.'base') THEN 
        JWR=10
	GO TO 800
        ELSE
        PRINT 240
  240   FORMAT(' NO SUCH SOLUTE NAME, TRY AGAIN')
        ENDIF
C
C
C A HGB SPECIES HAS BEEN DESIGNATED
C
  250   CONTINUE
	IF(VNAME2.EQ.'onh2') THEN 
        IWR=2+5*IC
        JWR=4
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh2') THEN 
        IWR=2+5*IC
        JWR=5
	GO TO 800
	ELSE IF(VNAME2.EQ.'onh3') THEN 
        IWR=2+5*IC
        JWR=6
	GO TO 800
	ELSE IF(VNAME2.EQ.'nh3') THEN 
        IWR=2+5*IC
        JWR=7
	GO TO 800
	ELSE IF(VNAME2.EQ.'onhco2') THEN 
        IWR=2+5*IC
        JWR=8
	GO TO 800
	ELSE IF(VNAME2.EQ.'nhco2') THEN 
        IWR=2+5*IC
        JWR=9
	GO TO 800
	ELSE IF(VNAME2.EQ.'oxh') THEN 
        IWR=2+5*IC
        JWR=2
	GO TO 800
	ELSE IF(VNAME2.EQ.'xh') THEN 
        IWR=2+6*IC
        JWR=3
	GO TO 800
	ELSE IF(VNAME2.EQ.'ox') THEN 
        IWR=2+6*IC
        JWR=4
	GO TO 800
	ELSE IF(VNAME2.EQ.'x') THEN 
        IWR=2+6*IC
        JWR=5
	GO TO 800
	ELSE IF(VNAME2.EQ.'oyh') THEN 
        IWR=2+6*IC
        JWR=6
	GO TO 800
	ELSE IF(VNAME2.EQ.'yh') THEN 
        IWR=2+6*IC
        JWR=7
	GO TO 800
	ELSE IF(VNAME2.EQ.'oy') THEN 
        IWR=2+6*IC
        JWR=8
	GO TO 800
	ELSE IF(VNAME2.EQ.'y') THEN 
        IWR=2+6*IC
        JWR=9
	GO TO 800
        ELSE
        PRINT 260
  260   FORMAT(' NO SUCH HGB NAME, TRY AGAIN')
        ENDIF
C
  800   RETURN
        END
