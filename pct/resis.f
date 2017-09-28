	PROGRAM RESIS
C  PROGRAM TO COMPUTE EPITHELIAL RESISTANCES FROM A VOLTAGE STEP EXPT.
C
C  THE EXPERIMENT STARTS WITH A STEADY STATE SOLUTION
C  AND A VOLTAGE INCREMENT IS FOLLOWED FOR A FEW TIME STEPS
C  THE CHANGE IN INTRAEPITHELIAL CURRENTS AND VOLTAGES ARE
C  READ AND THE LAST TIME STEP OF THE RECORD IS ANALYZED
C
	CHARACTER*13 IPTFIL(8),OPTFIL,IPTNAM
	REAL RESUL(105,9),
     1  VM(8,9),VI(8,9),VE(8,9),VS(8,9),
     1  IMI(8,9),IME(8,9),IIE(8,9),IIS(8,9),IES(8,9),
     1  DVM(8),DVI(8),DVE(8),DVS(8),
     1  DVMI(8),DVME(8),DVIE(8),DVIS(8),DVES(8),DVMIIS(8),
     1  DIMI(8),DIME(8),DIIE(8),DIIS(8),DIES(8),
     1  RMI(8),RME(8),RIE(8),RIS(8),RES(8),RTE(8)
C
    4   FORMAT(A13)
	PRINT 5
    5   FORMAT(' TYPE THE OUTPUT FILENAME:',$)
	READ 4, OPTFIL
C
	OPEN(2,FILE=OPTFIL)
C
	PRINT 1
    1   FORMAT(' TYPE THE NUMBER OF INPUT FILES ',$)
	READ *, NFIL
	DO 90 K=1,NFIL
C
	PRINT 2
    2   FORMAT(/,' TYPE THE INPUT FILENAME:',$)
	READ 4, IPTFIL(K)
    6   IPTNAM=IPTFIL(K)
	OPEN(1,FILE=IPTNAM)
C
C
  300   READ(1,310,END=900) (RESUL(1,J),J=1,5)
  310   FORMAT (73(/),10X,F9.4,////,6X,4D12.4)
	READ(1,315) (RESUL(2,J),J=1,4)
  315   FORMAT(///,8X,4F12.3)
	READ(1,320) (RESUL(3,J),J=1,4)
  320   FORMAT(///,8X,4F12.3)
	READ(1,325) ((RESUL(I,J),J=1,4),I=4,7)
  325   FORMAT(///,(9X,4F12.6))
	READ(1,330) (RESUL(8,J),J=1,4)
  330   FORMAT(10X,4D12.4)
	READ(1,335) ((RESUL(I,J),J=1,4),I=9,12)
  335   FORMAT(9X,4F12.6)
	READ(1,330) (RESUL(13,J),J=1,4)
	READ(1,335) (RESUL(14,J),J=1,4)
	READ(1,330) (RESUL(15,J),J=1,4)
	READ(1,335) (RESUL(16,J),J=1,4)
	READ(1,330) (RESUL(17,J),J=1,4)
	READ(1,335) (RESUL(18,J),J=1,4)
	READ(1,340) RESUL(19,1),(RESUL(19,J),J=3,4)
  340   FORMAT(9X,F12.6,12X,2F12.6)
	READ(1,335) ((RESUL(I,J),J=1,4),I=20,21)
	READ(1,345) ((RESUL(I,J),J=3,3),I=22,23)
  345   FORMAT(33X,1F12.6)
	READ(1,350) ((RESUL(I,J),J=1,3),I=24,38)
  350   FORMAT(///,(9X,3F12.4))
	READ(1,355) (RESUL(39,J),J=1,3)
  355   FORMAT(///,9X,3D12.4)
	READ(1,360) ((RESUL(I,J),J=1,3),I=40,54)
  360   FORMAT(///,(9X,3D12.4))
	READ(1,365) (RESUL(55,J),J=1,3)
  365   FORMAT(9X,3D12.4)
	READ(1,370) RESUL(56,1),RESUL(56,2),(RESUL(56,J),J=4,5)
  370   FORMAT(///,9X,2D12.4,12X,2D12.4)
	READ(1,375) (RESUL(I,1),RESUL(I,2),(RESUL(I,J),J=4,5),I=57,71)
  375   FORMAT(///,(9X,2D12.4,12X,2D12.4))
	READ(1,380) (RESUL(72,J),J=1,2)
  380   FORMAT(9X,2D12.4)
	READ(1,385) ((RESUL(I,J),J=1,3),I=73,87)
  385   FORMAT(///,(9X,3D12.4))
	READ(1,390) ((RESUL(I,J),J=1,3),I=88,102)
  390   FORMAT(///,(9X,3D12.4))
C
C
	ITIM=1
	VM(K,ITIM)=RESUL(2,1)
	VI(K,ITIM)=RESUL(2,3)
	VE(K,ITIM)=RESUL(2,2)
	VS(K,ITIM)=RESUL(2,4)
	IMI(K,ITIM)=RESUL(55,1)
	IME(K,ITIM)=RESUL(72,1)
	IES(K,ITIM)=RESUL(72,2)
	IIS(K,ITIM)=RESUL(55,3)
	IIE(K,ITIM)=IMI(K,ITIM)-IIS(K,ITIM)
C
C
  500   READ(1,510,END=900) (RESUL(1,J),J=1,5)
  510   FORMAT (/,10X,F9.4,////,6X,4D12.4)
	READ(1,515) (RESUL(2,J),J=1,4)
  515   FORMAT(///,8X,4F12.3)
	READ(1,520) (RESUL(3,J),J=1,4)
  520   FORMAT(///,8X,4F12.3)
	READ(1,525) ((RESUL(I,J),J=1,4),I=4,7)
  525   FORMAT(///,(9X,4F12.6))
	READ(1,530) (RESUL(8,J),J=1,4)
  530   FORMAT(10X,4D12.4)
	READ(1,535) ((RESUL(I,J),J=1,4),I=9,12)
  535   FORMAT(9X,4F12.6)
	READ(1,530) (RESUL(13,J),J=1,4)
	READ(1,535) (RESUL(14,J),J=1,4)
	READ(1,530) (RESUL(15,J),J=1,4)
	READ(1,535) (RESUL(16,J),J=1,4)
	READ(1,530) (RESUL(17,J),J=1,4)
	READ(1,535) (RESUL(18,J),J=1,4)
	READ(1,540) RESUL(19,1),(RESUL(19,J),J=3,4)
  540   FORMAT(9X,F12.6,12X,2F12.6)
	READ(1,535) ((RESUL(I,J),J=1,4),I=20,21)
	READ(1,545) ((RESUL(I,J),J=3,3),I=22,23)
  545   FORMAT(33X,1F12.6)
	READ(1,550) ((RESUL(I,J),J=1,3),I=24,38)
  550   FORMAT(///,(9X,3F12.4))
	READ(1,555) (RESUL(39,J),J=1,3)
  555   FORMAT(///,9X,3D12.4)
	READ(1,560) ((RESUL(I,J),J=1,3),I=40,54)
  560   FORMAT(///,(9X,3D12.4))
	READ(1,565) (RESUL(55,J),J=1,3)
  565   FORMAT(9X,3D12.4)
	READ(1,570) RESUL(56,1),RESUL(56,2),(RESUL(56,J),J=4,5)
  570   FORMAT(///,9X,2D12.4,12X,2D12.4)
	READ(1,575) (RESUL(I,1),RESUL(I,2),(RESUL(I,J),J=4,5),I=57,71)
  575   FORMAT(///,(9X,2D12.4,12X,2D12.4))
	READ(1,580) (RESUL(72,J),J=1,2)
  580   FORMAT(9X,2D12.4)
	READ(1,585) ((RESUL(I,J),J=1,3),I=73,87)
  585   FORMAT(///,(9X,3D12.4))
	READ(1,590) ((RESUL(I,J),J=1,3),I=88,102)
  590   FORMAT(///,(9X,3D12.4))
C
C
	ITIM=ITIM+1
	VM(K,ITIM)=RESUL(2,1)
	VI(K,ITIM)=RESUL(2,3)
	VE(K,ITIM)=RESUL(2,2)
	VS(K,ITIM)=RESUL(2,4)
	IMI(K,ITIM)=RESUL(55,1)
	IME(K,ITIM)=RESUL(72,1)
	IES(K,ITIM)=RESUL(72,2)
	IIS(K,ITIM)=RESUL(55,3)
	IIE(K,ITIM)=IMI(K,ITIM)-IIS(K,ITIM)
	GO TO 500
C
C
  900   CONTINUE
C
	CLOSE(UNIT=1)
C
	DVM(K)=VM(K,ITIM)-VM(K,1)
	DVI(K)=VI(K,ITIM)-VI(K,1)
	DVE(K)=VE(K,ITIM)-VE(K,1)
	DVS(K)=VS(K,ITIM)-VS(K,1)
	DIMI(K)=IMI(K,ITIM)-IMI(K,1)
	DIME(K)=IME(K,ITIM)-IME(K,1)
	DIES(K)=IES(K,ITIM)-IES(K,1)
	DIIS(K)=IIS(K,ITIM)-IIS(K,1)
	DIIE(K)=IIE(K,ITIM)-IIE(K,1)
C
C
	DVMI(K)=DVM(K)-DVI(K)
	DVME(K)=DVM(K)-DVE(K)
	DVIE(K)=DVI(K)-DVE(K)
	DVIS(K)=DVI(K)-DVS(K)
	DVES(K)=DVE(K)-DVS(K)
	DVMIIS(K)=DVMI(K)/DVIS(K)
	RMI(K)=DVMI(K)/DIMI(K)
	RME(K)=DVME(K)/DIME(K)
	RIE(K)=DVIE(K)/DIIE(K)
	RIS(K)=DVIS(K)/DIIS(K)
	RES(K)=DVES(K)/DIES(K)
	RTE(K)=(DVM(K)-DVS(K))/(DIMI(K)+DIME(K))
C
   90   CONTINUE
C
	WRITE(2,35) (IPTFIL(K),
     1  VM(K,1),VI(K,1),VE(K,1),VS(K,1),
     1  IMI(K,1),IME(K,1),IIE(K,1),IIS(K,1),IES(K,1),K=1,NFIL)
	WRITE(2,45) (IPTFIL(K),
     1  DVM(K),DVI(K),DVE(K),DVS(K),
     1  DVMI(K),DVME(K),DVIE(K),DVIS(K),DVES(K),DVMIIS(K),K=1,NFIL)
	WRITE(2,55) (IPTFIL(K),
     1  DIMI(K),DIME(K),DIIE(K),DIIS(K),DIES(K),
     1  RMI(K),RME(K),RIE(K),RIS(K),RES(K),RTE(K),K=1,NFIL)
C
   35   FORMAT(6X,'EXPT.',9X,
     1  'VM',8X,'VI',8X,'VE',8X,'VS',8X,10X,
     1  'IMI',7X,'IME',7X,'IIE',7X,'IIS',7X,'IES',//,
     1  7(1X,A13,4F10.3,10X,1P,5E10.2,0P,/))
   45   FORMAT(//,6X,'EXPT.',9X,
     1  'DVM',7X,'DVI',7X,'DVE',7X,'DVS',7X,10X,
     1  'DVMI',6X,'DVME',6X,'DVIE',6X,'DVIS',6X,'DVES',
     1  5X,'DVMIIS',//,
     1  7(1X,A13,4F10.3,10X,6F10.3,/))
   55   FORMAT(//,6X,'EXPT.',9X,
     1  'DIMI',6X,'DIME',6X,'DIIE',6X,'DIIS',6X,'DIES',6X,
     1  'RMI',7X,'RME',7X,'RIE',7X,'RIS',7X,'RES',7X,'RTE',//,
     1  7(1X,A13,1P,11E10.2,0P,/))
C
	STOP
	END
