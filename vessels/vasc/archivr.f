	SUBROUTINE ARCHIVR(IARCH,IVR,ISEG)
C
C PROGRAM TO ARCHIVE OR RETRIEVE A SEGMENTAL SOLUTION FOR EACH VESSEL
C ARCHIVAL VARIABLES IN ANTICIPATION OF MULTIPLE VESSELS
C THERE ARE TWO MEASURES OF DISTANCE ALONG THE VESSEL:
C    X IS AN INTEGER VARIABLE, AND WILL NOW BE ZEROED WITH EACH SEGMENT
C    DIST IS A REAL VARIABLE THAT TRACKS DISTANCE ALONG THE NEPHRON
C
	INTEGER SOLS,X,T,CHOP
C
C PARAMETERS
	DOUBLE PRECISION
     1   Z(15),ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15),
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
C 
C VARIABLES
	DOUBLE PRECISION
     1   VC(901,2),PC(901,2),CC(31,901,2),IMPC(901,2),
     1   LCHC(901),XC(15,901),
     1   VS(901,2),PS(901,2),CS(16,901,2),
     1   IMPS(901,2),LCHS(901),XS(15,901),
     1   CCS(15),SAT(901,2),CCN(901,2),CCX(901,2),CCY(901,2),
     1   SC(901,2),AC(901,2),FVC(901,2),FKC(31,901,2),
     1   QV,QC(31),FVCS(901,2),FKCS(16,901,2),
     1   FVCSE(901,2),FKCSE(16,901,2),FVCSI(901,2),FKCSI(16,901,2),
     1   HCTC(901,2),FBC(901,2),POX(901,2)
C
C KEY TO INDICES:
C  8 vessels: OM (IVR = 1,2) plus 5 IM (IVR = 3-7), plus 1 MR (IVR = 8)
C
C Segment numbers (ISEG)
C      1- OMDVR, MRDVR, OIDVR
C      2- IMDVR
C      3- IMAVR
C      4- OMAVR, MRAVR, OIAVR
C
	INTEGER IARCH,IVR,ISEG
C  IARCH = 0 (ARCHIVE) 
C        = 1 (RETRIEVE EVERYTHING) 
C        = 2 (INITIALIZE THE NEXT SEGMENT)
C        = 3 (RETRIEVE EVERYTHING EXCEPT X=1 VARIABLES)
C
        DOUBLE PRECISION RADVR(8),RDNUM(8),RANUM(8),FSCAL
C  RADVR- RATIO OF AVR NUMBER TO DVR NUMBER FOR EACH VESSEL
C  FSCAL- SCALAR EQUAL TO RADVR FOR THAT VESSEL:
C    IN "ARCHIVE" THE ASCENDING FLOWS CORRESPOND TO THE SUM OF AVR FLOWS
C    WHICH CORRESPOND TO THE SINGLE DVR THAT GAVE BIRTH TO THEM
C
	INTEGER 
     1   RVX(8,4), RVCHOP(8,4)
C PARAMETERS
	DOUBLE PRECISION
     1   RVEPSI(8,4),RVDIST(8,4),RVKHY(8,4),RVKDHY(8,4),
     1   RVDX(8,4),RVETA(8,4),RCL(8,4),RRC0(8,4),RMUC(8,4),
     1   RLPCS(8,4),RSCS(8,4,15),RLCS(8,4,15),RHCS(8,4,15),
     1   RLPCSE(8,4),RSCSE(8,4,15),RLCSE(8,4,15),RHCSE(8,4,15),
     1   RLPCSI(8,4),RSCSI(8,4,15),RLCSI(8,4,15),RHCSI(8,4,15)
C VARIABLES
	DOUBLE PRECISION
     1   RVVS(8,4,901,2),RVPS(8,4,901,2),RVCS(8,4,16,901,2),
     1   RVIMPS(8,4,901,2),RVLCHS(8,4,901),RVXS(8,4,15,901),
     1   RVC(8,4,901,2),RPC(8,4,901,2),RCC(8,4,31,901,2),
     1   RIMPC(8,4,901,2),RLCHC(8,4,901),RXC(8,4,15,901),
     1   RSAT(8,4,901,2),RCCN(8,4,901,2),
     1   RCCX(8,4,901,2),RCCY(8,4,901,2),RPOX(8,4,901,2),
     1   RSC(8,4,901,2),RAC(8,4,901,2),RFVC(8,4,901,2),
     1   RFKC(8,4,31,901,2),RFVCS(8,4,901,2),RFKCS(8,4,16,901,2),
     1   RFVCSE(8,4,901,2),RFKCSE(8,4,16,901,2),RFVCSI(8,4,901,2),
     1   RFKCSI(8,4,16,901,2),RHCTC(8,4,901,2),RFBC(8,4,901,2)
C
	COMMON/PAR/ SOLS,T,X,CHOP,
     1   Z,ZIMPC,RT,RTE,F,
     1   EPSI,SCALE,DT,RTAU,TIME,DIST,
     1   PKC,PKP,PKN,PKF,KHY,KDHY,
     1   LPCS,SCS,LCS,HCS,
     1   LPCSE,SCSE,LCSE,HCSE,
     1   LPCSI,SCSI,LCSI,HCSI,
     1   CL,DX,RC0,MUC,ETA,
     1   KXO,KXD,KYO,KYD,
     1   KCO,KCD,KZO,KZD,
     1   A0,A1,A2,P50,HILL,BOHR
	COMMON/VAR/
     1   VC,PC,CC,VS,PS,CS,IMPC,IMPS,
     1   LCHC,LCHS,XC,XS,CCS,
     1   SAT,CCN,CCX,CCY,HCTC,FBC,POX,
     1   SC,AC,FVC,FKC,QC,QV,FVCS,FKCS,
     1   FVCSE,FKCSE,FVCSI,FKCSI
C
	COMMON /ARCHVR/
     1   RVX,RVCHOP,RADVR,RDNUM,RANUM,
     1   RVEPSI,RVDIST,RVKHY,RVKDHY,
     1   RVDX,RVETA,RCL,RRC0,RMUC,
     1   RLPCS,RSCS,RLCS,RHCS,
     1   RLPCSE,RSCSE,RLCSE,RHCSE,
     1   RLPCSI,RSCSI,RLCSI,RHCSI,
     1   RVVS,RVPS,RVCS,RVIMPS,RVLCHS,RVXS,
     1   RVC,RPC,RCC,RIMPC,RLCHC,RXC,
     1   RSAT,RCCN,RCCX,RCCY,RSC,RAC,
     1   RFVC,RFKC,RFVCS,RFKCS,
     1   RFVCSE,RFKCSE,RFVCSI,
     1   RFKCSI,RHCTC,RFBC,RPOX
C	
C
	GOTO (100,200,300,400) IARCH+1
C
  100   CONTINUE
C
C GENERAL PARAMETERS
C
        RVX(IVR,ISEG) = X
        RVCHOP(IVR,ISEG) = CHOP
        RVEPSI(IVR,ISEG) = EPSI
        RVKHY(IVR,ISEG) = KHY
        RVKDHY(IVR,ISEG) = KDHY
        RVDX(IVR,ISEG) = DX
        RVETA(IVR,ISEG) = ETA
        RCL(IVR,ISEG) = CL
        RRC0(IVR,ISEG) = RC0
        RMUC(IVR,ISEG) = MUC
C
C VESSEL PERMEABILITIES
C
        RLPCS(IVR,ISEG) = LPCS
        RLPCSE(IVR,ISEG) = LPCSE
        RLPCSI(IVR,ISEG) = LPCSI
        DO 110 ISOL=1,SOLS
        RSCS(IVR,ISEG,ISOL) = SCS(ISOL)
        RLCS(IVR,ISEG,ISOL) = LCS(ISOL)
        RHCS(IVR,ISEG,ISOL) = HCS(ISOL)
        RSCSE(IVR,ISEG,ISOL) = SCSE(ISOL)
        RLCSE(IVR,ISEG,ISOL) = LCSE(ISOL)
        RHCSE(IVR,ISEG,ISOL) = HCSE(ISOL)
        RSCSI(IVR,ISEG,ISOL) = SCSI(ISOL)
        RLCSI(IVR,ISEG,ISOL) = LCSI(ISOL)
        RHCSI(IVR,ISEG,ISOL) = HCSI(ISOL)
  110   CONTINUE
C
C INTERSTITIAL AND VESSEL INTRINSIC VARIABLES
C
	DO 120 IX = 1,901
        RVVS(IVR,ISEG,IX,T) = VS(IX,T)
        RVPS(IVR,ISEG,IX,T) = PS(IX,T)
        RVIMPS(IVR,ISEG,IX,T) = IMPS(IX,T)
        RVLCHS(IVR,ISEG,IX) = LCHS(IX)
        RSAT(IVR,ISEG,IX,T) = SAT(IX,T)
        RPOX(IVR,ISEG,IX,T) = POX(IX,T)
	DO 119 ISOL = 1,SOLS
        RVCS(IVR,ISEG,ISOL,IX,T) = CS(ISOL,IX,T)
        RVXS(IVR,ISEG,ISOL,IX) = XS(ISOL,IX)
  119   CONTINUE
  120   CONTINUE
C
	DO 130 IX = 1,CHOP+1
        RVC(IVR,ISEG,IX,T) = VC(IX,T)
        RPC(IVR,ISEG,IX,T) = PC(IX,T)
        RIMPC(IVR,ISEG,IX,T) = IMPC(IX,T)
        RLCHC(IVR,ISEG,IX) = LCHC(IX)
        RCCN(IVR,ISEG,IX,T) = CCN(IX,T)
        RCCX(IVR,ISEG,IX,T) = CCX(IX,T)
        RCCY(IVR,ISEG,IX,T) = CCY(IX,T)
	DO 128 ISOL = 1,SOLS
        RCC(IVR,ISEG,ISOL,IX,T) = CC(ISOL,IX,T)
        RXC(IVR,ISEG,ISOL,IX) = XC(ISOL,IX)
  128   CONTINUE
	DO 129 ISOL = 1+SOLS,31
        RCC(IVR,ISEG,ISOL,IX,T) = CC(ISOL,IX,T)
  129   CONTINUE
  130   CONTINUE
C
C VESSEL FLUXES AND FLOWS
C
	DO 140 IX = 1,CHOP+1
        RSC(IVR,ISEG,IX,T) = SC(IX,T)
        RAC(IVR,ISEG,IX,T) = AC(IX,T)
        RFVCS(IVR,ISEG,IX,T) = FVCS(IX,T)
        RFVC(IVR,ISEG,IX,T) = FVC(IX,T)
        RFVCSE(IVR,ISEG,IX,T) = FVCSE(IX,T)
        RFVCSI(IVR,ISEG,IX,T) = FVCSI(IX,T)
        RHCTC(IVR,ISEG,IX,T) = HCTC(IX,T)
        RFBC(IVR,ISEG,IX,T) = FBC(IX,T)
	DO 138 ISOL = 1,SOLS
        RFKC(IVR,ISEG,ISOL,IX,T) = FKC(ISOL,IX,T)
        RFKCS(IVR,ISEG,ISOL,IX,T) = FKCS(ISOL,IX,T)
        RFKCSE(IVR,ISEG,ISOL,IX,T) = FKCSE(ISOL,IX,T)
        RFKCSI(IVR,ISEG,ISOL,IX,T) = FKCSI(ISOL,IX,T)
  138   CONTINUE
	DO 139 ISOL = 1+SOLS,31
        RFKC(IVR,ISEG,ISOL,IX,T) = FKC(ISOL,IX,T)
  139   CONTINUE
  140   CONTINUE
C
        RETURN
C
  200   CONTINUE
C
C GENERAL PARAMETERS
C
        X = RVX(IVR,ISEG)
        CHOP = RVCHOP(IVR,ISEG)
        EPSI = RVEPSI(IVR,ISEG)
        KHY = RVKHY(IVR,ISEG)
        KDHY = RVKDHY(IVR,ISEG)
        DX = RVDX(IVR,ISEG)
        ETA = RVETA(IVR,ISEG)
        CL = RCL(IVR,ISEG)
        RC0 = RRC0(IVR,ISEG)
        MUC = RMUC(IVR,ISEG)
C
C VESSEL PERMEABILITIES
C
        LPCS = RLPCS(IVR,ISEG)
        LPCSE = RLPCSE(IVR,ISEG)
        LPCSI = RLPCSI(IVR,ISEG)
        DO 210 ISOL=1,SOLS
        SCS(ISOL) = RSCS(IVR,ISEG,ISOL)
        LCS(ISOL) = RLCS(IVR,ISEG,ISOL)
        HCS(ISOL) = RHCS(IVR,ISEG,ISOL)
        SCSE(ISOL) = RSCSE(IVR,ISEG,ISOL)
        LCSE(ISOL) = RLCSE(IVR,ISEG,ISOL)
        HCSE(ISOL) = RHCSE(IVR,ISEG,ISOL)
        SCSI(ISOL) = RSCSI(IVR,ISEG,ISOL)
        LCSI(ISOL) = RLCSI(IVR,ISEG,ISOL)
        HCSI(ISOL) = RHCSI(IVR,ISEG,ISOL)
  210   CONTINUE
C
C INTERSTITIAL AND VESSEL INTRINSIC VARIABLES
C
	DO 220 IX = 1,901
        VS(IX,T) = RVVS(IVR,ISEG,IX,T)
        PS(IX,T) = RVPS(IVR,ISEG,IX,T)
        IMPS(IX,T) = RVIMPS(IVR,ISEG,IX,T)
        LCHS(IX) = RVLCHS(IVR,ISEG,IX)
        SAT(IX,T) = RSAT(IVR,ISEG,IX,T)
        POX(IX,T) = RPOX(IVR,ISEG,IX,T)
	DO 219 ISOL = 1,SOLS
        CS(ISOL,IX,T) = RVCS(IVR,ISEG,ISOL,IX,T)
        XS(ISOL,IX) = RVXS(IVR,ISEG,ISOL,IX)
  219   CONTINUE
  220   CONTINUE
C
	DO 230 IX = 1,CHOP+1
        VC(IX,T) = RVC(IVR,ISEG,IX,T)
        PC(IX,T) = RPC(IVR,ISEG,IX,T)
        IMPC(IX,T) = RIMPC(IVR,ISEG,IX,T)
        LCHC(IX) = RLCHC(IVR,ISEG,IX)
        CCN(IX,T) = RCCN(IVR,ISEG,IX,T)
        CCX(IX,T) = RCCX(IVR,ISEG,IX,T)
        CCY(IX,T) = RCCY(IVR,ISEG,IX,T)
	DO 228 ISOL = 1,SOLS
        CC(ISOL,IX,T) = RCC(IVR,ISEG,ISOL,IX,T)
        XC(ISOL,IX) = RXC(IVR,ISEG,ISOL,IX)
  228   CONTINUE
	DO 229 ISOL = 1+SOLS,31
        CC(ISOL,IX,T) = RCC(IVR,ISEG,ISOL,IX,T)
  229   CONTINUE
  230   CONTINUE
C
C VESSEL FLUXES AND FLOWS
C
	DO 240 IX = 1,CHOP+1
        SC(IX,T) = RSC(IVR,ISEG,IX,T)
        AC(IX,T) = RAC(IVR,ISEG,IX,T)
        FVC(IX,T) = RFVC(IVR,ISEG,IX,T)
        FVCS(IX,T) = RFVCS(IVR,ISEG,IX,T)
        FVCSE(IX,T) = RFVCSE(IVR,ISEG,IX,T)
        FVCSI(IX,T) = RFVCSI(IVR,ISEG,IX,T)
        HCTC(IX,T) = RHCTC(IVR,ISEG,IX,T)
        FBC(IX,T) = RFBC(IVR,ISEG,IX,T)
	DO 238 ISOL = 1,SOLS
        FKC(ISOL,IX,T) = RFKC(IVR,ISEG,ISOL,IX,T)
        FKCS(ISOL,IX,T) = RFKCS(IVR,ISEG,ISOL,IX,T)
        FKCSE(ISOL,IX,T) = RFKCSE(IVR,ISEG,ISOL,IX,T)
        FKCSI(ISOL,IX,T) = RFKCSI(IVR,ISEG,ISOL,IX,T)
  238   CONTINUE
	DO 239 ISOL = 1+SOLS,31
        FKC(ISOL,IX,T) = RFKC(IVR,ISEG,ISOL,IX,T)
  239   CONTINUE
  240   CONTINUE
C
	RETURN
C
C CONTROL PASSES HERE TO INITIALIZE THE NEXT SEGMENT FOR INTEGRATION
C
  300   CONTINUE
	IX = RVCHOP(IVR,ISEG)+1
C
C IDENTIFY THE TRANSITION POINT AT WHICH AVR BEGIN
        IF((IVR.EQ.1).AND.(ISEG.EQ.1)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.2).AND.(ISEG.EQ.1)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.3).AND.(ISEG.EQ.2)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.4).AND.(ISEG.EQ.2)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.5).AND.(ISEG.EQ.2)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.6).AND.(ISEG.EQ.2)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.7).AND.(ISEG.EQ.2)) THEN
        FSCAL=RADVR(IVR)
        ELSE IF((IVR.EQ.8).AND.(ISEG.EQ.1)) THEN
        FSCAL=RADVR(IVR)
        ELSE
        FSCAL=1.0
        ENDIF
C
C LUMINAL VARIABLES
	VC(1,T) = RVC(IVR,ISEG,IX,T)
	PC(1,T) = RPC(IVR,ISEG,IX,T)
	IMPC(1,T) =  RIMPC(IVR,ISEG,IX,T)
	LCHC(1) = RLCHC(IVR,ISEG,IX)
	FVC(1,T) = RFVC(IVR,ISEG,IX,T)/FSCAL
	DO 319 ISOL = 1,SOLS
	CC(ISOL,1,T) = RCC(IVR,ISEG,ISOL,IX,T)
	FKC(ISOL,1,T) = RFKC(IVR,ISEG,ISOL,IX,T)/FSCAL
  319   CONTINUE
	DO 320 ISOL = 1+SOLS,31
        CC(ISOL,1,T) = RCC(IVR,ISEG,ISOL,IX,T)
        FKC(ISOL,1,T) = RFKC(IVR,ISEG,ISOL,IX,T)/FSCAL
  320   CONTINUE
        HCTC(1,T) = RHCTC(IVR,ISEG,IX,T)
        RETURN
C
C RETRIEVE EVERYTHING EXCEPT X=1 VARIABLES
C
  400   CONTINUE
C
C GENERAL PARAMETERS
C
        X = RVX(IVR,ISEG)
        CHOP = RVCHOP(IVR,ISEG)
        EPSI = RVEPSI(IVR,ISEG)
        KHY = RVKHY(IVR,ISEG)
        KDHY = RVKDHY(IVR,ISEG)
        DX = RVDX(IVR,ISEG)
        ETA = RVETA(IVR,ISEG)
        CL = RCL(IVR,ISEG)
        RC0 = RRC0(IVR,ISEG)
        MUC = RMUC(IVR,ISEG)
C
C VESSEL PERMEABILITIES
C
        LPCS = RLPCS(IVR,ISEG)
        LPCSE = RLPCSE(IVR,ISEG)
        LPCSI = RLPCSI(IVR,ISEG)
        DO 410 ISOL=1,SOLS
        SCS(ISOL) = RSCS(IVR,ISEG,ISOL)
        LCS(ISOL) = RLCS(IVR,ISEG,ISOL)
        HCS(ISOL) = RHCS(IVR,ISEG,ISOL)
        SCSE(ISOL) = RSCSE(IVR,ISEG,ISOL)
        LCSE(ISOL) = RLCSE(IVR,ISEG,ISOL)
        HCSE(ISOL) = RHCSE(IVR,ISEG,ISOL)
        SCSI(ISOL) = RSCSI(IVR,ISEG,ISOL)
        LCSI(ISOL) = RLCSI(IVR,ISEG,ISOL)
        HCSI(ISOL) = RHCSI(IVR,ISEG,ISOL)
  410   CONTINUE
C
C INTERSTITIAL AND VESSEL INTRINSIC VARIABLES
C
	DO 420 IX = 1,901
        VS(IX,T) = RVVS(IVR,ISEG,IX,T)
        PS(IX,T) = RVPS(IVR,ISEG,IX,T)
        IMPS(IX,T) = RVIMPS(IVR,ISEG,IX,T)
        LCHS(IX) = RVLCHS(IVR,ISEG,IX)
        SAT(IX,T) = RSAT(IVR,ISEG,IX,T)
        POX(IX,T) = RPOX(IVR,ISEG,IX,T)
	DO 419 ISOL = 1,SOLS
        CS(ISOL,IX,T) = RVCS(IVR,ISEG,ISOL,IX,T)
        XS(ISOL,IX) = RVXS(IVR,ISEG,ISOL,IX)
  419   CONTINUE
  420   CONTINUE
        CCN(1,T) = RCCN(IVR,ISEG,1,T)
        CCX(1,T) = RCCX(IVR,ISEG,1,T)
        CCY(1,T) = RCCY(IVR,ISEG,1,T)
C
C
	RETURN
        END