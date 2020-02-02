	PROGRAM PCR
C TO INVOKE ALL OF THE SUBDIRECTORY PCR PROGRAMS
C
	INTEGER NSETS
C
	PRINT 60
   60   FORMAT(' NUMBER OF PARAMETER SETS NEEDED=',$)
	READ *, NSETS
C
        CALL OMDVRPCR(NSETS)
        CALL OMAVRPCR(NSETS)
        CALL IMDVRPCR(NSETS)
        CALL IMAVRPCR(NSETS)
        CALL OIDVRPCR(NSETS)
        CALL OIAVRPCR(NSETS)
        CALL MRDVRPCR(NSETS)
        CALL MRAVRPCR(NSETS)
C
	STOP
	END
C
        SUBROUTINE OMDVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='omdvr/omdparam.tem')
	OPEN (UNIT=20,FILE='omdvr/omdparam.dat')
	OPEN (UNIT=21,FILE='omdvr/omdparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE OMAVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='omavr/omaparam.tem')
	OPEN (UNIT=20,FILE='omavr/omaparam.dat')
	OPEN (UNIT=21,FILE='omavr/omaparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE IMDVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='imdvr/imdparam.tem')
	OPEN (UNIT=20,FILE='imdvr/imdparam.dat')
	OPEN (UNIT=21,FILE='imdvr/imdparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE IMAVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='imavr/imaparam.tem')
	OPEN (UNIT=20,FILE='imavr/imaparam.dat')
	OPEN (UNIT=21,FILE='imavr/imaparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE OIDVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='oidvr/oidparam.tem')
	OPEN (UNIT=20,FILE='oidvr/oidparam.dat')
	OPEN (UNIT=21,FILE='oidvr/oidparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE OIAVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='oiavr/oiaparam.tem')
	OPEN (UNIT=20,FILE='oiavr/oiaparam.dat')
	OPEN (UNIT=21,FILE='oiavr/oiaparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE MRDVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='mrdvr/mrdparam.tem')
	OPEN (UNIT=20,FILE='mrdvr/mrdparam.dat')
	OPEN (UNIT=21,FILE='mrdvr/mrdparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
C
        SUBROUTINE MRAVRPCR(NSETS)
C PROGRAM TO CREATE THE CPARAM.DAT FILE FOR EXECUTION OF CAP
C   PARAMETERS ARE READ FROM THE TEMPLATE FILE CPARAM.TEM 
C
	INTEGER SOLS,TAU,TLIM,CHOP
C
	DOUBLE PRECISION EPSI,RT,
     1   CL,DX,RC0,MUC,ETA,ZIMPC,
     1   LPCS,SCS(15),LCS(15),HCS(15),
     1   LPCSE,SCSE(15),LCSE(15),HCSE(15),
     1   LPCSI,SCSI(15),LCSI(15),HCSI(15)
C
C
	OPEN (UNIT=19,FILE='mravr/mraparam.tem')
	OPEN (UNIT=20,FILE='mravr/mraparam.dat')
	OPEN (UNIT=21,FILE='mravr/mraparam.cap')
C
	SOLS=15
	RT=0.193D+05
C
C
C GEOMETRIC AND MEMBRANE PROPERTIES ARE READ. 
C
	READ (19,30) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
   30   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	READ (19,26) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
   26   FORMAT (D12.4,/,(F6.3,D12.4))
C
C
	DO 200 KSETS=1,NSETS
C
	WRITE (20,130) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,ZIMPC,LPCSE,
     1  (SCSE(I),HCSE(I),I=1,SOLS)
  130   FORMAT (3I5,D9.1/,5D12.4,/,D12.4,/,(F6.3,D12.4))
C
	WRITE (20,126) LPCSI,
     1  (SCSI(I),HCSI(I),I=1,SOLS)
  126   FORMAT (D12.4,/,(F6.3,D12.4))
C
C REWRITE THE PARAMETERS WITH PARALLEL ELEMENTS COMBINED
C  FOR USE IN A CAPILLARY MODEL
C
        LPCS = (LPCSE + LPCSI)/RT
        DO 140 I=1,SOLS
        SCS(I) = (SCSE(I)*LPCSE + SCSI(I)*LPCSI)/(LPCSE+LPCSI)
        HCS(I) = HCSE(I) + HCSI(I)
  140   CONTINUE
C
	WRITE (21,150) TAU,TLIM,CHOP,EPSI,
     1  CL,RC0,MUC,ETA,
     1  ZIMPC,LPCS,
     1  (SCS(I),HCS(I),I=1,SOLS)
  150   FORMAT (3I5,D9.1/,4D12.4,/,2D12.4,/,(F6.3,D12.4))
C
  200   CONTINUE
	STOP
	END
