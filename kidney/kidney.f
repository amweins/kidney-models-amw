	PROGRAM KIDNEY
C
C 9/3/15 PROGRAM TO MERGE NEPHRON WITH VASC TO SOLVE FOR MEDULLARY INTERSITITAL PROFILE 
C  Interstitial concentrations and pressures are determined on a coarse grid to permit 
C  a tractable number of variables.  These grid values are refined by linear interpolation.
C At this level of the program, there will only be archived variables, plus the necessary
C  CS variables, which are updated on their grids.
C
C **FOR THIS PROGRAM TO WORK, VARIABLES OMCHOP AND IMCHOP MUST BE WELL-DEFINED
C    OMCHOP = CHOP(PST) + CHOP(sDHL) (=CHOP(lDHLu) = CHOP(AHLm) = CHOP(OMCD)
C    IMCHOP = CHOP(lDHLl) = CHOP(tAHL) = CHOP(IMCD)
C    MRCHOP = CHOP(AHLc) = CHOP(CCD)
C
C Boundary dat, specified and variable, will be set at the top level (kidney) and placed
C   in archive variables.
C                        Lumen (Init)            IS (Init)            IS (End)
C       SFPCT,JMPCT         KBSET            KBSET (Cortex)         KBSET (Cortex)
C       SFPST,JMPST           -              KGAM (OM - init)       KGAM (OM - mid)
C       sDHL,lDHLu            -              KGAM (OM - mid)        KGAM (OM - end)
C       lDHLl,tAHL            -              KGAM (OM - end)        KGAM (IM)
C       AHLm                  -              KGAM (OM - end)        KGAM (MR - end)
C       AHLc                  -              KGAM (MR - end)        KGAM (MR - init)
C       DCT,CNT               -              KBSET (Cortex)         KBSET (Cortex)
C       CCT                   -              KGAM (MR - init)       KGAM (MR - end)
C       OMCT                  -              KGAM (MR - end)        KGAM (OM - end)
C       IMCT                  -              KGAM (OM - end)        KGAM (IM - end)
C
C       OMDVR,OMAVR (IS)    KBSET            KGAM (OM - init)       KGAM (OM - end)
C       OIDVR,OIAVR         KBSET            KGAM (MR - end)        KGAM (OM - end)
C       IMDVR,IMAVR           -              KGAM (OM - end)        KGAM (IM)
C       MRDVR,MRAVR (IS)    KBSET            KGAM (MR - init)       KGAM (MR - end)
C
C 12/29/13  PROGRAM TO MERGE PN WITH DN TO CREATE A FULL NEPHRON
C Novel additions are 
C  Solving the 6 nephrons through dct,
C    then rapid mixing to enter a common CNT (no CO2 hydration/dehydration at the mixing point)
C  Entering CNT hydrostatic pressure will depend upon end-IMCD outlet pressure (nominally 0 mmHg)
C  Since CNT pressure depends upon CNT flow, it depends upon flow delivered by each nephron, 
C    so that nephrons are weakly coupled.
C  Since proximal nephron flux is pressure-dependent, initial PCT flow and pressure are 
C    subject to TGF and matching CNT inlet pressure, in all 13 global variables.
C  
C  Calculation plan: 
C  SFPCT -> SFPST -> SDHL -> AHLm -> AHLc -> DCT
C  JMPCT -> JMPST -> LDHLu -> LDHLl(2-6) -> tAHL(2-6) -> AHLm(2-6) -> AHLc(2-6) -> DCT(2-6)
C  Merged DCT -> CNT -> CCD -> OMCD -> IMCD
C 
C Carbonic Anhydrase
C   Clearly less in PST; there is some in sDHL; largely absent in IM.
C   Model parameters: 1% of full catalysis in PST, sDHL, and lDHLu; absent (0.01%) in lDHLl and tAHL
C   Will assume full catalysis in IS (in view of proximity to RBC).
C
C 12/22/12 - ARCHIVAL VARIABLES ADDED IN ANTICIPATION OF MULTIPLE NEPHRONS
C THERE ARE TWO MEASURES OF DISTANCE ALONG THE NEPHRON:
C    X IS AN INTEGER VARIABLE, AND WILL NOW BE ZEROED WITH EACH SEGMENT
C    DIST IS A REAL VARIABLE THAT TRACKS DISTANCE ALONG THE NEPHRON
C    ALTHOUGH THIS WILL BECOME AMBIGUOUS AS DIFFERENT LENGTHS COALESCE IN CNT
C
C 2/8/13 - TGF ADJUSTS SNGFR ACCORDING TO MACULA DENSA [CL-]i 
C    AND EARLY PT PRESSURE ACCORDING TO DCT PRESSURE.
C
C
C LOCAL VARIABLES
C
	INTEGER SOLS,T,X,COUNT,EXP,
     1    CHOP,OMCHOP,IMCHOP,MRCHOP,SLICE
C
C       SOLS-   NUMBER OF SOLUTES
C       T-      INDEX FOR OLD (=1) OR NEW (=2) TIME STEP
C       COUNT-  NUMBER OF EXPT. BEING SIMULATED
C       EXP-    TOTAL NUMBER OF EXPTS. IN THE RUN
C	X-	INDICATES SPATIAL STEP
C	CHOP-	SPATIAL CHOP OF TUBULE
C       OMCHOP,IMCHOP,MRCHOP - REGIONAL CHOP FOR OM, IM AND MR
C       SLICE-  OM+IM+MR DISCRETIZATION OF GUESSES 
C               (REFINED TO OMCHOP+IMCHOP+MRCHOP BY LINEAR INTERPOLATION)
C               GAMMA VECTOR IS BASED ON SLICE
C
        CHARACTER*8 SOL(16),ALABL(16)
        CHARACTER*12 NIL,VMAT(16,8)
C
C FOR THE BASELINE CONFIRGURATION:
C  SLICE DEMARCATIONS AT 0 - 7 MM FROM OM TO PAPILLA
C Interstitial unknows are pressure and concentration at 3 points in OM, 5 points
C  in IM, and two points in MR.  The x=0 point in OM abuts PST; there is a second
C  x=0 point in OM that is identified with MR=2mm, and this abuts AHL and OMCD.
C  Within MR, the x=0 point abuts AHLc and CCD.
C
C SPECIFIED BOUNDARY DATA
C
        DOUBLE PRECISION
     1   CSMR(15,201,2),CSOM(15,101,2),CSIM(15,501,2),
     1   PSMR(201,2),PSOM(101,2),PSIM(501,2),
     1   CSGAM(15,51),PSGAM(51),KEPSI,
     1   PS0(3),IMPM0(2),IMPS0(3),CM0(15,2),CS0(16,3), 
     1   FVC0(3),PC0(3),IMPC0(3),HCT0(3),CC0(16,3)
C
C       CSMR,CSOM,CSIM- CS CONCENTRATIONS AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C       PSMR,PSOM,PSIM- PS VALUES AT GRID POINTS (FOR EXPORT TO NEPHRON AND VASC)
C       CSGAM-  GUESS CS CONCENTRATIONS AT SLICE BOUNDARIES
C       PSGAM-  GUESS PS VALUES AT SLICE BOUNDARIES
C       KEPSI-   ERROR TOLERANCE FOR THE KIDNEY ITERATIONS
C
C	FVMO(1) and FVM0(2)- SFPCT and JMPCT(2) flows for all SF and JM(2) nephrons
C	PSO- 	SFPCT, JMPCT, and DCT/CNT vascular pressures (uniform)
C	IMPM0- filtered impermeant osmolality
C	IMPS0- vascular impermeant osmolality for SFPCT, JMPCT, and DCT/CNT (mM)
C	CM0-	glomerular filtrate composition
C	CS0-	vascular plasma composition for SFPCT, JMPCT, and DCT/CNT
C	CSMR-	cortical composition within medullary ray
C
C	FVC0-	capillary plasma flow for individual OMDVR, OIDVR, and MRDVR
C	PC0-	capillary pressures
C	IMPC0- 	capillary protein concentration (gm/dl)
C	HCT0-	hematocrits
C	CC0-	capillary solute composition
C
C VARIABLES UTILIZED IN THE GEN AND VGEN SUBROUTINES
C
        DOUBLE PRECISION
     1    OMDIST(901),IMDIST(901),MRDIST(901),
     1    OMJVMS(901),OMJKMS(16,901),
     1    IMJVMS(901),IMJKMS(16,901),
     1    MRJVMS(901),MRJKMS(16,901),
     1    MEDJVMS(8),MEDJKMS(16,8),SEGV(6,12),SEGK(6,12,16),
     1    OMJVCS(901),OMJKCS(16,901),
     1    IMJVCS(901),IMJKCS(16,901),
     1    MRJVCS(901),MRJKCS(16,901),
     1    MEDJVCS(8),MEDJKCS(16,8),VSEGV(8,4),VSEGK(8,4,16)
C
        INTEGER IG(250),IP(250),
     1   NSLC,KSLC(8),NTGF,NVOL,KVOL(15),NSOL,KSOL(15)
	DOUBLE PRECISION 
     1    PHI(250),GAMMA(250),DGAM(250),PPHI(250),
     1    PREMAX,MAXPHI,DELGAM(250),NDERIV(250,250),
     1    RPHI(250),RGAMMA(250),RDERIV(250,250)
C
C       PHI-    ERROR VECTOR- TO BE ZEROED VIA NEWTON ITERATION
C       GAMMA-  VARIABLE VECTOR
C       DGAM-   VARIABLE INCREMENT TO COMPUTE NUMERICAL DERIVATIVE
C       PPHI-   ERROR VECTOR AT GAMMA + DGAM
C       NDERIV- NUMERICAL DERIVATIVE OF PHI WITH RESPECT TO GAMMA
C       PREMAX- PRIOR MAXPHI TO ASSESS REUSE OF THE JACOBIAN
C       MAXPHI- MAXIMUM COMPONENT OF THE PHI VECTOR
C	DELGAM-	CORRECTIONS TO THE GAMMA VECTOR COMPUTED BY LES
C       RPHI-   REDUCED PHI WITH ENTRY J EQUAL TO PHI(IP(J))
C       RGAMMA- REDUCED GAMMA WITH ENTRY K EQUAL TO GAMMA(IG(K))
C       RDERIV- REDUCED JACOBIAN WITH ENTRY (J,K) EQUAL TO NDERIV(IP(J),IG(K))
C
C KEY TO INDICES:
C  6 nephrons: SF (INEPH = 1) plus 5 JM (INEPH = 2-6)
C    SFPCT and SFPST in series with SDHL
C    JMPCT and JMPST in series with LDHLu and LDHLl
C  LDHLu will be restricted to OM and LDHLl to IM
C
C Plan for 24000 SF nephrons (@30 nl/min -> 720 mul/min) 
C    12400 JM nephrons (@60 nl/min -> 720 mul/min) 
C     Of the 12000 JM, turns can parallel CD coning:
C    1 mm 6400
C    2 mm 3200
C    3 mm 1600
C    4 mm  800
C    5 mm  400
C
C Segment numbers (ISEG)
C      1- PCT
C      2- PST
C      3- sDHL, lDHLu
C      4- lDHLl
C      5- tAHL
C      6- AHLm
C      7- AHLc
C      8- DCT
C      9- CNT
C      10-CCD
C      11-OMCD
C      12-IMCD
C
C ARCHIVED NEPHRON VARIABLES
C
	INTEGER 
     1   RX(6,12), RCHOP(6,12)
C
C GENERAL PARAMETERS
	DOUBLE PRECISION
     1   RDIST(6,12), REPSI(6,12), RL0(6,12),
     1   RL(6,12,901,2), RKHY(6,12,5), RKDHY(6,12,5), RBCO2(6,12,3)
C
C LUMINAL AND PERITUBULAR PARAMETERS
	DOUBLE PRECISION
     1   RVM(6,12,901,2), RPM(6,12,901,2), RCM(6,12,15,901,2),
     1   RIMPM(6,12,901,2), RLCHM(6,12,901), RXM(6,12,15,901),
     1   RVS(6,12,901,2), RPS(6,12,901,2), RCS(6,12,15,901,2),
     1   RIMPS(6,12,901,2), RLCHS(6,12,901), RXS(6,12,15,901),
     1   RTL(6,12), RDX(6,12), RRM0(6,12),
     1   RMUM(6,12), RETA(6,12), RSM(6,12,901,2),
     1   RAM(6,12,901,2), RFVM(6,12,901,2), RFKM(6,12,16,901,2)
C
C INTERSPACE PARAMETERS
	DOUBLE PRECISION
     1   RAME(6,12), RAE0(6,12), RAE(6,12,901,2),
     1   RMUA(6,12), RCHVL0(6,12), RCHVL(6,12,901,2),
     1   RMUV(6,12), RLPME(6,12), RLPES(6,12),
     1   RSME(6,12,15), RSES(6,12,15), RHME(6,12,15),
     1   RHES(6,12,15), RVE(6,12,901,2), RPE(6,12,901,2),
     1   RCE(6,12,15,901,2), RLCHE(6,12,901), RXE(6,12,15,901),
     1   RFEVM(6,12,901,2), RFEKM(6,12,15,901,2), RFEVS(6,12,901,2),
     1   RFEKS(6,12,15,901,2), RCURE(6,12,901)
C
C CELL PARAMETERS
	DOUBLE PRECISION
     1   RAIE(6,12,3), RAMI(6,12,3), RAIS(6,12,3),
     1   RAI0(6,12,3), RCLVL0(6,12,3), RIMP0(6,12,3),
     1   RCLVL(6,12,3,901,2), RZIMP(6,12,3), RTBUF(6,12,3),
     1   RPKB(6,12,3), RCBUF(6,12,3,901,2), RHCBUF(6,12,3,901,2),
     1   RLPMI(6,12,3), RLPIS(6,12,3), RSMI(6,12,3,15),
     1   RSIS(6,12,3,15), RHMI(6,12,3,15), RHIS(6,12,3,15),
     1   RLMI(6,12,3,15,15), RLIS(6,12,3,15,15), RATMI(6,12,3,15,901),
     1   RATIS(6,12,3,15,901), RATIE(6,12,3,15,901), RVI(6,12,3,901,2),
     1   RPI(6,12,3,901,2), RCI(6,12,3,15,901,2), RIMP(6,12,3,901),
     1   RLCHI(6,12,3,901), RXI(6,12,3,15,901), RFIVM(6,12,3,901,2),
     1   RFIKM(6,12,3,15,901,2), RFIVS(6,12,3,901,2), 
     1   RFIKS(6,12,3,15,901,2), RCURI(6,12,3,901), 
     1   RJV(6,12,3,901,2), RJK(6,12,3,15,901,2)
C
C SPECIAL TRANSPORTERS
	CHARACTER*1 RISOFM(6,12)
	DOUBLE PRECISION
     1   RNP(6,12,3), RKNPN(6,12,3), RKNPK(6,12,3),
     1   RKNH4(6,12,3), RNPHK(6,12,3), RLHP(6,12,3),
     1   RXIHP(6,12,3), RXHP(6,12,3), RNAE1(6,12,3),
     1   RNTSC(6,12,3), RNNHE3(6,12,3), RJNAK(6,12,3,3,901,2),
     1   RJHK(6,12,3,901,2), RJHP(6,12,3,901,2), RJAE1(6,12,3,901,2),
     1   RJTSC(6,12,3,901,2), RJNHE3(6,12,3,3,901,2), RQIAMM(6,12),
     1   RNNKCC(6,12,3), RNKCL(6,12,3), RJNKCC(6,12,3,4,901,2),
     1   RJKCC(6,12,3,3,901,2)
C
C TORQUE AND COMPLIANCE VARIABLES
	DOUBLE PRECISION
     1   RTQM(6,12,501,2), RRM(6,12,501,2), RMUR(6,12), RSCALMI(6,12),
     1   RSCALIS(6,12), RVM0(6,12), RAM0(6,12), RTQM0(6,12),
     1   RLHP0(6,12,3), RNP0(6,12,3), RNNHE30(6,12,3), RLPMI0(6,12,3),
     1   RLPIS0(6,12,3), RHMI0(6,12,3,15), RHIS0(6,12,3,15),
     1   RLMI0(6,12,3,15,15), RLIS0(6,12,3,15,15), RQIAMM0(6,12),
     1   RMVL(6,12), RMVD(6,12), RRMT0(6,12)
C
C TGF VARIABLES: STARTING WITH THE BASELINE SNGFR, 
C   IT IS ADJUSTED ACCORDING TO MACULA DENSA CELL [CL-] 
C INITIAL PT PRESSURE IS ADJUSTED ACCORDING TO A DISTAL NEPHRON RESISTANCE
C
	DOUBLE PRECISION 
     1   FVM0(2),AFVM(6),DFVM(6),AFVM0(2),DFVM0(2),MDCL(2),DMDCL(2),
     1   TGGAM(15),TGPHI(15),TGPPHI(15),TGEPSI,TGDERIV(15,15),
     1   TGDGAM(15),TGDELG(15),MAXTGPHI,RTGPHI(15),DNR
C	
C	FVM0-	SNGFR READ FROM BOUND.DAT-INDEX FOR SF (1) OR JM (2)
C	AFVM-	FIXED COMPONENT OF FVM0
C	DFVM-	COMPONENT OF FVM0 THAT CAN BE MODULATED BY TGF
C	MDCL-	MACULA DENSA CL- REFERENCE CONCENTRATION
C	DMDCL-	MACULA DENSA CL- RANGE YIELDING DFVM
C	TGGAM-	TG FEEDBACK VARIABLE (CHANGE IN FV)
C               1-      SF FVM(1)/FVM0(1)
C               2-      SF PM(1)
C               2J-1-   JM FVM(1)/FVM0(2), J=2,6
C               2J-     JM PM(1)
C               13-     CNT PM(1)
C	TGPHI-	TG FEEDBACK RELATION
C               2J-1-   TGF RELATIONS, J=1,6
C         TGPHI(2J-1) = FVM(0) - AFVM - (DFVM / DMDCL) * [MDCL - CIAHL]
C               2J-     MATCHING PRESSURES AT CNT ENTRY 
C         TGPHI(2J) = PM-END-DCT = PM-CNT
C  xxx          13-     IMCD PM(L)=0 mmHg   - This did not work
C         TGPHI(13) = CNT PM(1) - DNR*FVMCNT
C	TGEPSI-	TG FEEDBACK ERROR
C	TGDERIV-	TG FEEDBACK JACOBIAN 
C	TGDGAM-	TG FEEDBACK VARIABLE INCREMENT
C	TGDELG-	TG FEEDBACK VARIABLE NEWTON CORRECTION
C
C VASC:
C Plan for 36000 OMDVR = 24000 sDVR + 12000 lDVR ~ about 3x JM glomeruli
C           7200 MRDVR (DVR from 24000 SF glom = 72000, so this is 10%)  
C Plan for 72000 OMAVR, 48000 IMAVR, 14400 MRAVR; AVR are 2x the DVR number
C  OM and IM DVR derive from 12000 JM nephrons; MR DVR from 2400 SF nephrons
C  Of the IM vessels, turns parallel those of the nephron
C         IMDVR       IMAVR
C    1 mm  6000       12000
C    2 mm  3000        6000
C    3 mm  1500        3000
C    4 mm   750        1500
C    5 mm   375         750
C
C PROGRAM TO BUILD A VAS RECTUM:
C   OUTER MEDULLARY 1,2, descending 1 or 2 mm; 
C   INNER MEDULLARY 3,7, descending 1,...,5 mm into IM;
C   AND MEDULLARY RAY, 8, 2 mm in length.
C  NOTE THAT OMDVR FEEDING IMDVR IS DIFFERENT FROM OMDVR WHICH TURNS IN OM
C
        DOUBLE PRECISION RADVR(8),RDNUM(8),RANUM(8)
C       RADVR- RATIO OF AVR NUMBER TO DVR NUMBER FOR EACH VESSEL
C       RDNUM- NUMBER OF EACH CLASS OF DVR
C       RANUM- NUMBER OF EACH CLASS OF AVR
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
C	
	DOUBLE PRECISION SW,PTOL
	INTEGER PR,PFL
C
C	SW-	ERROR SWITCH FOR LES
C	PTOL-	PIVOT TOLERANCE
C	PR-	NUMBER OF THE PIVOT ROW FOR EACH COLUMN
C	PFL-	PIVOT FLAG SHOWING ROWS FOR WHICH PIVOTS HAVE BEEN CHOSEN
C
	DIMENSION PR(250),PFL(250)
C
C
	COMMON /ARCH/
     1   RX, RCHOP,
     1   RDIST, REPSI, RL0, RL, RKHY, RKDHY, RBCO2,
     1   RVM, RPM, RCM, RIMPM, RLCHM, RXM,
     1   RVS, RPS, RCS, RIMPS, RLCHS, RXS,
     1   RTL, RDX, RRM0, RMUM, RETA, RSM, RAM,
     1   RFVM, RFKM,
     1   RAME, RAE0, RAE, RMUA, RCHVL0, RCHVL,
     1   RMUV, RLPME, RLPES, RSME, RSES,
     1   RHME, RHES, RVE, RPE,
     1   RCE, RLCHE, RXE, RFEVM, RFEKM,
     1   RFEVS, RFEKS, RCURE,
     1   RAIE, RAMI, RAIS, RAI0, RCLVL0, RIMP0,
     1   RCLVL, RZIMP, RTBUF, RPKB, RCBUF, RHCBUF,
     1   RLPMI, RLPIS, RSMI, RSIS, RHMI, RHIS,
     1   RLMI, RLIS, RATMI, RATIS, RATIE, RVI, RPI,
     1   RCI, RIMP, RLCHI, RXI, RFIVM, RFIKM,
     1   RFIVS, RFIKS, RCURI, RJV, RJK, RISOFM,
     1   RNP, RKNPN, RKNPK, RKNH4, RNPHK, RLHP,
     1   RXIHP, RXHP, RNAE1, RNTSC, RNNHE3,
     1   RJNAK, RJHK, RJHP, RJAE1, RJTSC, RJNHE3,
     1   RQIAMM, RNNKCC, RNKCL, RJNKCC, RJKCC,
     1   RTQM, RRM, RMUR, RSCALMI, RSCALIS, 
     1   RVM0, RAM0, RTQM0, RLHP0, RNP0, RNNHE30,
     1   RLPMI0, RLPIS0, RHMI0, RHIS0, RLMI0, 
     1   RLIS0, RQIAMM0, RMVL, RMVD, RRMT0
        COMMON /TGF/
     1   FVM0,AFVM,DFVM,AFVM0,DFVM0,MDCL,DMDCL,DNR
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
        COMMON /BOUND/
     1   OMCHOP,IMCHOP,MRCHOP,SLICE,
     1   CSMR,CSOM,CSIM,PSMR,PSOM,PSIM,
     1   TGEPSI,TGGAM,CSGAM,PSGAM,KEPSI,
     1   PS0,IMPM0,IMPS0,CM0,CS0, 
     1   FVC0,PC0,IMPC0,HCT0,CC0
C
        COMMON /KGEN/
     1    OMDIST,IMDIST,MRDIST,
     1    OMJVMS,OMJKMS,IMJVMS,IMJKMS,MRJVMS,MRJKMS,
     1    MEDJVMS,MEDJKMS,SEGV,SEGK,
     1    OMJVCS,OMJKCS,IMJVCS,IMJKCS,MRJVCS,MRJKCS,
     1    MEDJVCS,MEDJKCS,VSEGV,VSEGK
C
        DATA RVVS /57664*0.D0/,RVIMPS/57664*0.D0/
C
        OPEN (301,FILE='kguess.dat')
        OPEN (302,FILE='kparam.dat')
        OPEN (303,FILE='kguess.new')
        OPEN (304,FILE='kphi.dat')
        OPEN (305,FILE='kjac.dat')
C
	OPEN (10,FILE='nephron/lumen.dat')
	OPEN (11,FILE='nephron/blood.dat')
	OPEN (12,FILE='nephron/xy.co2')
	OPEN (13,FILE='nephron/tgf.dat')
	OPEN (14,FILE='nephron/guess.dat')
	OPEN (15,FILE='nephron/guess.new')
	OPEN (16,FILE='nephron/medjms.dat')
	OPEN (17,FILE='nephron/torque.dat')
	OPEN (18,FILE='nephron/xy.dat')
	OPEN (19,FILE='nephron/medjms.sum')
C
C
	OPEN (160,FILE='nephron/sfpct/sfpctparam.dat')
	OPEN (161,FILE='nephron/sfpct/sfpctresult.dat')
	OPEN (162,FILE='nephron/sfpct/sfpctbound.dat')
	OPEN (163,FILE='nephron/sfpct/sfpctptpick.dat')
	OPEN (164,FILE='nephron/sfpct/sfpctguess.dat')
	OPEN (165,FILE='nephron/sfpct/sfpctguess.new')
	OPEN (166,FILE='nephron/sfpct/sfpctptim.dat')
	OPEN (167,FILE='nephron/sfpct/sfpctfluxes.dat')
	OPEN (168,FILE='nephron/sfpct/sfpctxy.qo2')
C
	OPEN (150,FILE='nephron/sfpst/sfpstparam.dat')
	OPEN (151,FILE='nephron/sfpst/sfpstresult.dat')
	OPEN (152,FILE='nephron/sfpst/sfpstbound.dat')
	OPEN (153,FILE='nephron/sfpst/sfpstptpick.dat')
	OPEN (154,FILE='nephron/sfpst/sfpstguess.dat')
	OPEN (155,FILE='nephron/sfpst/sfpstguess.new')
	OPEN (156,FILE='nephron/sfpst/sfpstptim.dat')
	OPEN (157,FILE='nephron/sfpst/sfpstfluxes.dat')
	OPEN (158,FILE='nephron/sfpst/sfpstxy.qo2')
C
	OPEN (140,FILE='nephron/jmpct/jmpctparam.dat')
	OPEN (141,FILE='nephron/jmpct/jmpctresult.dat')
	OPEN (142,FILE='nephron/jmpct/jmpctbound.dat')
	OPEN (143,FILE='nephron/jmpct/jmpctptpick.dat')
	OPEN (144,FILE='nephron/jmpct/jmpctguess.dat')
	OPEN (145,FILE='nephron/jmpct/jmpctguess.new')
	OPEN (146,FILE='nephron/jmpct/jmpctptim.dat')
	OPEN (147,FILE='nephron/jmpct/jmpctfluxes.dat')
	OPEN (148,FILE='nephron/jmpct/jmpctxy.qo2')
C
	OPEN (130,FILE='nephron/jmpst/jmpstparam.dat')
	OPEN (131,FILE='nephron/jmpst/jmpstresult.dat')
	OPEN (132,FILE='nephron/jmpst/jmpstbound.dat')
	OPEN (133,FILE='nephron/jmpst/jmpstptpick.dat')
	OPEN (134,FILE='nephron/jmpst/jmpstguess.dat')
	OPEN (135,FILE='nephron/jmpst/jmpstguess.new')
	OPEN (136,FILE='nephron/jmpst/jmpstptim.dat')
	OPEN (137,FILE='nephron/jmpst/jmpstfluxes.dat')
	OPEN (138,FILE='nephron/jmpst/jmpstxy.qo2')
C
	OPEN (120,FILE='nephron/sdhl/sdhlparam.dat')
	OPEN (121,FILE='nephron/sdhl/sdhlresult.dat')
	OPEN (122,FILE='nephron/sdhl/sdhlbound.dat')
	OPEN (123,FILE='nephron/sdhl/sdhlptpick.dat')
	OPEN (124,FILE='nephron/sdhl/sdhlguess.dat')
	OPEN (125,FILE='nephron/sdhl/sdhlguess.new')
	OPEN (126,FILE='nephron/sdhl/sdhlptim.dat')
	OPEN (127,FILE='nephron/sdhl/sdhlfluxes.dat')
C
	OPEN (110,FILE='nephron/ldhlu/ldhluparam.dat')
	OPEN (111,FILE='nephron/ldhlu/ldhluresult.dat')
	OPEN (112,FILE='nephron/ldhlu/ldhlubound.dat')
	OPEN (113,FILE='nephron/ldhlu/ldhluptpick.dat')
	OPEN (114,FILE='nephron/ldhlu/ldhluguess.dat')
	OPEN (115,FILE='nephron/ldhlu/ldhluguess.new')
	OPEN (116,FILE='nephron/ldhlu/ldhluptim.dat')
	OPEN (117,FILE='nephron/ldhlu/ldhlufluxes.dat')
C
	OPEN (100,FILE='nephron/ldhll/ldhllparam.dat')
	OPEN (101,FILE='nephron/ldhll/ldhllresult.dat')
	OPEN (102,FILE='nephron/ldhll/ldhllbound.dat')
	OPEN (103,FILE='nephron/ldhll/ldhllptpick.dat')
	OPEN (104,FILE='nephron/ldhll/ldhllguess.dat')
	OPEN (105,FILE='nephron/ldhll/ldhllguess.new')
	OPEN (106,FILE='nephron/ldhll/ldhllptim.dat')
	OPEN (107,FILE='nephron/ldhll/ldhllfluxes.dat')
C
	OPEN (90,FILE='nephron/tahl/tahlparam.dat')
	OPEN (91,FILE='nephron/tahl/tahlresult.dat')
	OPEN (92,FILE='nephron/tahl/tahlbound.dat')
	OPEN (93,FILE='nephron/tahl/tahlptpick.dat')
	OPEN (94,FILE='nephron/tahl/tahlguess.dat')
	OPEN (95,FILE='nephron/tahl/tahlguess.new')
	OPEN (96,FILE='nephron/tahl/tahlptim.dat')
	OPEN (97,FILE='nephron/tahl/tahlfluxes.dat')
C
	OPEN (80,FILE='nephron/ahlm/ahlmparam.dat')
	OPEN (81,FILE='nephron/ahlm/ahlmresult.dat')
	OPEN (82,FILE='nephron/ahlm/ahlmbound.dat')
	OPEN (83,FILE='nephron/ahlm/ahlmptpick.dat')
	OPEN (84,FILE='nephron/ahlm/ahlmguess.dat')
	OPEN (85,FILE='nephron/ahlm/ahlmguess.new')
	OPEN (86,FILE='nephron/ahlm/ahlmptim.dat')
	OPEN (87,FILE='nephron/ahlm/ahlmfluxes.dat')
C
	OPEN (70,FILE='nephron/ahlc/ahlcparam.dat')
	OPEN (71,FILE='nephron/ahlc/ahlcresult.dat')
	OPEN (72,FILE='nephron/ahlc/ahlcbound.dat')
	OPEN (73,FILE='nephron/ahlc/ahlcptpick.dat')
	OPEN (74,FILE='nephron/ahlc/ahlcguess.dat')
	OPEN (75,FILE='nephron/ahlc/ahlcguess.new')
	OPEN (76,FILE='nephron/ahlc/ahlcptim.dat')
	OPEN (77,FILE='nephron/ahlc/ahlcfluxes.dat')
C
	OPEN (60,FILE='nephron/dct/dcparam.dat')
	OPEN (61,FILE='nephron/dct/dcresult.dat')
	OPEN (62,FILE='nephron/dct/dcbound.dat')
	OPEN (63,FILE='nephron/dct/dcptpick.dat')
	OPEN (66,FILE='nephron/dct/dcptim.dat')
	OPEN (64,FILE='nephron/dct/dcguess.dat')
	OPEN (65,FILE='nephron/dct/dcguess.new')
	OPEN (67,FILE='nephron/dct/dcfluxes.dat')
C
	OPEN (50,FILE='nephron/cnt/cnparam.dat')
	OPEN (51,FILE='nephron/cnt/cnresult.dat')
	OPEN (52,FILE='nephron/cnt/cnbound.dat')
	OPEN (53,FILE='nephron/cnt/cnptpick.dat')
	OPEN (54,FILE='nephron/cnt/cnguess.dat')
	OPEN (55,FILE='nephron/cnt/cnguess.new')
	OPEN (56,FILE='nephron/cnt/cnptim.dat')
	OPEN (57,FILE='nephron/cnt/cnfluxes.dat')
C
	OPEN (40,FILE='nephron/cct/ccparam.dat')
	OPEN (41,FILE='nephron/cct/ccresult.dat')
	OPEN (42,FILE='nephron/cct/ccbound.dat')
	OPEN (43,FILE='nephron/cct/ccptpick.dat')
	OPEN (44,FILE='nephron/cct/ccguess.dat')
	OPEN (45,FILE='nephron/cct/ccguess.new')
	OPEN (46,FILE='nephron/cct/ccptim.dat')
	OPEN (47,FILE='nephron/cct/ccfluxes.dat')
C
	OPEN (30,FILE='nephron/ompct/omparam.dat')
	OPEN (31,FILE='nephron/ompct/omresult.dat')
	OPEN (32,FILE='nephron/ompct/ombound.dat')
	OPEN (33,FILE='nephron/ompct/omptpick.dat')
	OPEN (34,FILE='nephron/ompct/omguess.dat')
	OPEN (35,FILE='nephron/ompct/omguess.new')
	OPEN (36,FILE='nephron/ompct/omptim.dat')
	OPEN (37,FILE='nephron/ompct/omfluxes.dat')
C
	OPEN (20,FILE='nephron/imct/imparam.dat')
	OPEN (21,FILE='nephron/imct/imresult.dat')
	OPEN (22,FILE='nephron/imct/imbound.dat')
	OPEN (23,FILE='nephron/imct/imptpick.dat')
	OPEN (24,FILE='nephron/imct/imguess.dat')
	OPEN (25,FILE='nephron/imct/imguess.new')
	OPEN (26,FILE='nephron/imct/imptim.dat')
	OPEN (27,FILE='nephron/imct/imfluxes.dat')
C
C
	OPEN (290,FILE='vasc/vrnum.dat')
	OPEN (291,FILE='vasc/vasc.dat')
	OPEN (292,FILE='vasc/xy.fv')
        OPEN (293,FILE='vasc/medjcs.dat')
        OPEN (294,FILE='vasc/medjcs.sum')
C
	OPEN (UNIT=201,FILE='vasc/omdvr/omdparam.dat')
	OPEN (UNIT=202,FILE='vasc/omdvr/omdbound.dat')
	OPEN (UNIT=203,FILE='vasc/omdvr/omdresult.dat')
	OPEN (UNIT=204,FILE='vasc/omdvr/omdend.dat')
C
	OPEN (UNIT=211,FILE='vasc/omavr/omaparam.dat')
	OPEN (UNIT=212,FILE='vasc/omavr/omabound.dat')
	OPEN (UNIT=213,FILE='vasc/omavr/omaresult.dat')
	OPEN (UNIT=214,FILE='vasc/omavr/omaend.dat')
C
	OPEN (UNIT=221,FILE='vasc/imdvr/imdparam.dat')
	OPEN (UNIT=222,FILE='vasc/imdvr/imdbound.dat')
	OPEN (UNIT=223,FILE='vasc/imdvr/imdresult.dat')
	OPEN (UNIT=224,FILE='vasc/imdvr/imdend.dat')
C
	OPEN (UNIT=231,FILE='vasc/imavr/imaparam.dat')
	OPEN (UNIT=232,FILE='vasc/imavr/imabound.dat')
	OPEN (UNIT=233,FILE='vasc/imavr/imaresult.dat')
	OPEN (UNIT=234,FILE='vasc/imavr/imaend.dat')
C
	OPEN (UNIT=241,FILE='vasc/mrdvr/mrdparam.dat')
	OPEN (UNIT=242,FILE='vasc/mrdvr/mrdbound.dat')
	OPEN (UNIT=243,FILE='vasc/mrdvr/mrdresult.dat')
	OPEN (UNIT=244,FILE='vasc/mrdvr/mrdend.dat')
C
	OPEN (UNIT=251,FILE='vasc/mravr/mraparam.dat')
	OPEN (UNIT=252,FILE='vasc/mravr/mrabound.dat')
	OPEN (UNIT=253,FILE='vasc/mravr/mraresult.dat')
	OPEN (UNIT=254,FILE='vasc/mravr/mraend.dat')
C
	OPEN (UNIT=261,FILE='vasc/oidvr/oidparam.dat')
	OPEN (UNIT=262,FILE='vasc/oidvr/oidbound.dat')
	OPEN (UNIT=263,FILE='vasc/oidvr/oidresult.dat')
	OPEN (UNIT=264,FILE='vasc/oidvr/oidend.dat')
C
	OPEN (UNIT=271,FILE='vasc/oiavr/oiaparam.dat')
	OPEN (UNIT=272,FILE='vasc/oiavr/oiabound.dat')
	OPEN (UNIT=273,FILE='vasc/oiavr/oiaresult.dat')
	OPEN (UNIT=274,FILE='vasc/oiavr/oiaend.dat')
C
	SOL(1)='  NA '
	SOL(2)='  K  '
	SOL(3)='  CL '
	SOL(4)=' HCO3'
	SOL(5)='H2CO3'
	SOL(6)=' CO2 '
	SOL(7)=' HPO4'
	SOL(8)='H2PO4'
	SOL(9)=' UREA'
	SOL(10)=' NH3 '
	SOL(11)=' NH4+'
	SOL(12)='  H+ '
	SOL(13)='HCO2-'
	SOL(14)='H2CO2'
	SOL(15)=' GLUC'
	SOL(16)=' PROT'
C
	T=1
	TIME=0.D0
C
C ASSIGN MODEL PARAMETERS, INCLUDING CORTICAL IS PRESSURE AND SOLUTE CONCENTRATIONS
C  CORTICAL CONDITIONS ARE READ AS PSMR(1,1) AND CSMR(I,1,1)
C    WHICH ARE ALSO EQUAL TO PSOM(1,1) AND CSOM(I,1,1)
C
        CALL KBSET
C
        SOLS=15
        SLICE=8
        NN=13+SLICE+SLICE*SOLS
	PTOL=1.D-15
	M1=250
	IERR=0
        NIL='    --     '
        DO 2 ISOL=1,SOLS+1
        DO 2 ISLC=1,SLICE
    2   VMAT(ISOL,ISLC)=NIL
C
        CALL KGAMINIT(GAMMA)
        DO 10 K=1,NN
   10   DGAM(K)=0.1D-1*GAMMA(K)
C
C KGAM PROVIDES TRANSLATION FROM SLICE GUESSES TO BOUNDARY DATA
C  FOR MEDULLARY AND MEDULLARY RAY STRUCTURES
C  PRESSURE AND FLOW CHANGES OF TGGAM WILL BE TRANSLATED WITHIN NEPHRON
C
        CALL KGAM(GAMMA,1)
        CALL KIS
C
C THERE WAS ENOUGH INFORMATION TO DETERMINE ALL OF THE INTERSTITIAL CONDITIONS
C  THROUGHOUT THE WHOLE KIDNEY MODEL (TUBULES AND VESSELS).
C THIS WAS DONE IN KGAM, RATHER THAN IN INDIVIDUAL TUBULE AND VESSEL ROUTINES
C  --AVOIDS INCONSISTENCIES AMONG THE 20 SUBROUTINES
C BOUNDARY DATA WILL BE ESTABLISHED IN THE FLOW.F ROUTINE, RATHER THAN IN *NEWT.F
C
C VASCULAR PARAMETERS
C
        READ (290,25) (RADVR(I),RDNUM(I),I=1,8)
   25   FORMAT(2D12.4)
        DO 26 I=1,8
   26   RANUM(I)=RDNUM(I)*RADVR(I)
C
C EXPERIMENTS WILL CYCLE FROM THIS POINT.
C
C	PRINT 30
   30   FORMAT(' NUMBER OF EXPERIMENTS = ',$)
C	READ *, EXP
        EXP=1
C
C ALLOW FOR 'SUBPROBLEMS" IN WHICH A SUBSET OF THE NN GAMMA AND PHI ARE SELECTED
	CALL KREDUCE(IN,IG,IP,NSLC,KSLC,NTGF,NVOL,KVOL,NSOL,KSOL)
C  
	DO 600 COUNT=1,EXP
C
C TG FEEDBACK SETTINGS: (Set outside proximal tubule since they involve AHL)
        READ(13,35) TGEPSI,DNR,
     1   (AFVM0(J),DFVM0(J),MDCL(J),DMDCL(J),J=1,2)
   35   FORMAT(2D12.4,/,4D12.4,/,4D12.4)
C
        CALL NEPHRON(TGPHI,COUNT,0)
        CALL GEN(1)
C
C CHECK THE VASCULAR FLUXES
C
        CALL VASC(COUNT,0)
        CALL VGEN(1)
C
C CHECK THE ERROR VECTOR
C
        CALL KERR(TGPHI,PHI)
C
        DO 39 ISOL=1,SOLS
   39   ALABL(ISOL)=SOL(ISOL)
        ALABL(4)=' TCO2'
        ALABL(5)=' PKC '
        ALABL(6)='CO2HY'
        ALABL(7)=' TPO4'
        ALABL(8)=' PKP '
        ALABL(10)=' TAMM'
        ALABL(11)=' PKN '
        ALABL(13)='TFORM'
        ALABL(14)=' PKF '
C
        ITER=1
        WRITE(304,40) ITER
   40   FORMAT(/,2X,'ITER=',I3)
        PRINT 40, ITER
C
        IF(NTGF.EQ.0) THEN
        PRINT 102, (NIL,ITGF=1,13)
        WRITE(304,102) (NIL,ITGF=1,13)
  102   FORMAT(7X,8A12,/,7X,5A12)
        ELSE
        PRINT 104, (TGGAM(ITGF),ITGF=1,13)
        WRITE(304,104) (TGGAM(ITGF),ITGF=1,13)
  104   FORMAT(7X,8D12.4,/,7X,5D12.4)
        ENDIF
C
C FIRST PRINT ALL OF THE IG(XX) VARIABLES INTO THE VMAT MATRIX
C   VOLUME ALONG THE FIRST ROW, THEN SOLUTES IN ROWS 2-16
C
        IF(NSLC.EQ.0) GO TO 150
        IF(NVOL.EQ.0) GO TO 120
        DO 114 IVOL=1,NVOL
  114   WRITE(VMAT(1,KVOL(IVOL)),105) GAMMA(13+KVOL(IVOL))
  105   FORMAT(D12.4)
C
  120   IF(NSOL.EQ.0) GO TO 150
        DO 130 ISLC=1,NSLC
        DO 125 ISOL=1,NSOL
  125   WRITE(VMAT(1+KSOL(ISOL),KSLC(ISLC)),105) 
     1    GAMMA(IG(NTGF+NVOL+(ISLC-1)*NSOL+ISOL))
  130   CONTINUE
C
  150   CONTINUE
        PRINT 152, (ISLC, ISLC=1,SLICE)
        WRITE(304,152) (ISLC, ISLC=1,SLICE)
  152   FORMAT(/,2X,'SLICE',4X,8(I3,9X),/)
        PRINT 154, (VMAT(1,ISLC), ISLC=1,SLICE)
        WRITE(304,154) (VMAT(1,ISLC), ISLC=1,SLICE)
  154   FORMAT(' VOL ',2X,8A12)
        DO 155 ISOL=1,SOLS
        PRINT 156, SOL(ISOL),(VMAT(ISOL+1,ISLC), ISLC=1,SLICE)
        WRITE(304,156) SOL(ISOL),(VMAT(ISOL+1,ISLC), ISLC=1,SLICE)
  155   CONTINUE
  156   FORMAT(A5,2X,8A12)
C
        WRITE(304,41) (TGPHI(K),K=1,13)
   41   FORMAT(/,7X,8D12.4,/,7X,5D12.4)
        WRITE(304,42) (ISLC,ISLC=1,SLICE)
   42   FORMAT(/,2X,'SLICE',4X,8(I3,9X),/)
        WRITE(304,43) (PHI(K),K=14,13+SLICE)
   43   FORMAT(' VOL ',2X,8D12.4)
        DO 45 ISOL=1,SOLS
   45   WRITE(304,46) ALABL(ISOL),
     1   (PHI(21+(ISLC-1)*SOLS+ISOL),ISLC=1,SLICE)
   46   FORMAT(A5,2X,8D12.4)
C
        PRINT 41, (TGPHI(K),K=1,13)
        PRINT 42, (K,K=1,SLICE)
        PRINT 43, (PHI(K),K=14,13+SLICE)
        DO 48 ISOL=1,SOLS
   48   PRINT 46, ALABL(ISOL),
     1   (PHI(21+(ISLC-1)*SOLS+ISOL),ISLC=1,SLICE)
   49   FORMAT(5X,8D12.4)
C
C	
	MAXPHI=DABS(PHI(IP(1)))
	MAXJ=IP(1)
	DO 110 J=1,IN
	IF (MAXPHI.LT.DABS(PHI(IP(J)))) MAXJ=IP(J)
  110   IF (MAXPHI.LT.DABS(PHI(IP(J)))) MAXPHI=DABS(PHI(IP(J)))
        WRITE(304,112) MAXJ, MAXPHI
	PRINT 112, MAXJ, MAXPHI
  112   FORMAT (1X,I5,D12.4)
	IF(MAXPHI.LT.KEPSI) GO TO 401
C
C THE PHI ARE TOO LARGE.
C   WE WILL DETERMINE THE DERIVATIVE OF THE PHI WITH RESPECT TO
C THE VARIABLES ,GAMMA, AND STORE THESE IN THE ARRAY NDERIV( , ).
C
  199   ITER=ITER+1
	DO  200 K=1,IN
	GAMMA(IG(K))=GAMMA(IG(K))+DGAM(IG(K))
        CALL KGAM(GAMMA,1)
        CALL KIS
C
	CALL NEPHRON(TGPHI,COUNT,2)
        CALL GEN(0)
        CALL VASC(COUNT,2)
        CALL VGEN(0)
        CALL KERR(TGPHI,PPHI)
C
	DO 210 J=1,NN
  210   NDERIV(J,IG(K))=(PPHI(J)-PHI(J))/DGAM(IG(K))
C
	GAMMA(IG(K))=GAMMA(IG(K))-DGAM(IG(K))
        CALL KGAM(GAMMA,1)
        CALL KIS
        PRINT 201
  201   FORMAT('.',$)
  200   CONTINUE
C
C THE JACOBIAN HAS BEEN COMPUTED NUMERICALLY. THE CURRENT GUESSES
C SIT IN GAMMA AND THE ERRORS IN PHI.  IT REMAINS TO COMPUTE
C NDERIV-1(PHI) AT THE CORRECTION TO BE SUBTRACTED FROM GAMMA.
C  FIRST FORMULATE THE REDUCED PROBLEM
C
  219   DO 220 K=1,IN
  220   RGAMMA(K)=GAMMA(IG(K))
        DO 222 J=1,IN
  222   RPHI(J)=PHI(IP(J))
        DO 225 K=1,IN
        DO 225 J=1,IN
  225   RDERIV(J,K)=NDERIV(IP(J),IG(K))
C
        DO 227 J=1,IN
  227   WRITE(305,228) (RDERIV(J,K),K=1,IN)
  228   FORMAT(6D20.12)
C
	CALL LES8(RDERIV,RPHI,PR,PFL,IN,SW,PTOL,DELGAM,M1)
	IF (SW) 230,600,230
C
  230   CONTINUE
        DSCAL=1.D0-0.00D0**(ITER-1)
C
        DO 232 K=1,IN
  232   RGAMMA(K)=RGAMMA(K)-DSCAL*DELGAM(K)
C  232   RGAMMA(K)=RGAMMA(K)-DELGAM(K)
C
C        PRINT 49, (DELGAM(K),K=1,IN)
C        PRINT 49, (RGAMMA(K),K=1,IN)

        DO 235 K=1,IN
  235   GAMMA(IG(K))=RGAMMA(K)
        CALL KGAM(GAMMA,1)
        CALL KIS
C
C RETEST THIS CORRECTION
C
	CALL NEPHRON(TGPHI,COUNT,2)
        CALL GEN(1)
        CALL VASC(COUNT,2)
        CALL VGEN(1)
        CALL KERR(TGPHI,PHI)
C
        WRITE(304,40) ITER
        PRINT 40, ITER
        IF(NTGF.EQ.0) THEN
        PRINT 102, (NIL,ITGF=1,13)
        WRITE(304,102) (NIL,ITGF=1,13)
        ELSE
        PRINT 104, (TGGAM(ITGF),ITGF=1,13)
        WRITE(304,104) (TGGAM(ITGF),ITGF=1,13)
        ENDIF
C
C FIRST PRINT ALL OF THE IG(XX) VARIABLES INTO THE VMAT MATRIX
C   VOLUME ALONG THE FIRST ROW, THEN SOLUTES IN ROWS 2-16
C
        IF(NSLC.EQ.0) GO TO 350
        IF(NVOL.EQ.0) GO TO 320
        DO 314 IVOL=1,NVOL
  314   WRITE(VMAT(1,KVOL(IVOL)),105) GAMMA(13+KVOL(IVOL))
C
  320   IF(NSOL.EQ.0) GO TO 350
        DO 330 ISLC=1,NSLC
        DO 325 ISOL=1,NSOL
  325   WRITE(VMAT(1+KSOL(ISOL),KSLC(ISLC)),105) 
     1    GAMMA(IG(NTGF+NVOL+(ISLC-1)*NSOL+ISOL))
  330   CONTINUE
C
  350   CONTINUE
        PRINT 152, (ISLC, ISLC=1,SLICE)
        WRITE(304,152) (ISLC, ISLC=1,SLICE)
        PRINT 154, (VMAT(1,ISLC), ISLC=1,SLICE)
        WRITE(304,154) (VMAT(1,ISLC), ISLC=1,SLICE)
        DO 355 ISOL=1,SOLS
        PRINT 156, SOL(ISOL),(VMAT(ISOL+1,ISLC), ISLC=1,SLICE)
        WRITE(304,156) SOL(ISOL),(VMAT(ISOL+1,ISLC), ISLC=1,SLICE)
  355   CONTINUE
C
        WRITE(304,41) (TGPHI(K),K=1,13)
        WRITE(304,42) (K,K=1,SLICE)
        WRITE(304,43) (PHI(K),K=14,13+SLICE)
        DO 245 ISOL=1,SOLS
  245   WRITE(304,46) ALABL(ISOL),
     1   (PHI(21+(ISLC-1)*SOLS+ISOL),ISLC=1,SLICE)
C
        PRINT 41, (TGPHI(K),K=1,13)
        PRINT 42, (K,K=1,SLICE)
        PRINT 43, (PHI(K),K=14,13+SLICE)
        DO 248 ISOL=1,SOLS
  248   PRINT 46, ALABL(ISOL),
     1   (PHI(21+(ISLC-1)*SOLS+ISOL),ISLC=1,SLICE)
C
        PREMAX=MAXPHI
	MAXPHI=DABS(PHI(IP(1)))
	MAXJ=IP(1)
	DO 250 J=1,IN
	IF (MAXPHI.LT.DABS(PHI(IP(J)))) MAXJ=IP(J)
  250   IF (MAXPHI.LT.DABS(PHI(IP(J)))) MAXPHI=DABS(PHI(IP(J)))
        WRITE(304,112) MAXJ, MAXPHI
	PRINT 112, MAXJ, MAXPHI
	IF(MAXPHI.LT.KEPSI) GO TO 401
C
	IF (ITER.LT.20) THEN
        IF (MAXPHI.LE.0.1D0*PREMAX) GO TO 219
        IF (MAXPHI.GT.0.1D0*PREMAX) GO TO 199
        ELSE
        ENDIF
C
C IF THE SPATIAL ITERATION FAILS TO CONVERGE THE PROGRAM STOPS.
	PRINT 280
  280   FORMAT (' TOO MANY ITERATIONS IN KIDNEY')
	STOP
C
C THE SPATIAL STEP CONVERGED,AND THE PROGRAM PICKS UP AT THIS POINT:
C
  401   CONTINUE
        CALL RESULT(COUNT)
        CALL VRESULT(COUNT)
C
        WRITE (15,420) (TGGAM(J),J=1,13)
  420   FORMAT(D30.20)
        WRITE (303,428) (PSGAM(K),K=1,SLICE+2)
        DO 472 J=1,11
  472   WRITE (303,428) (CSGAM(J,K),K=1,SLICE+2)
        DO 474 J=13,15
  474   WRITE (303,428) (CSGAM(J,K),K=1,SLICE+2)
  428   FORMAT(10D16.8)
C
  600   CONTINUE
	STOP
	END
