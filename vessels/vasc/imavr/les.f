C
      SUBROUTINE LES(A,B,PR,PFL,M,SW,PTOL,C,M1)
C     THIS IS A LINEAR EQUATION SOLVER ROUTINE.  IT USES GAUSSIAN
C     ELIMINATION WITH PARTIAL PIVOTING AND BACK SUBSTITUTION TO SOLVE
C     THE SYSTEM OF EQUATIONS AX=B.  IF A SMALL PIVOT OCCURS THEN SW IS
C     SET TO 0 AND A SMALL PIVOT DIAGNOSTIC IS PRINTED.  A IS THE M X M
C     MATRIX, B IS THE RIGHT HAND SIDE.  PR AND PFL ARE TWO M
C     VECTORS OF INTEGERS.  M IS THE SIZE OF THE MATRIX.  SW IS 1 IF
C     NORMAL RETURN.  PTOL IS PIVOT TOLERANCE.  THE SUBROUTINE RETURNS
C     WITH THE SOLUTION IN C.  THE INITIAL CONTENTS OF A AND B ARE
C     CHANGED BY THE SUBROUTINE.  PR CONTAINS THE NUMBER OF THE PIVOT
C     ROW FOR EACT COLUMN.  PFL IS A PIVOT FLAG VECTOR AND INDICATES
C     THE ROWS IN WHICH THE PIVOTS HAVE ALREADY BEEN CHOSEN.  IN BACK
C     SUBSTITUTION PFL HAS ZEROS IN PIVOTS FOR WHICH THE BACK
C     SUBSTITUTION HAS ALREADY BEEN DONE.
C
      REAL*8 SW,PTOL,T2,T3,T4
        REAL*8 A,B,C
	DIMENSION A(M1,M1),B(M1),C(M1)
      INTEGER PR(1),PFL(1)
      SW=1.0
      DO 10 I=1,M
      PR(I)=0
 10   PFL(I)=0
C     REDUCE A TO A PERMUTATION OF UPPER TRIANGULAR FROM USING
C     PARTIAL PIVOTING
C
      DO 50 K=1,M
      KCOL=K
      T2=PTOL
      IR=0
C
C     FIND THE ROW NUMBER 'IR' OF MAX A(I,K) IN NON-PIVOT ROW
C
      DO 15 I=1,M
      IF (PFL(I).NE.0) GOTO 15
      T3=ABS(A(I,K))
      IF (T3.LE.T2) GOTO 15
      T2=T3
      IR=I
 15   CONTINUE
C
C     IF ALL A(I,K) TOO SMALL, PRINT DIAGNOSTIC, SET SWITCH TO ZERO
C
      IF (IR.NE.0) GOTO 20
C     WRITE (16,30) KCOL
        PRINT 30,KCOL
 30   FORMAT (' ','MATRIX SINGULAR, COLUMN',I4,'ZERO')
      SW=0.000
      GOTO 999
C     RECORD PIVOT ROW IR IN PR(K) AND SET PFL(IR)=1
C
 20   PR(K)=IR
      PFL(IR)=1
C
C     PERFORM PIVOTING WITH A(IR,K) AS THE PIVOT ON A,B
C
      KPLO=K+1
      T3=1.000/A(IR,K)
      IF (KPLO.GT.M) GOTO 26
      DO 25 J=KPLO,M
 25   A(IR,J)=A(IR,J)*T3
 26   B(IR)=B(IR)*T3
      A(IR,K)=1.000
C
C     NOW MAKE ALL ELEMENTS OF K-TH COLUMN IN NON-FLAGGED ROWS ZERO
C     BY ROW OPERATIONS.  SKIP PIVOT ROW IR.  SKIP FLAGGED ROWS
C
      DO 50 I=1,M
      IF (PFL(I).NE.0) GOTO 50
      IF (KPLO.GT.M) GOTO 46
      DO 45 J=KPLO,M
 45   A(I,J)=A(I,J)-A(I,K)*A(IR,J)
 46   B(I)=B(I)-A(I,K)*B(IR)
      A(I,K)=0.000
 50   CONTINUE
C
C     BACK SUBSTITUTION ON A,B TO GET THE SOLUTION.  KC IS COMPLEMENT
C     OF COLUMN NUMBERS, VARIES FROM M TO 1.
C
      DO 70 K=1,M
      KC=M-K+1
      IR=PR(KC)
      PFL(IR)=0
      DO 60 I=1,M
      IF (PFL(I).EQ.0) GOTO 60
      B(I)=B(I)-B(IR)*A(I,KC)
 60   CONTINUE
      C(KC)=B(IR)
 70   CONTINUE
      DO 80 J=1,M
 80   B(J)=C(J)
999   RETURN
      END
