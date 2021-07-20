
!========================================================================
!
!                          T o m o f a s t - x
!                        -----------------------
!
!           Authors: Vitaliy Ogarko, Jeremie Giraud, Roland Martin.
!
!               (c) 2021 The University of Western Australia.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

! --------------------------------------------------------------------------
! Program to calculate the first zeroes (root abscissas) of the first
! kind Bessel function of integer order N using the subroutine ROOTJ.
! --------------------------------------------------------------------------
! SAMPLE RUN:
!
! (Calculate the first 10 zeroes of 1st kind Bessel function of order 2).
!
! Zeroes of Bessel Function of order:  2
!
! Number of calculated zeroes: 10
!
! Table of root abcissas (5 items per line)
! 0.51356223D+01 0.84172442D+01 0.11619841D+02 0.14795952D+02 0.17959819D+02
! 0.21116997D+02 0.24270112D+02 0.27420574D+02 0.30569204D+02 0.33716520D+02
!
! Table of error codes (5 items per line)
!   0   0   0   0   0
!   0   0   0   0   0
!
! --------------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                                  F90 Release 1.0 By J-P Moreau, Paris
!                                           (www.jpmoreau.fr)
!---------------------------------------------------------------------------
!PROGRAM TROOTJ
!
!REAL*8  JZ(10)
!INTEGER IE(10)
!
!  N=2
!  NR=10
!
!  CALL ROOTJ(N,NR,JZ,IE)
!
!  print *,' '
!  write(*,10) N
!  write(*,20) NR
!  print *,' '
!  print *,'Table of root abcissas (5 items per line)'
!  write(*,50) (JZ(i),i=1,NR)
!  print *,' '
!  print *,'Table of error codes (5 items per line)'
!  write(*,60) (IE(i),i=1,NR)
!  print *,' '
!
!  stop
!
!10 format(/' Zeroes of Bessel Function of order: ', I2)
!20 format(/' Number of calculated zeroes: ', I2)
!
!50 format(5D15.8)
!60 format(5I4)
!
!END

module trootj

implicit none

contains

      SUBROUTINE ROOTJ(N,NK,JZERO,IER)
! ----------------------------------------------------------------------
!     CALCULATE THE FIRST NK ZEROES OF BESSEL FUNCTION J(N,X)
!
!     INPUTS:
!       N    ORDER OF FUNCTION J (INTEGER >= 0)                  I*4
!       NK   NUMBER OF FIRST ZEROES  (INTEGER > 0)               I*4
!     OUTPUTS:
!       JZERO(NK)  TABLE OF FIRST ZEROES (ABCISSAS)              R*8
!       IER(NK)    TABLE OF ERROR CODES (MUST BE ZEROES)         I*4
!
!     REFERENCE :
!     ABRAMOWITZ M. & STEGUN IRENE A.
!     HANDBOOK OF MATHEMATICAL FUNCTIONS
! ----------------------------------------------------------------------
      INTEGER N,NK,IERROR,K,NITMX
      REAL(KIND=8) ZEROJ,B0,B1,B2,B3,B5,B7,T0,T1,T3,T5,T7,PI,FN,FK,  &
      C1,C2,C3,C4,C5,F1,F2,F3,TOL,ERRJ,JZERO(NK)
      INTEGER IER(NK)
      DATA PI/3.14159265358979D0/,TOL/1.D-8/,NITMX/10/
      DATA C1,C2,C3,C4,C5 /1.8557571D0,1.033150D0,.00397D0,.0908D0,.043D0/
      FN = FLOAT(N)

!     FIRST ZERO
      IF (N==0) THEN

        ZEROJ = C1+C2-C3-C4+C5
!       WRITE(*,'(1X,A,I5,E15.6)') 'K=1,N=0,ZEROJ',K,ZEROJ

        CALL  SECANT(N,NITMX,TOL,ZEROJ,IERROR)

        IER(1)=IERROR
        JZERO(1)=ZEROJ
      ELSE

        F1 = FN**(1.D0/3.D0)
        F2 = F1*F1*FN
        F3 = F1*FN*FN
        ZEROJ = FN+C1*F1+(C2/F1)-(C3/FN)-(C4/F2)+(C5/F3)
!       WRITE(*,'(1X,A,I5,E15.6)') 'K=1,NE.0,ZEROJ',K,ZEROJ
        CALL       SECANT(N,NITMX,TOL,ZEROJ,IERROR)
        IER(1)=IERROR
        JZERO(1)=ZEROJ

      endif
        T0 = 4.D0*FN*FN
        T1 = T0-1.D0
        T3 = 4.D0*T1*(7.D0*T0-31.D0)
        T5 = 32.D0*T1*((83.D0*T0-982.D0)*T0+3779.D0)
        T7 = 64.D0*T1*(((6949.D0*T0-153855.D0)*T0+1585743.D0)*T0  &
                      -6277237.D0)

!    OTHER ZEROES
      DO K = 2,NK
        JZERO(K) = 0.D0
        FK = FLOAT(K)

!    MAC MAHON'S SERIES FOR K>>N

        B0 = (FK+.5D0*FN-.25D0)*PI
        B1 = 8.D0*B0
        B2 = B1*B1
        B3 = 3.D0*B1*B2
        B5 = 5.D0*B3*B2
        B7 = 7.D0*B5*B2
        ZEROJ = B0-(T1/B1)-(T3/B3)-(T5/B5)-(T7/B7)

        ERRJ=ABS(BESSJ(N,ZEROJ))
!      WRITE(*,'(1X,A,2I5,2E15.6)') 'N,K,ZEROJ,ERRJ',N,K,ZEROJ,ERRJ

!    IMPROVE SOLUTION USING SUBROUTINE SECANT

        IF (ERRJ>TOL) CALL SECANT(N,NITMX,TOL,ZEROJ,IERROR)
        JZERO(K)=ZEROJ
        IER(K)=IERROR
      enddo
      RETURN
      END SUBROUTINE ROOTJ
! ----------------------------------------------------------------------
      SUBROUTINE SECANT(N,NITMX,TOL,ZEROJ,IER)
      INTEGER N,NITMX,IER,IT,NEV,NTRY
      REAL(KIND=8) TOL,ZEROJ,P0,P1,Q0,Q1,DP,P,C(2)
      DATA C/.95D0,.999D0/
      NTRY=1
      IER=0

    5 P0 = C(NTRY)*ZEROJ

      P1 = ZEROJ
      NEV = 2
      Q0 = BESSJ(N,P0)
      Q1 = BESSJ(N,P1)
!     WRITE(*,'(1X,A,I5,4E15.6)') 'NTRY,P0,Q0,P1,Q1',NTRY,P0,Q0,P1,Q1
      DO IT = 1,NITMX
        IF(Q1==Q0) GO TO 15
        P = P1-Q1*(P1-P0)/(Q1-Q0)
        DP = P-P1
!       WRITE(*,'(1X,A,I5,4E15.6)') 'IT,P,DP',IT,P,DP
        IF (IT==1) GO TO 10
        IF (ABS(DP)<TOL) GO TO 20

   10   NEV = NEV+1
        P0 = P1
        Q0 = Q1
        P1 = P
        Q1 = BESSJ(N,P1)
      enddo
   15 NTRY=NTRY+1
      IF(NTRY<=2) GO TO 5
      IER= NTRY
   20 ZEROJ = P
      END SUBROUTINE SECANT

! ----------------------------------------------------------------------
      FUNCTION BESSJ (N,X)

!     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
!     OF ORDER N, INTEGER FOR ANY REAL X. WE USE HERE THE CLASSICAL
!     RECURRENT FORMULA, WHEN  X > N. FOR X < N, THE MILLER'S ALGORITHM
!     IS USED TO AVOID OVERFLOWS.
!     REFERENCE :
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.

      INTEGER N,IACC,J,JSUM,M
      REAL(KIND=8) BIGNO,BIGNI
      PARAMETER (IACC = 40,BIGNO = 1.D10, BIGNI = 1.D-10)
      !REAL(KIND=8) X,BESSJ,BESSJ0,BESSJ1,TOX,BJM,BJ,BJP,SUM
      REAL(KIND=8) X,BESSJ,TOX,BJM,BJ,BJP,SUM
      IF (N==0) THEN
      BESSJ = BESSJ0(X)
      RETURN
      endif
      IF (N==1) THEN
      BESSJ = BESSJ1(X)
      RETURN
      endif
      IF (X==0.) THEN
      BESSJ = 0.
      RETURN
      endif
      TOX = 2./X
      IF (X>FLOAT(N)) THEN
      BJM = BESSJ0(X)
      BJ  = BESSJ1(X)
      DO 11 J = 1,N-1
      BJP = J*TOX*BJ-BJM
      BJM = BJ
      BJ  = BJP
   11 CONTINUE
      BESSJ = BJ
      ELSE
      M = 2*((N+INT(SQRT(FLOAT(IACC*N))))/2)
      BESSJ = 0.
      JSUM = 0
      SUM = 0.
      BJP = 0.
      BJ  = 1.
      DO 12 J = M,1,-1
      BJM = J*TOX*BJ-BJP
      BJP = BJ
      BJ  = BJM
      IF (ABS(BJ)>BIGNO) THEN
      BJ  = BJ*BIGNI
      BJP = BJP*BIGNI
      BESSJ = BESSJ*BIGNI
      SUM = SUM*BIGNI
      endif
      IF (JSUM/=0) SUM = SUM+BJ
      JSUM = 1-JSUM
      IF (J==N) BESSJ = BJP
   12 CONTINUE
      SUM = 2.*SUM-BJ
      BESSJ = BESSJ/SUM
      endif
      RETURN
      END FUNCTION BESSJ
! ----------------------------------------------------------------------
      FUNCTION BESSJ0 (X)

!     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
!     OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
!     BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
!     REFERENCES :
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL(KIND=8) X,BESSJ0,AX,FR,FS,Z,FP,FQ,XX
      REAL(KIND=8) Y,P1,P2,P3,P4,P5,R1,R2,R3,R4,R5,R6,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,-.1098628627D-2,.2734510407D-4,  &
      -.2073370639D-5,.2093887211D-6 /
      DATA Q1,Q2,Q3,Q4,Q5 /-.1562499995D-1,.1430488765D-3,   &
      -.6911147651D-5,.7621095161D-6,-.9349451520D-7 /
      DATA R1,R2,R3,R4,R5,R6 /57568490574.D0,-13362590354.D0,    &
      651619640.7D0,-11214424.18D0,77392.33017D0,-184.9052456D0 /
      DATA S1,S2,S3,S4,S5,S6 /57568490411.D0,1029532985.D0,  &
      9494680.718D0,59272.64853D0,267.8532712D0,1.D0 /
      IF(X==0.D0) GO TO 1
      AX = ABS (X)
      IF (AX<8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ0 = FR/FS
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-.785398164
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ0 = SQRT(.636619772/AX)*(FP*COS(XX)-Z*FQ*SIN(XX))
      endif
      RETURN
    1 BESSJ0 = 1.D0
      RETURN
      END FUNCTION BESSJ0
! ----------------------------------------------------------------------
      FUNCTION BESSJ1 (X)

!     THIS FUNCTION RETURNS THE VALUE OF THE FIRST KIND BESSEL FUNCTION
!     OF ORDER 0 FOR ANY REAL X. WE USE HERE THE POLYNOMIAL APPROXIMATION
!     BY SERIES OF CHEBYSHEV POLYNOMIALS FOR 0<X<8 AND 0<8/X<1.
!     REFERENCES :
!     M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
!     C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
!     VOL.5, 1962.

      REAL(KIND=8) X,BESSJ1,AX,FR,FS,Z,FP,FQ,XX
      REAL(KIND=8) Y,P1,P2,P3,P4,P5,P6,R1,R2,R3,R4,R5,R6,Q1,Q2,Q3,Q4,Q5,S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5 /1.D0,.183105D-2,-.3516396496D-4,  &
      .2457520174D-5,-.240337019D-6 /,P6 /.636619772D0 /
      DATA Q1,Q2,Q3,Q4,Q5 /.04687499995D0,-.2002690873D-3,   &
      .8449199096D-5,-.88228987D-6,.105787412D-6 /
      DATA R1,R2,R3,R4,R5,R6 /72362614232.D0,-7895059235.D0, &
      242396853.1D0,-2972611.439D0,15704.48260D0,-30.16036606D0 /
      DATA S1,S2,S3,S4,S5,S6 /144725228442.D0,2300535178.D0, &
      18583304.74D0,99447.43394D0,376.9991397D0,1.D0 /

      AX = ABS(X)
      IF (AX<8.) THEN
      Y = X*X
      FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))))
      FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))))
      BESSJ1 = X*(FR/FS)
      ELSE
      Z = 8./AX
      Y = Z*Z
      XX = AX-2.35619491
      FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)))
      FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))
      BESSJ1 = SQRT(P6/AX)*(COS(XX)*FP-Z*SIN(XX)*FQ)*SIGN(S6,X)
      endif
      RETURN
      END FUNCTION BESSJ1

! ----------------------------------------------------------------------
!     Calculate modified Bessel function of order N=0 using polynomial approximation.
!     Handbook of Mathematical Functions by Milton Abramowitz and Irene Stegun, page 378.
!     F90 Release 1.2 By J-P Moreau, Paris.
      FUNCTION BESSI0(X)
      IMPLICIT NONE
      REAL(KIND=8) X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,  &
      Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,1.2067492D0,  &
      0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, &
      0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,  &
      0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X)<3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      endif
      RETURN
      END FUNCTION BESSI0

end module trootj

