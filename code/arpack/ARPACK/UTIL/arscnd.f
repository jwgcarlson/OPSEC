      SUBROUTINE SECOND( T )
*
      REAL       T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the CPU_TIME intrinsic.
*
*     .. Local Scalars ..
      REAL               T1
*     ..
*     .. Local Arrays ..
*     ..
*     .. External Functions ..
*     ..
*     .. Executable Statements ..
*
      CALL CPU_TIME(T1)
      T = T1

      RETURN
*
*     End of SECOND
*
      END
