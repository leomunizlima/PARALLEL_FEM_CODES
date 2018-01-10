* COPYRIGHT (c) 1990 AEA Technology and
* Council for the Central Laboratory of the Research Councils
C Original date Aug 2001
C  August 2001: threadsafe version of MC45
C
C 12th July 2004 Version 1.0.0. Version numbering added.
C 4th September 2007 Version 1.1.0 Sizes of control and info arrays
C                    for MA48 calls updated

      SUBROUTINE MC75I(ICNTL)
      INTEGER ICNTL(5)

C  initialize optional user controls
      ICNTL(1) = 6
      ICNTL(2) = 0
      ICNTL(3) = 0
      ICNTL(4) = 0
      ICNTL(5) = 0

      RETURN
      END
      SUBROUTINE MC75A(N,NZ,LA,A,IRN,JCN,COND,LIW,IW,LW,W,ICNTL,INFO)
C
C THIS SUBROUTINE COMPUTES THE CLASSICAL CONDITION NUMBER
C AND THE SKEEL'S CONDITION NUMBER, IN INFINITY NORM.
C
C THE PARAMETERS ARE...
C
C N     IS EQUAL TO THE ORDER OF THE MATRIX. IT IS NOT ALTERED BY THE
C       SUBROUTINE.
C NZ    IS THE NUMBER OF ENTRIES IN A. IT IS NOT ALTERED BY THE
C       SUBROUTINE.
C LA    INTEGER SPECIFYING THE LENGTH OF A, IRN, JCN..  IT IS CLOSELY
C       TIED TO WORKSPACE REQUIREMENTS OF MA48 AND 4*NZ IS SUGGESTED
C       AS GOOD FIRST CHOICE.
C A     IS A  REAL ARRAY  LENGTH LA.  HOLDS NON-ZEROS OF MATRIX
C       ON ENTRY AND NON-ZEROS OF FACTORS ON EXIT.  REORDERED BY
C       MC20A/AD AND MC23A/AD AND ALTERED BY M280A/AD.
C IRN   INTEGER ARRAY OF LENGTH LA THAT HOLDS ROW INDICES ON INPUT.
C       USED AS WORKSPACE BY MA48A/AD.
C JCN   INTEGER ARRAY OF LENGTH LA. IT HOLDS COLUMN INDICES ON ENTRY.
C       IT HOLDS THE COLUMN INDICES OF DECOMPOSED MATRIX ON EXIT.
C       REORDERED BY MC20A/AD AND MC23A/AD AND ALTERED BY MA48A/AD.
C COND  IS A REAL ARRAY OF LENGTH 2.  COND(1) IS THE CLASSICAL CONDITION
C       NUMBER, COND(2) IS THE SKEEL'S CONDITION NUMBER. IF THE MATRIX
C       IS SINGULAR COND(2) IS SET TO ZERO. IF IFLAG IS LESS THAN ZERO
C       COND(1) AND COND(2) ARE SET TO ZERO.
C LIW   INTEGER SPECIFYING THE LENGTH OF IW.  IT IS CLOSELY TIED TO
C       WORKSPACE REQUIREMENTS OF MA48 AND MUST BE AT LEAST 19*N+7.
C IW    IS AN INTEGER ARRAY OF LENGTH LIW. IT IS USED AS WORKSPACE
C LW    INTEGER SPECIFYING THE LENGTH OF W.  IT IS CLOSELY TIED TO
C       WORKSPACE REQUIREMENTS OF MA48 AND MUST BE AT LEAST 7*N.
C W     IS A REAL ARRAY OF LENGTH LW. IT IS USED AS WORKSPACE.
C ICNTL IS AN INTEGER ARRAY, LENGTH 5, FOR USER OPTIONS
C   ICNTL(1) INTEGER  DEFAULT VALUE 6 (LINE PRINTER).  UNIT NUMBER
C     FOR ERROR MESSAGES.
C   ICNTL(2) INTEGER  DEFAULT VALUE 6 (LINE PRINTER).  UNIT NUMBER
C     FOR INFORMATION MESSAGES.
C   ICNTL(3) TO ICNTL(5) NOT USED.
C INFO  IS AN INTEGER ARRAY, LENGTH 5, RETURNS INFORMATION.
C   INFO(1) SET ON RETURN TO A POSITIVE OR ZERO VALUE TO INDICATE
C       SUCCESS, OR NEGATIVE TO INDICATE AN ERROR.
C   INFO(2) SUPPLEMENTARY INFORMATION TO INFO(1).
C   INFO(3) INTEGER.  MINIMUM LENGTH OF ARRAYS IRN, JCN AND A, I.E.
C           VALUE OF LA, FOR SUCCESS ON FUTURE RUNS.
C   INFO(4) INTEGER   ESTIMATED RANK OF MATRIX.
C   INFO(5) NOT USED
C
C
C     .. Parameters ..
      REAL ZERO,U
      PARAMETER (ZERO=0.0,U=0.1)
C     ..
C     .. Scalar Arguments ..
      INTEGER LA,N,NZ,LIW,LW
C     ..
C     .. Array Arguments ..
      REAL A(LA),COND(2),W(LW)
      INTEGER IRN(LA),IW(LIW),JCN(LA),ICNTL(5),INFO(5)
C     ..
C     .. Local Arrays ..
      REAL ERROR(3)
      INTEGER KEEP71(5)
C      MA48 control and info arrays
      REAL CNTL48(10),RINF48(10)
      INTEGER ICNT48(20),INFO48(20)
C     ..
C     .. Local Scalars ..
      REAL DANRM
      INTEGER I,IDUMMY,K,KASE,IW71
      INTEGER IIW,LIIW,IKEEP,LIKEEP,MINIW,IW75B,IWRHS,IWX,IW48,MINW
C     ..
C     .. External Functions ..
      INTEGER ISAMAX
      EXTERNAL ISAMAX
C     ..
C     .. External Subroutines ..
      EXTERNAL MA48A,MA48B,MA48C,MC71A,MC75B
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Executable Statements ..
C
C
C   CHECK DATA AND INITIALIZATIONS
C
      INFO(1) = 0
      IF (N.LE.0) THEN
        INFO(1) = -7
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9700) INFO(1)
 9700   FORMAT (1X,'Error in MC75: no.',I5,' (N LE 0)')
        GO TO 1010
      END IF

      IF (NZ.LE.0) THEN
        INFO(1) = -6
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9600) INFO(1)
 9600   FORMAT (1X,'Error in MC75: no.',I5,' (NZ LE 0)')
        GO TO 1010
      END IF

      IF (LA.LT.NZ) THEN
        INFO(1) = -5
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9500) INFO(1)
 9500   FORMAT (1X,'Error in MC75: no.',I5,' (LA LT NZ)')
        GO TO 1010
      END IF

      COND(1) = ZERO
      COND(2) = ZERO

C
C  PARTITION INTEGER WORKSPACE IW
C
      IIW = 0
      LIIW = 9*N
      IKEEP = IIW + LIIW
      LIKEEP = 10*N + 7
      MINIW = IKEEP + LIKEEP

      IF (MINIW.GT.LIW) THEN
        INFO(1) = -4
        INFO(2) = MINIW
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9400) INFO(1)
 9400   FORMAT (1X,'Error in MC75: no.',I5,' (IW array too short)')
        GO TO 1010
      END IF
C
C  PARTITION REAL WORKSPACE W
C
      IW75B = 0
      IW71  = IW75B + N
      IWRHS = IW71 + N
      IWX   = IWRHS + N
      IW48  = IWX + N
      MINW  = IW48 + 4*N

      IF (MINW.GT.LW) THEN
        INFO(1) = -3
        INFO(2) = MINW
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9300) INFO(1)
 9300   FORMAT (1X,'Error in MC75: no.',I5,' (W array too short)')
        GO TO 1010
      END IF
C
C  COMPUTE THE VECTOR W(I,1) = SUM (|A(I,J)|) FROM J=1 TO N
C
      DO 140 I = 1,N
        W(IW75B+I) = ZERO
  140 CONTINUE
      DO 150 K = 1,NZ
        I = IRN(K)
        W(IW75B+I) = W(IW75B+I) + ABS(A(K))
  150 CONTINUE
      DANRM = ABS(W(IW75B+ISAMAX(N,W(IW75B+1),1)))
C
C  INITIALIZE MA48
C
      CALL MA48I(CNTL48,ICNT48)
      ICNT48(1) = ICNTL(2)
      ICNT48(2) = ICNTL(2)
      ICNT48(7) = 0
C
C  ANALYSE THE LU DECOMPOSITION
C
      CALL MA48A(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48,
     +           ICNT48,IW(IIW+1),INFO48,RINF48)
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
C
      IF (INFO48(1).EQ.1) THEN
C        nonzero's index out of range or duplicate entries
        IF (INFO48(12).GT.0) THEN
C          there are indices out of range
          INFO(1) = -8
        ELSE
C          all indices in range therefore must be duplicate entries
          INFO(1) = 1
        ENDIF
      ELSE IF (INFO48(1).EQ.-3) THEN
C        LA too small, minimum value in INFO(3), recommended in INFO(4)
        INFO(1) = -1
      ELSE IF (INFO48(1).LT.0) THEN
C        unexpected MA48A/AD errors
        INFO(1) = -9
        INFO(2) = INFO48(1)
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9009) INFO(1)
 9009   FORMAT (1X,'Error in MC75: no.',I5,' (MA48A/AD error)')
        GO TO 1010
      END IF

      IF (INFO(1).LT.0) GO TO 9999
C
C  COMPUTE THE LU DECOMPOSITION
C
      CALL MA48B(N,N,NZ,1,LA,A,IRN,JCN,IW(IKEEP+1),CNTL48,
     +           ICNT48,W(IW48+1),IW(IIW+1),INFO48,RINF48)
C
      INFO(3) = INFO48(3)
      INFO(4) = INFO48(4)
      INFO(5) = INFO48(5)
      IF (INFO48(1).EQ.-3) THEN
C        LA too small required value in INFO(4)
        INFO(1) = -2
      ELSE IF (INFO48(1).EQ.2) THEN
C        A is rank defficient, rank in INFO(5)
        INFO(1) = 2
      ELSE IF (INFO48(1).LT.0) THEN
C        unexpected MA48B/BD errors
        INFO(1) = -10
        INFO(2) = INFO48(1)
        IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9010) INFO(1)
 9010   FORMAT (1X,'Error in MC75: no.',I5,' (MA48B/BD error)')
        GO TO 1010
      END IF

      IF (INFO(1).LT.0) GO TO 9999
C
C  COMPUTE THE CONDITION NUMBERS
C
      KASE = 0
      IF (INFO(1).LE.1) THEN
        DO 170 IDUMMY = 1,100
          CALL MC71A(N,KASE,W(IWX+1),COND(2),W(IW71+1),IW(IIW+1),KEEP71)
          IF (KASE.LE.0) GO TO 160
          IF (KASE.EQ.1) CALL MC75B(N,W(IWX+1),W(IW75B+1))
          DO 165 I = 1,N
            W(IWRHS+I) = W(IWX+I)
  165     CONTINUE
          CALL MA48C(N,N,KASE.EQ.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48,
     +               ICNT48,W(IWRHS+1),W(IWX+1),ERROR,
     +               W(IW48+1),IW(IIW+1),INFO48)
          IF (INFO48(1).LT.0) THEN
C             unexpected MA48C/CD errors
             INFO(1) = -11
             INFO(2) = INFO48(1)
             IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9011) INFO(1)
 9011        FORMAT (1X,'Error in MC75: no.',I5,' (MA48C/CD error)')
             GO TO 1010
          END IF
          IF (KASE.EQ.2) CALL MC75B(N,W(IWX+1),W(IW75B+1))
  170   CONTINUE
      END IF

  160 KASE = 0
      DO 190 IDUMMY = 1,100
        CALL MC71A(N,KASE,W(IWX+1),COND(1),W(IW71+1),IW(IIW+1),KEEP71)
        IF (KASE.LE.0) GO TO 180
        DO 175 I = 1,N
          W(IWRHS+I) = W(IWX+I)
  175   CONTINUE
        CALL MA48C(N,N,KASE.EQ.1,1,LA,A,IRN,IW(IKEEP+1),CNTL48,
     +             ICNT48,W(IWRHS+1),W(IWX+1),ERROR,
     +             W(IW48+1),IW(IIW+1),INFO48)
        IF (INFO48(1).LT.0) THEN
C             unexpected MA48C/CD errors
           INFO(1) = -11
           INFO(2) = INFO48(1)
           IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9011) INFO(1)
           GO TO 1010
        END IF
  190 CONTINUE
  180 COND(1) = COND(1)*DANRM
      GO TO 1010

 9999 IF (ICNTL(1).GT.0) WRITE (ICNTL(1),FMT=9000) INFO(1)
 9000 FORMAT (1X,'Error in MC75: no.',I5)

 1010 RETURN
      END
C
      SUBROUTINE MC75B(N,R,W)
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL R(*),W(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Executable Statements ..
      DO 100 I = 1,N
        R(I) = R(I)*W(I)
  100 CONTINUE
      RETURN
      END
