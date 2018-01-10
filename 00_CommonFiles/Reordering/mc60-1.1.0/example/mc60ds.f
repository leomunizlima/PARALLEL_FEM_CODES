      INTEGER LIRN,MAXN
      PARAMETER (LIRN=20, MAXN=5)
      INTEGER CASE,N,NNZ,IRN(LIRN),ICPTR(MAXN+1),ICNTL(2),
     *      IW(3*MAXN+1),INFO(4),NSUP,I,SVAR(MAXN),VARS(MAXN),
     *      PERM(MAXN),PERMSV(MAXN),JCNTL(2),PAIR(2,MAXN/2)
      DOUBLE PRECISION WEIGHT(2),W(MAXN),RINFO(4)
      CHARACTER ALG

C Set parameter values
      ICNTL(1) = 0
      ICNTL(2) = 6
      JCNTL(1) = 0
      JCNTL(2) = 0
      WEIGHT(1) = 2.0
      WEIGHT(2) = 1.0

      DO 20 CASE = 1,2
C Read in data for the lower-triangular part
        READ (5,*) N,NNZ
        READ (5,*) (IRN(I),I = 1,NNZ)
        READ (5,*) (ICPTR(I),I = 1,N+1)

C Construct pattern of whole matrix
        CALL MC60AD(N,LIRN,IRN,ICPTR,ICNTL,IW,INFO)

C Check for an error return
        IF (INFO(1).NE.0) THEN
          WRITE(6,'(A,2I3)')' MC60AD failed with INFO=',INFO
          STOP
        END IF

        READ(5,*) ALG
        IF (ALG.EQ.'S') THEN
C Work with supervariables
          CALL MC60BD(N,LIRN,IRN,ICPTR,NSUP,SVAR,VARS,IW)
          WRITE(6,'(A,5I4)') 'The number of supervariables is', NSUP
          CALL MC60CD(N,NSUP,LIRN,IRN,ICPTR,VARS,JCNTL,PERMSV,WEIGHT,
     *            PAIR,INFO,IW,W)
          CALL MC60FD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERMSV,IW,RINFO)
          CALL MC60DD(N,NSUP,SVAR,VARS,PERMSV,PERM,IW)
        ELSE
C Work with variables
          NSUP = N
          DO 10 I = 1,N
            VARS(I) = 1
   10     CONTINUE
          CALL MC60CD(N,NSUP,LIRN,IRN,ICPTR,VARS,JCNTL,PERM,WEIGHT,
     *            PAIR,INFO,IW,W)
          CALL MC60FD(N,NSUP,LIRN,IRN,ICPTR,VARS,PERM,IW,RINFO)
        END IF

        WRITE(6,'(A,5I4)') 'The chosen permutation is', PERM
        WRITE(6,'(A,F4.0,/)') 'The profile is', RINFO(1)

   20 CONTINUE

       END
