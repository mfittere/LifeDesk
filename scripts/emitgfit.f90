PROGRAM EMITGFIT
IMPLICIT NONE
INTEGER        IARGC
INTEGER        hbuf_4(30031), cnts(2), nstep, flag, leng, IC, i, j, k, m
INTEGER(1)     hbuf_1(120124)
EQUIVALENCE   (hbuf_4, hbuf_1)
REAL(8)        hist(1001,5,3), emitt(2), tfact, val
LOGICAL        equilibr
CHARACTER(160) fname
CHARACTER(80)  str

IF (IARGC() == 0) STOP "Filename (LIFETRAC's output file) was not specified"
IF (IARGC() /= 1) STOP "Must be only one command line argument"
CALL GETARG (1, fname)
i = LEN_TRIM (fname)
IF (fname(i-3:i) /= '.ltr') fname(i+1:) = '.ltr'
OPEN (1, FILE=fname, STATUS='old', IOSTAT=IC)
IF (IC /= 0) THEN
  PRINT *, 'Error opening file: '//TRIM(fname)
  STOP
END IF
flag = 0; hist = 0; nstep = 0; emitt = 0; tfact = 1

DO
  READ (1, '(A)', IOSTAT=IC) str
  IF (IC /= 0.AND.flag < 2) EXIT
  IF (flag == 2.AND.(str == ' '.OR.INDEX(str,':') > 0)) THEN
    flag = 0
    CALL UNZIP_7 (hbuf_4, i, leng, k)
    IF (k == 0) CALL UNZIP_D (hbuf_4, leng, i, j, k)
    IF (k == 0.AND.j == 4) CALL UNZIP_0 (hbuf_4, i, leng, k)
    IF (k /= 0.OR.j /= 4.OR.i /= 30031.OR.hbuf_4(i) /= hbuf_4(1)) THEN
      PRINT *, 'File: '//TRIM(fname)//' Error reading histograms, step:', nstep
      STOP
    ELSE IF (equilibr) THEN
      hist = hist + tfact * RESHAPE (hbuf_4(:15015), (/1001,5,3/))
    ELSE
      hist = RESHAPE (hbuf_4(:15015), (/1001,5,3/))
      CALL GFIT
    END IF
    IF (IC /= 0) EXIT
  END IF
  IF (flag < 2) CALL UPCASE (str)
  IF (flag == 1.AND.str == 'HISTOGRAMS:') THEN
    flag = 2; leng = 0
  ELSE IF (flag == 2) THEN
    i = LEN_TRIM (str)
    hbuf_1(leng+1:leng+i) = TRANSFER (str(:i), 0_1, i)
    leng = leng + i
  ELSE IF (nstep == 0.AND.str(:10) == 'BOUNDARY:') THEN
    equilibr = (INDEX(str,'(PART)') == 0)
  ELSE IF (str(:5) == 'STEP_'.AND.str(10:) == ':') THEN
    i = VERIFY (str(:10), '_:', BACK=.TRUE.)
    READ (str(6:i), *) nstep
    flag = 1
  ELSE IF (flag == 1.AND.equilibr.AND.str(:12) == 'TURN_COUNT:') THEN
    IF (nstep == 1) THEN
      i = INDEX (str,  '(INIT)=')
      j = INDEX (str(i+1:), ' ')
      READ (str(i+7:i+j), *) k
    ELSE
      tfact = tfact * cnts(2) / cnts(1)
    END IF
    i = INDEX (str,  '(ALL)=')
    j = INDEX (str(i+1:), ' ')
    READ (str(i+6:i+j), *) cnts(1)
    IF (nstep == 1) cnts(1) = cnts(1) - k
    i = INDEX (str,  '(OUT)=')
    j = INDEX (str(i+1:), ' ')
    READ (str(i+6:i+j), *) cnts(2)
  ELSE IF (nstep == 1.AND.str(:6) == '|EMIT|') THEN
    str = ADJUSTL (str(7:))
    DO i = 1, 2
      DO j = 1, 2
        m = INDEX (str, ' ')
        IF (str(1:1) == '(') THEN
          str = str(2:)
          m = m - 2
        END IF
        READ (str(:m-1), *, IOSTAT=IC) val
        IF (IC /= 0) THEN
          PRINT *, 'File: '//TRIM(fname)//'  Error reading emittances'
          STOP
        END IF
        str = ADJUSTL (str(m+1:))
        IF (j == 1) emitt(i) = val
        IF (j == 2) emitt(i) = emitt(i) / val
      END DO
    END DO
  END IF
END DO

CLOSE (1)
IF (equilibr.AND.ANY(hist /= 0)) THEN
  hist = hist + tfact * RESHAPE (hbuf_4(15016:30030), (/1001,5,3/))
  CALL GFIT
END IF
CONTAINS
!------------------------------------------------------------------------------

SUBROUTINE GFIT
INTEGER pass, i, j, k
REAL(8) data(1001), sigm(2), sig2, pos1, pos2, expp, dval, dold, dmax

DO i = 1, 2
  data = hist(:,5,i) / SUM (hist(:,5,i)) * 10
  DO k = 1000, 1, -1
    IF (data(k) /= 0) EXIT
  END DO
  sigm(i) = 1; dmax = 0.2_8
  DO pass = 1, 20
    dold = dval; dval = 0
    sig2 = 2 * sigm(i) ** 2
    DO j = 1, k
      pos1 = (j - 0.5_8) * 0.1_8
      pos2 = pos1 ** 2
      expp = pos1 * EXP (-pos2/sig2)
      dval = dval + expp * (sig2 - pos2) * (expp - sig2 * data(j) / 2)
    END DO
    IF (pass > 1.AND.dval*dold < 0) dmax = dmax / 2
    sigm(i) = sigm(i) + SIGN (MIN(dmax,ABS(dval)/20), dval)
    IF (pass > 3.AND.MIN(ABS(dval),ABS(dval-dold)) < 0.001) EXIT
  END DO
END DO
PRINT *, REAL(emitt*sigm**2), nstep
END SUBROUTINE GFIT
END
