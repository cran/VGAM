      SUBROUTINE vdigami(D, X, P, GPLOG, GP1LOG, PSIP, PSIP1, PSIDP,
     *     PSIDP1, IFAULT, TMAX)
C
C     ALGORITHM AS 187  APPL. STATIST. (1982) VOL.31, NO.3
C
C     Computes derivatives of the incomplete gamma integral for positive
C     parameters, X, P, using a series expansion if P > X or X <= 1, and
C     a continued fraction expansion otherwise.
C
C     Calculation of D(4) in line 60 corrected 5 October 1993.
C
C     N.B. The user must input values of the incomplete gamma, digamma
C          and trigamma functions.  These can be obtained using AS 239
C          (or 32), AS 103 and AS 121 respectively.
C
C
C
C
C  20130214; adapted by T. W. Yee to handle DOUBLE PRECISION arguments.
C  And declarations of *all* variables.
C  And a wrapper function written to call this subroutine.
C  TMAX is now input.
C  Seems to work but more testing is required.
C
C  20141108; A, C, CP, CPP, DSP, DSPP, DFP, DFPP, F, S, TMAXP etc. now 
C  declared, by T. W. Yee.
C  ABS() changed to DABS() too.
C
C
      DOUBLE PRECISION X, P, GPLOG, GP1LOG, PSIP, PSIP1, PSIDP, PSIDP1
      DOUBLE PRECISION TMAX
      INTEGER          IFAULT
C
      DOUBLE PRECISION A, AN, B, C, CP, CPC, CPP, DSP, DSPP, DFP, DFPP
      DOUBLE PRECISION F, PM1, S, S0, XLOG, TERM, TMAXP
C
C
C
C
C
      INTEGER          I, I2
      DOUBLE PRECISION PN(6), D(6), DP(6), DPP(6), ZERO, ONE, TWO
C     DATA TMAX/100.0/
      DATA E, OFLO, VSMALL/1.D-6, 1.D30, 1.D-30/
      DATA ZERO/0.0/, ONE/1.0/, TWO/2.0/
C
      IFAULT = 0
C
C     Derivatives with respect to X
C
      PM1 = P - ONE
      XLOG = DLOG(X)
      D(1) = DEXP(-GPLOG + PM1*XLOG - X)
      D(2) = D(1) * (PM1/X - ONE)
      D(5) = D(1) * (XLOG - PSIP)
C
C     Derivatives with respect to P
C
      IF (X .GT. ONE .AND. X .GE. P) GO TO 30
C
C     Series expansion
C
      F = DEXP(P*XLOG - GP1LOG - X)
      DFP = F * (XLOG - PSIP1)
      DFPP = DFP*DFP/F - F*PSIDP1
C
      TMAXP = TMAX + P
      C = ONE
      S = ONE
      CP = ZERO
      CPP = ZERO
      DSP = ZERO
      DSPP = ZERO
      A = P
    1 A = A + ONE
      CPC = CP / C
      CP = CPC - ONE/A
      CPP = CPP/C - CPC*CPC + ONE/A**2
      C = C*X/A
      CP = CP*C
      CPP = CPP*C + CP*CP/C
      S = S + C
      DSP = DSP + CP
      DSPP = DSPP + CPP
      IF (A .GT. TMAXP) GO TO 1001
      IF (C .GT. E*S) GO TO 1
      D(6) = S*F
      D(3) = S*DFP + F*DSP
      D(4) = S*DFPP + TWO*DFP*DSP + F*DSPP
      RETURN
C
C     Continued fraction expansion
C
   30 F = DEXP(P*XLOG - GPLOG - X)
      DFP = F * (XLOG - PSIP)
      DFPP = DFP*DFP/F - F*PSIDP
C
      A = PM1
      B = X + ONE - A
      TERM = ZERO
      PN(1) = ONE
      PN(2) = X
      PN(3) = X + ONE
      PN(4) = X * B
      S0 = PN(3) / PN(4)
      DO 31 I = 1, 4
        DP(I) = ZERO
        DPP(I) = ZERO
   31 CONTINUE
      DP(4) = -X
C
   32 A = A - ONE
      B = B + TWO
      TERM = TERM + ONE
      AN = A*TERM
      PN(5) = B*PN(3) + AN*PN(1)
      PN(6) = B*PN(4) + AN*PN(2)
      DP(5) = B*DP(3) - PN(3) + AN*DP(1) + PN(1)*TERM
      DP(6) = B*DP(4) - PN(4) + AN*DP(2) + PN(2)*TERM
      DPP(5) = B*DPP(3) + AN*DPP(1) + TWO*(TERM*DP(1) - DP(3))
      DPP(6) = B*DPP(4) + AN*DPP(2) + TWO*(TERM*DP(2) - DP(4))
C
      IF (DABS(PN(6)) .LT. VSMALL) GO TO 35
      S = PN(5) / PN(6)
      C = DABS(S - S0)
      IF (C*P .GT. E) GO TO 34
      IF (C .LE. E*S) GO TO 42
C
   34 S0 = S
   35 DO 36 I = 1, 4
        I2 = I + 2
        DP(I) = DP(I2)
        DPP(I) = DPP(I2)
        PN(I) = PN(I2)
   36 CONTINUE
C
      IF (TERM .GT. TMAX) GO TO 1001
      IF (DABS(PN(5)) .LT. OFLO) GO TO 32
      DO 41 I = 1, 4
        DP(I) = DP(I) / OFLO
        DPP(I) = DPP(I) / OFLO
        PN(I) = PN(I) / OFLO
   41 CONTINUE
      GO TO 32
C
   42 D(6) = ONE - F*S
      DSP = (DP(5) - S*DP(6)) / PN(6)
      DSPP = (DPP(5) - S*DPP(6) - TWO*DSP*DP(6)) / PN(6)
      D(3) = -F*DSP - S*DFP
      D(4) = -F*DSPP - TWO*DSP*DFP - S*DFPP
      RETURN
C
C     Set fault indicator
C
 1001 IFAULT = 1
      RETURN
      END








