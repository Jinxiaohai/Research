      FUNCTION INVFLV(IART)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
      IF (IART .EQ. -6) THEN
         INVFLV = -1114
         RETURN
      END IF
      IF (IART .EQ. -7) THEN
         INVFLV = -2114
         RETURN
      END IF
      IF (IART .EQ. -8) THEN
         INVFLV = -2214
         RETURN
      END IF
      IF (IART .EQ. -9) THEN
         INVFLV = -2224
         RETURN
      END IF
      IF (IART .EQ. -1) THEN
         INVFLV = -2212
         RETURN
      END IF
      IF (IART .EQ. -2) THEN
         INVFLV = -2112
         RETURN
      END IF
      IF (IART .EQ. 0) THEN
         INVFLV = 221
         RETURN
      END IF
      IF (IART .EQ. 1) THEN
         INVFLV = 2212
         RETURN
      END IF
      IF (IART .EQ. 2) THEN
         INVFLV = 2112
         RETURN
      END IF
      IF (IART .EQ. 3) THEN
         INVFLV = -211
         RETURN
      END IF
      IF (IART .EQ. 4) THEN
         INVFLV = 111
         RETURN
      END IF
      IF (IART .EQ. 5) THEN
         INVFLV = 211
         RETURN
      END IF
      IF (IART .EQ. 6) THEN
         INVFLV = 1114
         RETURN
      END IF
      IF (IART .EQ. 7) THEN
         INVFLV = 2114
         RETURN
      END IF
      IF (IART .EQ. 8) THEN
         INVFLV = 2214
         RETURN
      END IF
      IF (IART .EQ. 9) THEN
         INVFLV = 2224
         RETURN
      END IF
      IF (IART .EQ. 14) THEN
         INVFLV = 3122
         RETURN
      END IF
      IF (IART .EQ. -14) THEN
         INVFLV = -3122
         RETURN
      END IF 
      IF (IART .EQ. 15) THEN
         INVFLV = 3112
         RETURN
      END IF
      IF (IART .EQ. -15) THEN
         INVFLV = -3112
         RETURN
      END IF 
      IF (IART .EQ. 16) THEN
         INVFLV = 3212
         RETURN
      END IF
      IF (IART .EQ. -16) THEN
         INVFLV = -3212
         RETURN
      END IF
      IF (IART .EQ. 17) THEN
         INVFLV = 3222
         RETURN
      END IF
      IF (IART .EQ. -17) THEN
         INVFLV = -3222
         RETURN
      END IF 
      IF (IART .EQ. 21) THEN
            INVFLV = -321
         RETURN
      END IF
      IF (IART .EQ. 23) THEN
            INVFLV = 321
         RETURN
      END IF
      IF (IART .EQ. 22) THEN
         INVFLV = 130
         RETURN
      ENDIF
      IF (IART .EQ. 24) THEN
         INVFLV = 310
         RETURN
      ENDIF
      IF (IART .EQ. 25) THEN
         INVFLV = -213
         RETURN
      END IF
      IF (IART .EQ. 26) THEN
         INVFLV = 113
         RETURN
      END IF
      IF (IART .EQ. 27) THEN
         INVFLV = 213
         RETURN
      END IF
      IF (IART .EQ. 28) THEN
         INVFLV = 223
         RETURN
      END IF
      IF (IART .EQ. 29) THEN
         INVFLV = 333
         RETURN
      END IF
      IF (IART .EQ. 30) THEN
         INVFLV = 323
         IF (RANART(NSEED).GT.0.5) INVFLV = 313
         RETURN
      END IF
      IF (IART .EQ. -30) THEN
         INVFLV = -323
         IF (RANART(NSEED).GT.0.5) INVFLV = -313
         RETURN
      END IF
      IF (IART .EQ. 31) THEN
         INVFLV = 331
         RETURN
      END IF
      IF (IART .EQ. 32) THEN
         INVFLV = 777
         RETURN
      END IF
      IF (IART .EQ. 40) THEN
         INVFLV = 3312
         RETURN
      END IF                   
      IF (IART .EQ. -40) THEN
         INVFLV = -3312
         RETURN
      END IF
      IF (IART .EQ. 41) THEN
         INVFLV = 3322
         RETURN
      END IF
      IF (IART .EQ. -41) THEN
         INVFLV = -3322
         RETURN
      END IF
      IF (IART .EQ. 45) THEN
         INVFLV = 3334
         RETURN
      END IF
      IF (IART .EQ. -45) THEN
         INVFLV = -3334
         RETURN
      END IF
      IF (IART .EQ. 44) THEN
         INVFLV = 6666
         RETURN
      END IF
      IF (IART .EQ. 42) THEN
         INVFLV = 42
         RETURN
      ELSEIF (IART .EQ. -42) THEN         
         INVFLV = -42
         RETURN
      END IF
      INVFLV = IART - 10000
      RETURN
      END
