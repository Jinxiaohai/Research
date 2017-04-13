       REAL FUNCTION RHOMAS(DMAX,ISEED)
      COMMON/RNDF77/NSEED
      SAVE   
          DMIN = 0.28
          IF(DMAX.LT.0.77) THEN
          FM=FRHO(DMAX)
          ELSE
          FM=1.
          ENDIF
          IF(FM.EQ.0.)FM=1.E-06
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FRHO(DM)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 10
          if(dm.gt.1.07) goto 10
       RHOMAS=DM
       RETURN
       END
