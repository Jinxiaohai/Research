       REAL FUNCTION RMASS(DMAX,ISEED)
      COMMON/RNDF77/NSEED
      SAVE   
          DMIN = 1.078
          IF(DMAX.LT.1.232) THEN
          FM=FDELTA(DMAX)
          ELSE
          FM=1.
          ENDIF
          IF(FM.EQ.0.)FM=1.E-06
          NTRY1=0
10        DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
          NTRY1=NTRY1+1
          IF((RANART(NSEED) .GT. FDELTA(DM)/FM).AND.
     1    (NTRY1.LE.10)) GOTO 10
          if(dm.gt.1.47) goto 10
       RMASS=DM
       RETURN
       END
