        SUBROUTINE ATTFLV(ID,IDQ,IDQQ)
      COMMON/RNDF77/NSEED
        SAVE   
        IF(ABS(ID).LT.100) THEN
                NSIGN=1
                IDQ=ID/100
                IDQQ=-ID/10+IDQ*10
                IF(ABS(IDQ).EQ.3) NSIGN=-1
                IDQ=NSIGN*IDQ
                IDQQ=NSIGN*IDQQ
                IF(IDQ.LT.0) THEN
                        ID0=IDQ
                        IDQ=IDQQ
                        IDQQ=ID0
                ENDIF
                RETURN
        ENDIF
        IDQ=2
        IF(ABS(ID).EQ.2112) IDQ=1
        IDQQ=2101
        X=RANART(NSEED)
        IF(X.LE.0.5) GO TO 30
        IF(X.GT.0.666667) GO TO 10
        IDQQ=2103
        GO TO 30
10        IDQ=1
        IDQQ=2203
        IF(ABS(ID).EQ.2112) THEN
                IDQ=2
                IDQQ=1103
        ENDIF
30        IF(ID.LT.0) THEN
                ID00=IDQQ
                IDQQ=-IDQ
                IDQ=-ID00
        ENDIF
        RETURN
        END        
