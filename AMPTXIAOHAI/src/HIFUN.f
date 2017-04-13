        SUBROUTINE HIFUN(I,XMIN,XMAX,FHB)
        COMMON/HIJHB/RR(10,201),XX(10,201)
        EXTERNAL FHB
        SAVE   
        FNORM=GAUSS1(FHB,XMIN,XMAX,0.001)
        DO 100 J=1,201
                XX(I,J)=XMIN+(XMAX-XMIN)*(J-1)/200.0
                XDD=XX(I,J)
                RR(I,J)=GAUSS1(FHB,XMIN,XDD,0.001)/FNORM
100        CONTINUE
        RETURN
        END
