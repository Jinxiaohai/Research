        SUBROUTINE CRSJET
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL HIPR1(100),HINT1(100)
        COMMON/HPARNT/HIPR1,IHPR2(50),HINT1,IHNT2(50)
        COMMON/NJET/N,ipcrs
        COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
        COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
        COMMON/BVEG3/F,TI,TSI
        COMMON/SEDVAX/NUM1
        EXTERNAL FJET,FJETRG
        SAVE   
        NDIM=3
        ipcrs=0
        CALL VEGAS(FJET,AVGI,SD,CHI2A)
        HINT1(14)=sngl(AVGI)/2.5682
        IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
                ipcrs=1
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(15)=sngl(AVGI)/2.5682
        ENDIF
        IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
                ipcrs=2
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(16)=sngl(AVGI)/2.5682
        ENDIF
        IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
                ipcrs=3
                CALL VEGAS(FJET,AVGI,SD,CHI2A)
                HINT1(17)=sngl(AVGI)/2.5682
        ENDIF
        IF(IHPR2(3).NE.0) THEN
           ipcrs=0
           CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
           HINT1(61)=sngl(AVGI)/2.5682
           IF(IHPR2(6).EQ.1 .AND. IHNT2(1).GT.1) THEN
              ipcrs=1
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(62)=sngl(AVGI)/2.5682
           ENDIF
           IF(IHPR2(6).EQ.1 .AND. IHNT2(3).GT.1) THEN
              ipcrs=2
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(63)=sngl(AVGI)/2.5682
           ENDIF
           IF(IHPR2(6).EQ.1.AND.IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
              ipcrs=3
              CALL VEGAS(FJETRG,AVGI,SD,CHI2A)
              HINT1(64)=sngl(AVGI)/2.5682
           ENDIF
        ENDIF
        RETURN
        END
