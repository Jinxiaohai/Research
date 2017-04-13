        SUBROUTINE JETINI(JP,JT,itrig)
        CHARACTER BEAM*16,TARG*16
        DIMENSION XSEC0(8,0:200),COEF0(8,200,20),INI(8),
     &                MINT44(8),MINT45(8)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
        COMMON/HPINT/MINT4,MINT5,ATCO(200,20),ATXS(0:200)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)    
        COMMON/PYSUBS/MSEL,MSUB(200),KFIN(2,-40:40),CKIN(200)
        COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
        COMMON/PYINT1/MINT(400),VINT(400)
        COMMON/PYINT2/ISET(200),KFPR(200,2),COEF(200,20),ICOL(40,4,2)
        COMMON/PYINT5/NGEN(0:200,3),XSEC(0:200,3)
        SAVE
        DATA INI/8*0/,ilast/-1/
        IHNT2(11)=JP
        IHNT2(12)=JT
        IF(IHNT2(5).NE.0 .AND. IHNT2(6).NE.0) THEN
           itype=1
        ELSE IF(IHNT2(5).NE.0 .AND. IHNT2(6).EQ.0) THEN
           itype=1
           IF(NFT(JT,4).EQ.2112) itype=2
        ELSE IF(IHNT2(5).EQ.0 .AND. IHNT2(6).NE.0) THEN
           itype=1
           IF(NFP(JP,4).EQ.2112) itype=2
        ELSE
           IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2212) THEN
              itype=1
           ELSE IF(NFP(JP,4).EQ.2212 .AND. NFT(JT,4).EQ.2112) THEN
              itype=2
           ELSE IF(NFP(JP,4).EQ.2112 .AND. NFT(JT,4).EQ.2212) THEN
              itype=3
           ELSE
              itype=4
           ENDIF
        ENDIF
        IF(itrig.NE.0) GO TO 160
        IF(itrig.EQ.ilast) GO TO 150
        MSTP(2)=2
        MSTP(33)=1
        PARP(31)=HIPR1(17)
        MSTP(51)=3
        MSTP(61)=1
        MSTP(71)=1
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.2) MSTP(61)=0
        IF(IHPR2(2).EQ.0.OR.IHPR2(2).EQ.1) MSTP(71)=0
        MSTP(81)=0
        MSTP(82)=1
        MSTP(111)=0
        IF(IHPR2(10).EQ.0) MSTP(122)=0
        PARP(81)=HIPR1(8)
        CKIN(5)=HIPR1(8)
        CKIN(3)=HIPR1(8)
        CKIN(4)=HIPR1(9)
        IF(HIPR1(9).LE.HIPR1(8)) CKIN(4)=-1.0
        CKIN(9)=-10.0
        CKIN(10)=10.0
        MSEL=0
        DO 100 ISUB=1,200
           MSUB(ISUB)=0
 100    CONTINUE
        MSUB(11)=1
        MSUB(12)=1
        MSUB(13)=1
        MSUB(28)=1
        MSUB(53)=1
        MSUB(68)=1
        MSUB(81)=1
        MSUB(82)=1
        DO 110 J=1,MIN(8,MDCY(21,3))
 110    MDME(MDCY(21,2)+J-1,1)=0
        ISEL=4
        IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
        MDME(MDCY(21,2)+ISEL-1,1)=1
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
 150    IF(INI(itype).NE.0) GO TO 800
        GO TO 400
 160    itype=4+itype
        IF(itrig.EQ.ilast) GO TO 260
        PARP(81)=ABS(HIPR1(10))-0.25
        CKIN(5)=ABS(HIPR1(10))-0.25
        CKIN(3)=ABS(HIPR1(10))-0.25
        CKIN(4)=ABS(HIPR1(10))+0.25
        IF(HIPR1(10).LT.HIPR1(8)) CKIN(4)=-1.0
        MSEL=0
        DO 101 ISUB=1,200
           MSUB(ISUB)=0
 101    CONTINUE
        IF(IHPR2(3).EQ.1) THEN
           MSUB(11)=1
           MSUB(12)=1
           MSUB(13)=1
           MSUB(28)=1
           MSUB(53)=1
           MSUB(68)=1
           MSUB(81)=1
           MSUB(82)=1
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
           DO 102 J=1,MIN(8,MDCY(21,3))
 102           MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
        ELSE IF(IHPR2(3).EQ.2) THEN
           MSUB(14)=1
           MSUB(18)=1
           MSUB(29)=1
        ELSE IF(IHPR2(3).EQ.3) THEN
           CKIN(3)=MAX(0.0,HIPR1(10))
           CKIN(5)=HIPR1(8)
           PARP(81)=HIPR1(8)
           MSUB(81)=1
           MSUB(82)=1
           DO 105 J=1,MIN(8,MDCY(21,3))
 105           MDME(MDCY(21,2)+J-1,1)=0
           ISEL=4
           IF(HINT1(1).GE.20.0 .and. IHPR2(18).EQ.1) ISEL=5
           MDME(MDCY(21,2)+ISEL-1,1)=1
        ENDIF
260        IF(INI(itype).NE.0) GO TO 800
400        INI(itype)=1
        IF(IHPR2(10).EQ.0) MSTP(122)=0
        IF(NFP(JP,4).EQ.2212) THEN
                BEAM='P'
        ELSE IF(NFP(JP,4).EQ.-2212) THEN
                BEAM='P~'
        ELSE IF(NFP(JP,4).EQ.2112) THEN
                BEAM='N'
        ELSE IF(NFP(JP,4).EQ.-2112) THEN
                BEAM='N~'
        ELSE IF(NFP(JP,4).EQ.211) THEN
                BEAM='PI+'
        ELSE IF(NFP(JP,4).EQ.-211) THEN
                BEAM='PI-'
        ELSE IF(NFP(JP,4).EQ.321) THEN
                BEAM='PI+'
        ELSE IF(NFP(JP,4).EQ.-321) THEN
                BEAM='PI-'
        ELSE
                WRITE(6,*) 'unavailable beam type', NFP(JP,4)
        ENDIF
        IF(NFT(JT,4).EQ.2212) THEN
                TARG='P'
        ELSE IF(NFT(JT,4).EQ.-2212) THEN
                TARG='P~'
        ELSE IF(NFT(JT,4).EQ.2112) THEN
                TARG='N'
        ELSE IF(NFT(JT,4).EQ.-2112) THEN
                TARG='N~'
        ELSE IF(NFT(JT,4).EQ.211) THEN
                TARG='PI+'
        ELSE IF(NFT(JT,4).EQ.-211) THEN
                TARG='PI-'
        ELSE IF(NFT(JT,4).EQ.321) THEN
                TARG='PI+'
        ELSE IF(NFT(JT,4).EQ.-321) THEN
                TARG='PI-'
        ELSE
                WRITE(6,*) 'unavailable target type', NFT(JT,4)
        ENDIF
        IHNT2(16)=1
        CALL PYINIT('CMS',BEAM,TARG,HINT1(1))
        MINT4=MINT(44)
        MINT5=MINT(45)
        MINT44(itype)=MINT(44)
        MINT45(itype)=MINT(45)
        ATXS(0)=XSEC(0,1)
        XSEC0(itype,0)=XSEC(0,1)
        DO 500 I=1,200
                ATXS(I)=XSEC(I,1)
                XSEC0(itype,I)=XSEC(I,1)
                DO 500 J=1,20
                        ATCO(I,J)=COEF(I,J)
                        COEF0(itype,I,J)=COEF(I,J)
500        CONTINUE
        IHNT2(16)=0
        RETURN
800        MINT(44)=MINT44(itype)
        MINT(45)=MINT45(itype)
        MINT4=MINT(44)
        MINT5=MINT(45)
        XSEC(0,1)=XSEC0(itype,0)
        ATXS(0)=XSEC(0,1)
        DO 900 I=1,200
                XSEC(I,1)=XSEC0(itype,I)
                ATXS(I)=XSEC(I,1)
        DO 900 J=1,20
                COEF(I,J)=COEF0(itype,I,J)
                ATCO(I,J)=COEF(I,J)
900        CONTINUE
        ilast=itrig
        MINT(11)=NFP(JP,4)
        MINT(12)=NFT(JT,4)
        RETURN
        END
