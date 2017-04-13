      SUBROUTINE VEGAS(FXN,AVGI,SD,CHI2A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
      COMMON/BVEG3/F,TI,TSI   
      EXTERNAL FXN
      DIMENSION D(50,10),DI(50,10),XIN(50),R(50),DX(10),DT(10),X(10)
     1   ,KG(10),IA(10)
      REAL QRAN(10)
      SAVE   
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1.D0/,MDS/-1/
      NDO=1
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
      ENTRY VEGAS1(FXN,AVGI,SD,CHI2A)
      IT=0
      SI=0.d0
      SI2=SI
      SWGT=SI
      SCHI=SI
      ENTRY VEGAS2(FXN,AVGI,SD,CHI2A)
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=int((real(NCALL)/2.)**(1./real(NDIM)))
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE/CALLS
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.d0
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6 I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      NDO=ND
8     CONTINUE
      ENTRY VEGAS3(FXN,AVGI,SD,CHI2A)
9     IT=IT+1
      TI=0.d0
      TSI=TI
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      D(I,J)=TI
10    DI(I,J)=TI
11    FB=0.d0
      F2B=FB
      K=0
12    K=K+1
      CALL ARAN9(QRAN,NDIM)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=dble(float(KG(J))-QRAN(J))*DXG+ONE
      IA(J)=int(XN)
      IF(IA(J).GT.1) GO TO 13
      XO=XI(IA(J),J)
      RC=(XN-IA(J))*XO
      GO TO 14
13    XO=XI(IA(J),J)-XI(IA(J)-1,J)
      RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
14    X(J)=XL(J)+RC*DX(J)
      WGT=WGT*XO*XND
15    CONTINUE
      F=WGT
      F=F*FXN(X,WGT)
      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      DI(IA(J),J)=DI(IA(J),J)+F
16    IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
      IF(K.LT.NPG) GO TO 12
      F2B=DSQRT(F2B*NPG)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
17    D(IA(J),J)=D(IA(J),J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/(TSI+1.0d-37)
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SWGT=SWGT+1.0D-37
      SI2=SI2+1.0D-37
      SCHI=SCHI+TI2*WGT
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/dble(float(IT)-.999)
      SD=DSQRT(ONE/SD)
      IF(NPRN.EQ.0) GO TO 21
      TSI=DSQRT(TSI)
21    DO 23 J=1,NDIM
      XO=D(1,J)
      XN=D(2,J)
      D(1,J)=(XO+XN)/2.d0
      DT(J)=D(1,J)
      DO 22 I=2,NDM
      D(I,J)=XO+XN
      XO=XN
      XN=D(I+1,J)
      D(I,J)=(D(I,J)+XN)/3.d0
22    DT(J)=DT(J)+D(I,J)
      D(ND,J)=(XN+XO)/2.d0
23    DT(J)=DT(J)+D(ND,J)
      DO 28 J=1,NDIM
      RC=0.d0
      DO 24 I=1,ND
      R(I)=0.d0
      IF (DT(J).GE.1.0D18) THEN
       WRITE(6,*) '************** A SINGULARITY >1.0D18'
      END IF    
      IF(D(I,J).LE.1.0D-18) GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.d0
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/(R(K)+1.0d-30)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
      IF(IT.LT.ITMX.AND.ACC*DABS(AVGI).LT.SD) GO TO 9
      RETURN
      END
