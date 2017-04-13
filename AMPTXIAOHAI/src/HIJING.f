        SUBROUTINE HIJING(FRAME,BMIN0,BMAX0)
        PARAMETER (MAXPTN=400001)
        PARAMETER (MAXSTR=150001)
        PARAMETER (MAXIDL=4001)
        DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
        DOUBLE PRECISION  GX5, GY5, GZ5, FT5, PX5, PY5, PZ5, E5, XMASS5
        DOUBLE PRECISION  ATAUI, ZT1, ZT2, ZT3
        DOUBLE PRECISION  xnprod,etprod,xnfrz,etfrz,
     & dnprod,detpro,dnfrz,detfrz
        DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
        CHARACTER FRAME*8
        DIMENSION SCIP(300,300),RNIP(300,300),SJIP(300,300),JTP(3),
     &                        IPCOL(90000),ITCOL(90000)
        COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
        COMMON/hjcrdn/YP(3,300),YT(3,300)
        COMMON/HJGLBR/NELT,NINTHJ,NELP,NINP
        COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
        COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
        COMMON/HSTRNG/NFP(300,15),PP(300,15),NFT(300,15),PT(300,15)
        COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &                PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &                PJPM(300,500),NTJ(300),KFTJ(300,500),
     &                PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &                PJTE(300,500),PJTM(300,500)
        COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &       K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &       PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
        COMMON/HJJET4/NDR,IADR(MAXSTR,2),KFDR(MAXSTR),PDR(MAXSTR,5)
        common/xydr/rtdr(MAXSTR,2)
      COMMON/RNDF77/NSEED
        COMMON/LUJETS/N,K(9000,5),P(9000,5),V(9000,5)   
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON /ARPRC/ ITYPAR(MAXSTR),
     &       GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &       PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &       XMAR(MAXSTR)
        COMMON /PARA1/ MUL
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
        COMMON /prec2/GX5(MAXPTN),GY5(MAXPTN),GZ5(MAXPTN),FT5(MAXPTN),
     &       PX5(MAXPTN), PY5(MAXPTN), PZ5(MAXPTN), E5(MAXPTN),
     &       XMASS5(MAXPTN), ITYP5(MAXPTN)
        COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
        COMMON /ilist8/ LSTRG1(MAXPTN), LPART1(MAXPTN)
        COMMON /SREC1/ NSP, NST, NSI
        COMMON /SREC2/ATAUI(MAXSTR),ZT1(MAXSTR),ZT2(MAXSTR),ZT3(MAXSTR)
        COMMON /frzout/ xnprod(30),etprod(30),xnfrz(30),etfrz(30),
     & dnprod(30),detpro(30),dnfrz(30),detfrz(30)
      common/anim/nevent,isoft,isflag,izpc
      DOUBLE PRECISION PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
        COMMON /NOPREC/ NNOZPC, ITYPN(MAXIDL),
     &       GXN(MAXIDL), GYN(MAXIDL), GZN(MAXIDL), FTN(MAXIDL),
     &       PXN(MAXIDL), PYN(MAXIDL), PZN(MAXIDL), EEN(MAXIDL),
     &       XMN(MAXIDL)
        common /lastt/itimeh,bimp
        COMMON /AREVT/ IAEVT, IARUN, MISS
        common/phidcy/iphidcy,pttrig,ntrig,maxmiss,ipi0dcy
        common /para7/ ioscar,nsmbbbar,nsmmeson
        common /phiHJ/iphirp,phiRP
        common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1       xstrg0(MAXPTN),ystrg0(MAXPTN),
     2       xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
        SAVE   
        BMAX=MIN(BMAX0,HIPR1(34)+HIPR1(35))
        BMIN=MIN(BMIN0,BMAX)
        IF(IHNT2(1).LE.1 .AND. IHNT2(3).LE.1) THEN
                BMIN=0.0
                BMAX=2.5*SQRT(HIPR1(31)*0.1/HIPR1(40))
        ENDIF
        YP(1,1)=0.0
        YP(2,1)=0.0
        YP(3,1)=0.0
        IF(IHNT2(1).LE.1) GO TO 14
        DO 10 KP=1,IHNT2(1)
5        R=HIRND(1)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
        PHI=RANART(NSEED)*2.0*HIPR1(40)
        YP(1,KP)=R*SX*COS(PHI)
        YP(2,KP)=R*SX*SIN(PHI)
        YP(3,KP)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 10
        DO 8  KP2=1,KP-1
                DNBP1=(YP(1,KP)-YP(1,KP2))**2
                DNBP2=(YP(2,KP)-YP(2,KP2))**2
                DNBP3=(YP(3,KP)-YP(3,KP2))**2
                DNBP=DNBP1+DNBP2+DNBP3
                IF(DNBP.LT.HIPR1(29)*HIPR1(29)) GO TO 5
8        CONTINUE
10        CONTINUE
        if(IHNT2(1).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
           R=R/2.
           YP(1,1)=R*SX*COS(PHI)
           YP(2,1)=R*SX*SIN(PHI)
           YP(3,1)=R*CX
           YP(1,2)=-YP(1,1)
           YP(2,2)=-YP(2,1)
           YP(3,2)=-YP(3,1)
        endif
        DO 12 I=1,IHNT2(1)-1
        DO 12 J=I+1,IHNT2(1)
        IF(YP(3,I).GT.YP(3,J)) GO TO 12
        Y1=YP(1,I)
        Y2=YP(2,I)
        Y3=YP(3,I)
        YP(1,I)=YP(1,J)
        YP(2,I)=YP(2,J)
        YP(3,I)=YP(3,J)
        YP(1,J)=Y1
        YP(2,J)=Y2
        YP(3,J)=Y3
12        CONTINUE
14        YT(1,1)=0.0
        YT(2,1)=0.0
        YT(3,1)=0.0
        IF(IHNT2(3).LE.1) GO TO 24
        DO 20 KT=1,IHNT2(3)
15        R=HIRND(2)
        X=RANART(NSEED)
        CX=2.0*X-1.0
        SX=SQRT(1.0-CX*CX)
        PHI=RANART(NSEED)*2.0*HIPR1(40)
        YT(1,KT)=R*SX*COS(PHI)
        YT(2,KT)=R*SX*SIN(PHI)
        YT(3,KT)=R*CX
        IF(HIPR1(29).EQ.0.0) GO TO 20
        DO 18  KT2=1,KT-1
                DNBT1=(YT(1,KT)-YT(1,KT2))**2
                DNBT2=(YT(2,KT)-YT(2,KT2))**2
                DNBT3=(YT(3,KT)-YT(3,KT2))**2
                DNBT=DNBT1+DNBT2+DNBT3
                IF(DNBT.LT.HIPR1(29)*HIPR1(29)) GO TO 15
18        CONTINUE
20        CONTINUE
        if(IHNT2(3).EQ.2) then
           rnd1=max(RANART(NSEED),1.0e-20)
           rnd2=max(RANART(NSEED),1.0e-20)
           rnd3=max(RANART(NSEED),1.0e-20)
           R=-(log(rnd1)*4.38/2.0+log(rnd2)*0.85/2.0
     &          +4.38*0.85*log(rnd3)/(4.38+0.85))
           X=RANART(NSEED)
           CX=2.0*X-1.0
           SX=SQRT(1.0-CX*CX)
           PHI=RANART(NSEED)*2.0*HIPR1(40)
           R=R/2.
           YT(1,1)=R*SX*COS(PHI)
           YT(2,1)=R*SX*SIN(PHI)
           YT(3,1)=R*CX
           YT(1,2)=-YT(1,1)
           YT(2,2)=-YT(2,1)
           YT(3,2)=-YT(3,1)
        endif
        DO 22 I=1,IHNT2(3)-1
        DO 22 J=I+1,IHNT2(3)
        IF(YT(3,I).LT.YT(3,J)) GO TO 22
        Y1=YT(1,I)
        Y2=YT(2,I)
        Y3=YT(3,I)
        YT(1,I)=YT(1,J)
        YT(2,I)=YT(2,J)
        YT(3,I)=YT(3,J)
        YT(1,J)=Y1
        YT(2,J)=Y2
        YT(3,J)=Y3
22        CONTINUE
24        MISS=-1
50        MISS=MISS+1
        IF(MISS.GT.maxmiss) THEN
           WRITE(6,*) 'infinite loop happened in  HIJING'
           STOP
        ENDIF
        itest=0
        NATT=0
        JATT=0
        EATT=0.0
        CALL HIJINI
        NLOP=0
60        NT=0
        NP=0
        N0=0
        N01=0
        N10=0
        N11=0
        NELT=0
        NINTHJ=0
        NELP=0
        NINP=0
        NSG=0
        NCOLT=0
        BB=SQRT(BMIN**2+RANART(NSEED)*(BMAX**2-BMIN**2))
        PHI=0.
        if(iphirp.eq.1) PHI=2.0*HIPR1(40)*RANART(NSEED)
        phiRP=phi
        BBX=BB*COS(PHI)
        BBY=BB*SIN(PHI)
        HINT1(19)=BB
        HINT1(20)=PHI
        DO 70 JP=1,IHNT2(1)
        DO 70 JT=1,IHNT2(3)
           SCIP(JP,JT)=-1.0
           B2=(YP(1,JP)+BBX-YT(1,JT))**2+(YP(2,JP)+BBY-YT(2,JT))**2
           R2=B2*HIPR1(40)/HIPR1(31)/0.1
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)
     &          /1.2**2/REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)
     &          /1.2**2/REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(18)=HINT1(14)-APHX1*HINT1(15)
     &                        -APHX2*HINT1(16)+APHX1*APHX2*HINT1(17)
           IF(IHPR2(14).EQ.0.OR.
     &          (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) THEN
              GS=1.0-EXP(-(HIPR1(30)+HINT1(18))*ROMG(R2)/HIPR1(31))
              RANTOT=RANART(NSEED)
              IF(RANTOT.GT.GS) GO TO 70
              GO TO 65
           ENDIF
           GSTOT0=2.0*(1.0-EXP(-(HIPR1(30)+HINT1(18))
     &             /HIPR1(31)/2.0*ROMG(0.0)))
           R2=R2/GSTOT0
           GS=1.0-EXP(-(HIPR1(30)+HINT1(18))/HIPR1(31)*ROMG(R2))
           GSTOT=2.0*(1.0-SQRT(1.0-GS))
           RANTOT=RANART(NSEED)*GSTOT0
           IF(RANTOT.GT.GSTOT) GO TO 70
           IF(RANTOT.GT.GS) THEN
              CALL HIJCSC(JP,JT)
              GO TO 70
           ENDIF
 65           SCIP(JP,JT)=R2
           RNIP(JP,JT)=RANTOT
           SJIP(JP,JT)=HINT1(18)
           NCOLT=NCOLT+1
           IPCOL(NCOLT)=JP
           ITCOL(NCOLT)=JT
70        CONTINUE
        bimp=bb
        write(6,*) '#impact parameter,nlop,ncolt=',bimp,nlop,ncolt
        IF(NCOLT.EQ.0) THEN
           NLOP=NLOP+1
           IF(NLOP.LE.20.OR.
     &           (IHNT2(1).EQ.1.AND.IHNT2(3).EQ.1)) GO TO 60
           RETURN
        ENDIF
        IF(IHPR2(3).NE.0) THEN
           NHARD=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
           NHARD=MIN(NHARD,NCOLT)
           JPHARD=IPCOL(NHARD)
           JTHARD=ITCOL(NHARD)
        ENDIF
        IF(IHPR2(9).EQ.1) THEN
                NMINI=1+INT(RANART(NSEED)*(NCOLT-1)+0.5)
                NMINI=MIN(NMINI,NCOLT)
                JPMINI=IPCOL(NMINI)
                JTMINI=ITCOL(NMINI)
        ENDIF
        DO 200 JP=1,IHNT2(1)
        DO 200 JT=1,IHNT2(3)
        IF(SCIP(JP,JT).EQ.-1.0) GO TO 200
                NFP(JP,11)=NFP(JP,11)+1
                NFT(JT,11)=NFT(JT,11)+1
        IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).GT.1) THEN
                NP=NP+1
                N01=N01+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).LE.1) THEN
                NT=NT+1
                N10=N10+1
        ELSE IF(NFP(JP,5).LE.1 .AND. NFT(JT,5).LE.1) THEN
                NP=NP+1
                NT=NT+1
                N0=N0+1
        ELSE IF(NFP(JP,5).GT.1 .AND. NFT(JT,5).GT.1) THEN
                N11=N11+1
        ENDIF
        JOUT=0
        NFP(JP,10)=0
        NFT(JT,10)=0
        IF(IHPR2(8).EQ.0 .AND. IHPR2(3).EQ.0) GO TO 160
        IF(NFP(JP,6).LT.0 .OR. NFT(JT,6).LT.0) GO TO 160
        R2=SCIP(JP,JT)
        HINT1(18)=SJIP(JP,JT)
        TT=ROMG(R2)*HINT1(18)/HIPR1(31)
        TTS=HIPR1(30)*ROMG(R2)/HIPR1(31)
        NJET=0
        IF(IHPR2(3).NE.0 .AND. JP.EQ.JPHARD .AND. JT.EQ.JTHARD) THEN
           CALL JETINI(JP,JT,1)
           CALL HIJHRD(JP,JT,0,JFLG,0)
           HINT1(26)=HINT1(47)
           HINT1(27)=HINT1(48)
           HINT1(28)=HINT1(49)
           HINT1(29)=HINT1(50)
           HINT1(36)=HINT1(67)
           HINT1(37)=HINT1(68)
           HINT1(38)=HINT1(69)
           HINT1(39)=HINT1(70)
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                                JFLG.GE.3) IASG(NSG,3)=1
           IHNT2(9)=IHNT2(14)
           IHNT2(10)=IHNT2(15)
           DO 105 I05=1,5
              HINT1(20+I05)=HINT1(40+I05)
              HINT1(30+I05)=HINT1(50+I05)
 105           CONTINUE
           JOUT=1
           IF(IHPR2(8).EQ.0) GO TO 160
           RRB1=MIN((YP(1,JP)**2+YP(2,JP)**2)/1.2**2
     &                /REAL(IHNT2(1))**0.6666667,1.0)
           RRB2=MIN((YT(1,JT)**2+YT(2,JT)**2)/1.2**2
     &                /REAL(IHNT2(3))**0.6666667,1.0)
           APHX1=HIPR1(6)*4.0/3.0*(IHNT2(1)**0.3333333-1.0)
     &           *SQRT(1.0-RRB1)
           APHX2=HIPR1(6)*4.0/3.0*(IHNT2(3)**0.3333333-1.0)
     &           *SQRT(1.0-RRB2)
           HINT1(65)=HINT1(61)-APHX1*HINT1(62)
     &                        -APHX2*HINT1(63)+APHX1*APHX2*HINT1(64)
           TTRIG=ROMG(R2)*HINT1(65)/HIPR1(31)
           NJET=-1
           XR1=-ALOG(EXP(-TTRIG)+RANART(NSEED)*(1.0-EXP(-TTRIG)))
 106           NJET=NJET+1
           XR1=XR1-ALOG(max(RANART(NSEED),1.0e-20))
           IF(XR1.LT.TTRIG) GO TO 106
           XR=0.0
 107           NJET=NJET+1
           XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
           IF(XR.LT.TT-TTRIG) GO TO 107
           NJET=NJET-1
           GO TO 112
        ENDIF
        IF(IHPR2(9).EQ.1.AND.JP.EQ.JPMINI.AND.JT.EQ.JTMINI) GO TO 110
        IF(IHPR2(8).GT.0 .AND.RNIP(JP,JT).LE.EXP(-TT)*
     &                 (1.0-EXP(-TTS))) GO TO 160
110        XR=-ALOG(EXP(-TT)+RANART(NSEED)*(1.0-EXP(-TT)))
111        NJET=NJET+1
        XR=XR-ALOG(max(RANART(NSEED),1.0e-20))
        IF(XR.LT.TT) GO TO 111
112        NJET=MIN(NJET,IHPR2(8))
        IF(IHPR2(8).LT.0)  NJET=ABS(IHPR2(8))
        DO 150 ijet=1,NJET
           CALL JETINI(JP,JT,0)
           CALL HIJHRD(JP,JT,JOUT,JFLG,1)
           IF(JFLG.EQ.0) GO TO 160
           IF(JFLG.LT.0) THEN
              IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJHRD'
              GO TO 50
           ENDIF
           JOUT=JOUT+1
           IF(ABS(HINT1(46)).GT.HIPR1(11).AND.JFLG.EQ.2) NFP(JP,7)=1
           IF(ABS(HINT1(56)).GT.HIPR1(11).AND.JFLG.EQ.2) NFT(JT,7)=1
           IF(MAX(ABS(HINT1(46)),ABS(HINT1(56))).GT.HIPR1(11).AND.
     &                        JFLG.GE.3) IASG(NSG,3)=1
 150        CONTINUE
 160        CONTINUE
        CALL HIJSFT(JP,JT,JOUT,IERROR)
        IF(IERROR.NE.0) THEN
           IF(IHPR2(10).NE.0) WRITE(6,*) 'error occured in HIJSFT'
           GO TO 50
        ENDIF
        JATT=JATT+JOUT
200        CONTINUE
           call minijet_out(BB,phiRP)
           if(pttrig.gt.0.and.ntrig.eq.0) goto 50
        DO 201 JP=1,IHNT2(1)
           IF(NFP(JP,5).GT.2) THEN
              NINP=NINP+1
           ELSE IF(NFP(JP,5).EQ.2.OR.NFP(JP,5).EQ.1) THEN
              NELP=NELP+1
           ENDIF
 201    continue
        DO 202 JT=1,IHNT2(3)
           IF(NFT(JT,5).GT.2) THEN
              NINTHJ=NINTHJ+1
           ELSE IF(NFT(JT,5).EQ.2.OR.NFT(JT,5).EQ.1) THEN
              NELT=NELT+1
           ENDIF
 202    continue
        IF((IHPR2(8).NE.0.OR.IHPR2(3).NE.0).AND.IHPR2(4).GT.0.AND.
     &                        IHNT2(1).GT.1.AND.IHNT2(3).GT.1) THEN
                DO 271 I=1,IHNT2(1)
                        IF(NFP(I,7).EQ.1) CALL QUENCH(I,1)
271                CONTINUE
                DO 272 I=1,IHNT2(3)
                        IF(NFT(I,7).EQ.1) CALL QUENCH(I,2)
272                CONTINUE
                DO 273 ISG=1,NSG
                        IF(IASG(ISG,3).EQ.1) CALL QUENCH(ISG,3)
273                CONTINUE
        ENDIF
        if(isoft.eq.1) then
           isflag=1
        NSP = IHNT2(1)
        NST = IHNT2(3)
        NSI = NSG
        ISTR = 0
        NPAR = 0
        DO 1008 I = 1, IHNT2(1)
           ISTR = ISTR + 1
           DO 1007 J = 1, NPJ(I)
              IF (KFPJ(I, J) .EQ. 21) THEN
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = KFPJ(I, J)
              GX0(NPAR) = dble(YP(1, I)+0.5*BB*cos(phiRP))
              GY0(NPAR) = dble(YP(2, I)+0.5*BB*sin(phiRP))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PJPX(I, J))
              PY0(NPAR) = dble(PJPY(I, J))
              PZ0(NPAR) = dble(PJPZ(I, J))
              XMASS0(NPAR) = dble(PJPM(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
              END IF
 1007      CONTINUE
 1008   CONTINUE
        DO 1010 I = 1, IHNT2(3)
           ISTR = ISTR + 1
           DO 1009 J = 1, NTJ(I)
              IF (KFTJ(I, J) .EQ. 21) THEN
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = KFTJ(I, J)
              GX0(NPAR) = dble(YT(1, I)-0.5*BB*cos(phiRP))
              GY0(NPAR) = dble(YT(2, I)-0.5*BB*sin(phiRP))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PJTX(I, J))
              PY0(NPAR) = dble(PJTY(I, J))
              PZ0(NPAR) = dble(PJTZ(I, J))
              XMASS0(NPAR) = dble(PJTM(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
              END IF
 1009      CONTINUE
 1010   CONTINUE
        DO 1012 I = 1, NSG
           ISTR = ISTR + 1
           DO 1011 J = 1, NJSG(I)
              IF (K2SG(I, J) .EQ. 21) THEN
              NPAR = NPAR + 1
              LSTRG0(NPAR) = ISTR
              LPART0(NPAR) = J
              ITYP0(NPAR) = K2SG(I, J)
              GX0(NPAR) = 0.5d0 * 
     1             dble(YP(1, IASG(I, 1)) + YT(1, IASG(I, 2)))
              GY0(NPAR) = 0.5d0 * 
     2             dble(YP(2, IASG(I, 1)) + YT(2, IASG(I, 2)))
              GZ0(NPAR) = 0d0
              FT0(NPAR) = 0d0
              PX0(NPAR) = dble(PXSG(I, J))
              PY0(NPAR) = dble(PYSG(I, J))
              PZ0(NPAR) = dble(PZSG(I, J))
              XMASS0(NPAR) = dble(PMSG(I, J))
              E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1             +PZ0(NPAR)**2+XMASS0(NPAR)**2)
              END IF
 1011      CONTINUE
 1012   CONTINUE
        MUL = NPAR
        CALL HJANA1
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
        CALL ZPCMN
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        DO 1013 I = 1, MUL
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
 1013   CONTINUE
 210    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))
 211    format(I6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))
 395    format(3I8,f10.4,4I5)
        itest=itest+1
        DO 1014 I = 1, MUL
           IF (LSTRG1(I) .LE. NSP) THEN
              NSTRG = LSTRG1(I)
              NPART = LPART1(I)
              KFPJ(NSTRG, NPART) = ITYP5(I)
              PJPX(NSTRG, NPART) = sngl(PX5(I))
              PJPY(NSTRG, NPART) = sngl(PY5(I))
              PJPZ(NSTRG, NPART) = sngl(PZ5(I))
              PJPE(NSTRG, NPART) = sngl(E5(I))
              PJPM(NSTRG, NPART) = sngl(XMASS5(I))
           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
              NSTRG = LSTRG1(I) - NSP
              NPART = LPART1(I)
              KFTJ(NSTRG, NPART) = ITYP5(I)
              PJTX(NSTRG, NPART) = sngl(PX5(I))
              PJTY(NSTRG, NPART) = sngl(PY5(I))
              PJTZ(NSTRG, NPART) = sngl(PZ5(I))
              PJTE(NSTRG, NPART) = sngl(E5(I))
              PJTM(NSTRG, NPART) = sngl(XMASS5(I))
           ELSE
              NSTRG = LSTRG1(I) - NSP - NST
              NPART = LPART1(I)
              K2SG(NSTRG, NPART) = ITYP5(I)
              PXSG(NSTRG, NPART) = sngl(PX5(I))
              PYSG(NSTRG, NPART) = sngl(PY5(I))
              PZSG(NSTRG, NPART) = sngl(PZ5(I))
              PESG(NSTRG, NPART) = sngl(E5(I))
              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
           END IF
 1014   CONTINUE
        CALL HJANA2
        elseif(isoft.eq.2) then
        NSP = IHNT2(1)
        NST = IHNT2(3)
        NSI = NSG
        NPAR=0
        ISTR=0
        MSTJ(1)=0
        IHPR2(1)=0
        isflag=0
        IF(IHPR2(20).NE.0) THEN
           DO 320 NTP=1,2
              DO 310 jjtp=1,IHNT2(2*NTP-1)
                 ISTR = ISTR + 1
                 CALL HIJFRG(jjtp,NTP,IERROR)
                 if(NTP.eq.1) then
                    NPJ(jjtp)=MAX0(N-2,0)
                 else
                    NTJ(jjtp)=MAX0(N-2,0)
                 endif
                 do 300 ii=1,N
                 NPAR = NPAR + 1
                 LSTRG0(NPAR) = ISTR
                 LPART0(NPAR) = II
                 ITYP0(NPAR) = K(II,2)
                 GZ0(NPAR) = 0d0
                 FT0(NPAR) = 0d0
                 PX0(NPAR) = dble(P(II,1))
                 PY0(NPAR) = dble(P(II,2))
                 PZ0(NPAR) = dble(P(II,3))
                 XMASS0(NPAR) = dble(P(II,5))
                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
                 IF (NTP .EQ. 1) THEN
                    GX0(NPAR) = dble(YP(1, jjtp)+0.5*BB*cos(phiRP))
                    GY0(NPAR) = dble(YP(2, jjtp)+0.5*BB*sin(phiRP))
                    IITYP=ITYP0(NPAR)
                    nstrg=LSTRG0(NPAR)
                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (II.eq.1.or.II.eq.N)) then
                       PP(nstrg,6)=sngl(PX0(NPAR))
                       PP(nstrg,7)=sngl(PY0(NPAR))
                       PP(nstrg,14)=sngl(XMASS0(NPAR))
                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(II.eq.1.or.II.eq.N)) then
                       PP(nstrg,8)=sngl(PX0(NPAR))
                       PP(nstrg,9)=sngl(PY0(NPAR))
                       PP(nstrg,15)=sngl(XMASS0(NPAR))
                    else
                       NPART = LPART0(NPAR)-1
                       KFPJ(NSTRG, NPART) = ITYP0(NPAR)
                       PJPX(NSTRG, NPART) = sngl(PX0(NPAR))
                       PJPY(NSTRG, NPART) = sngl(PY0(NPAR))
                       PJPZ(NSTRG, NPART) = sngl(PZ0(NPAR))
                       PJPE(NSTRG, NPART) = sngl(E0(NPAR))
                       PJPM(NSTRG, NPART) = sngl(XMASS0(NPAR))
                    endif
                 ELSE
                    GX0(NPAR) = dble(YT(1, jjtp)-0.5*BB*cos(phiRP))
                    GY0(NPAR) = dble(YT(2, jjtp)-0.5*BB*sin(phiRP))
                    IITYP=ITYP0(NPAR)
                    nstrg=LSTRG0(NPAR)-NSP
                    if(IITYP.eq.2112.or.IITYP.eq.2212) then
                    elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (II.eq.1.or.II.eq.N)) then
                       PT(nstrg,6)=sngl(PX0(NPAR))
                       PT(nstrg,7)=sngl(PY0(NPAR))
                       PT(nstrg,14)=sngl(XMASS0(NPAR))
                    elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(II.eq.1.or.II.eq.N)) then
                       PT(nstrg,8)=sngl(PX0(NPAR))
                       PT(nstrg,9)=sngl(PY0(NPAR))
                       PT(nstrg,15)=sngl(XMASS0(NPAR))
                    else
                       NPART = LPART0(NPAR)-1
                       KFTJ(NSTRG, NPART) = ITYP0(NPAR)
                       PJTX(NSTRG, NPART) = sngl(PX0(NPAR))
                       PJTY(NSTRG, NPART) = sngl(PY0(NPAR))
                       PJTZ(NSTRG, NPART) = sngl(PZ0(NPAR))
                       PJTE(NSTRG, NPART) = sngl(E0(NPAR))
                       PJTM(NSTRG, NPART) = sngl(XMASS0(NPAR))
                    endif
                 END IF
 300          continue
 310          continue
 320       continue
           DO 330 ISG=1,NSG
              ISTR = ISTR + 1
              CALL HIJFRG(ISG,3,IERROR)
              NJSG(ISG)=N
              do 1001 ii=1,N
                 NPAR = NPAR + 1
                 LSTRG0(NPAR) = ISTR
                 LPART0(NPAR) = II
                 ITYP0(NPAR) = K(II,2)
                 GX0(NPAR)=0.5d0*
     1                dble(YP(1,IASG(ISG,1))+YT(1,IASG(ISG,2)))
                 GY0(NPAR)=0.5d0*
     2                dble(YP(2,IASG(ISG,1))+YT(2,IASG(ISG,2)))
                 GZ0(NPAR) = 0d0
                 FT0(NPAR) = 0d0
                 PX0(NPAR) = dble(P(II,1))
                 PY0(NPAR) = dble(P(II,2))
                 PZ0(NPAR) = dble(P(II,3))
                 XMASS0(NPAR) = dble(P(II,5))
                 E0(NPAR) = dsqrt(PX0(NPAR)**2+PY0(NPAR)**2
     1                +PZ0(NPAR)**2+XMASS0(NPAR)**2)
 1001         continue
 330       continue
        endif
        MUL = NPAR
        CALL HJANA1
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
        CALL ZPCMN
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        itest=itest+1
        DO 1015 I = 1, MUL
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
 1015   CONTINUE
        do 1004 nmom=1,5
           do 1002 nstrg=1,nsp
              PP(nstrg,nmom)=0.
 1002      continue
           do 1003 nstrg=1,nst
              PT(nstrg,nmom)=0.
 1003      continue
 1004   continue
        DO 1005 I = 1, MUL
           IITYP=ITYP5(I)
           IF (LSTRG1(I) .LE. NSP) THEN
              NSTRG = LSTRG1(I)
              if(IITYP.eq.2112.or.IITYP.eq.2212) then
                 PP(nstrg,1)=sngl(PX5(I))
                 PP(nstrg,2)=sngl(PY5(I))
                 PP(nstrg,3)=sngl(PZ5(I))
                 PP(nstrg,4)=sngl(E5(I))
                 PP(nstrg,5)=sngl(XMASS5(I))
              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
                 PP(nstrg,6)=sngl(PX5(I))
                 PP(nstrg,7)=sngl(PY5(I))
                 PP(nstrg,14)=sngl(XMASS5(I))
                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
              elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NPJ(NSTRG)+2))) then
                 PP(nstrg,8)=sngl(PX5(I))
                 PP(nstrg,9)=sngl(PY5(I))
                 PP(nstrg,15)=sngl(XMASS5(I))
                 PP(nstrg,1)=PP(nstrg,1)+sngl(PX5(I))
                 PP(nstrg,2)=PP(nstrg,2)+sngl(PY5(I))
                 PP(nstrg,3)=PP(nstrg,3)+sngl(PZ5(I))
                 PP(nstrg,4)=PP(nstrg,4)+sngl(E5(I))
                 PP(nstrg,5)=sqrt(PP(nstrg,4)**2-PP(nstrg,1)**2
     1                -PP(nstrg,2)**2-PP(nstrg,3)**2)
              else
                 NPART = LPART1(I)-1
                 KFPJ(NSTRG, NPART) = ITYP5(I)
                 PJPX(NSTRG, NPART) = sngl(PX5(I))
                 PJPY(NSTRG, NPART) = sngl(PY5(I))
                 PJPZ(NSTRG, NPART) = sngl(PZ5(I))
                 PJPE(NSTRG, NPART) = sngl(E5(I))
                 PJPM(NSTRG, NPART) = sngl(XMASS5(I))
              endif
           ELSE IF (LSTRG1(I) .LE. NSP + NST) THEN
              NSTRG = LSTRG1(I) - NSP
              if(IITYP.eq.2112.or.IITYP.eq.2212) then
                 PT(nstrg,1)=sngl(PX5(I))
                 PT(nstrg,2)=sngl(PY5(I))
                 PT(nstrg,3)=sngl(PZ5(I))
                 PT(nstrg,4)=sngl(E5(I))
                 PT(nstrg,5)=sngl(XMASS5(I))
              elseif((IITYP.eq.1.or.IITYP.eq.2).and.
     1 (LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
                 PT(nstrg,6)=sngl(PX5(I))
                 PT(nstrg,7)=sngl(PY5(I))
                 PT(nstrg,14)=sngl(XMASS5(I))
                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
              elseif((IITYP.eq.1103.or.IITYP.eq.2101
     1 .or.IITYP.eq.2103.or.IITYP.eq.2203.
     2 .or.IITYP.eq.3101.or.IITYP.eq.3103.
     3 .or.IITYP.eq.3201.or.IITYP.eq.3203.or.IITYP.eq.3303)
     4 .and.(LPART1(I).eq.1.or.LPART1(I).eq.(NTJ(NSTRG)+2))) then
                 PT(nstrg,8)=sngl(PX5(I))
                 PT(nstrg,9)=sngl(PY5(I))
                 PT(nstrg,15)=sngl(XMASS5(I))
                 PT(nstrg,1)=PT(nstrg,1)+sngl(PX5(I))
                 PT(nstrg,2)=PT(nstrg,2)+sngl(PY5(I))
                 PT(nstrg,3)=PT(nstrg,3)+sngl(PZ5(I))
                 PT(nstrg,4)=PT(nstrg,4)+sngl(E5(I))
                 PT(nstrg,5)=sqrt(PT(nstrg,4)**2-PT(nstrg,1)**2
     1                -PT(nstrg,2)**2-PT(nstrg,3)**2)
              else
                 NPART = LPART1(I)-1
                 KFTJ(NSTRG, NPART) = ITYP5(I)
                 PJTX(NSTRG, NPART) = sngl(PX5(I))
                 PJTY(NSTRG, NPART) = sngl(PY5(I))
                 PJTZ(NSTRG, NPART) = sngl(PZ5(I))
                 PJTE(NSTRG, NPART) = sngl(E5(I))
                 PJTM(NSTRG, NPART) = sngl(XMASS5(I))
              endif
           ELSE
              NSTRG = LSTRG1(I) - NSP - NST
              NPART = LPART1(I)
              K2SG(NSTRG, NPART) = ITYP5(I)
              PXSG(NSTRG, NPART) = sngl(PX5(I))
              PYSG(NSTRG, NPART) = sngl(PY5(I))
              PZSG(NSTRG, NPART) = sngl(PZ5(I))
              PESG(NSTRG, NPART) = sngl(E5(I))
              PMSG(NSTRG, NPART) = sngl(XMASS5(I))
           END IF
 1005   CONTINUE
        MSTJ(1)=1
        IHPR2(1)=1
        isflag=1
        HIPR1(1)=0.94
        CALL HJANA2
        elseif(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
        isflag=0
        IF(IHPR2(20).NE.0) THEN
           DO 560 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 551                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  551
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
                nsbstR=0
                DO 560 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 560
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   GXAR(NATT) = 0.5 * (YP(1, IASG(ISG, 1)) +
     &                YT(1, IASG(ISG, 2)))
                   GYAR(NATT) = 0.5 * (YP(2, IASG(ISG, 1)) +
     &                YT(2, IASG(ISG, 2)))
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
                   istrg0(NATT)=ISG
 560            CONTINUE
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 600 NTP=1,2
           DO 600 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
 581                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  581
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 590 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 590
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   IF (NTP .EQ. 1) THEN
                      GXAR(NATT) = YP(1, jjtp)+0.5*BB*cos(phiRP)
                      GYAR(NATT) = YP(2, jjtp)+0.5*BB*sin(phiRP)
                   ELSE
                      GXAR(NATT) = YT(1, jjtp)-0.5*BB*cos(phiRP)
                      GYAR(NATT) = YT(2, jjtp)-0.5*BB*sin(phiRP)
                   END IF
                   GZAR(NATT) = 0.
                   FTAR(NATT) = 0.
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
                   xstrg0(NATT)=dble(GXAR(NATT))
                   ystrg0(NATT)=dble(GYAR(NATT))
                   istrg0(NATT)=NTP*10000+jjtp
 590                CONTINUE 
 600           CONTINUE
        ENDIF
        if(NDR.ge.1) then
        DO 650 I=1,NDR
                NATT=NATT+1
                KATT(NATT,1)=KFDR(I)
                KATT(NATT,2)=40
                KATT(NATT,3)=0
                PATT(NATT,1)=PDR(I,1)
                PATT(NATT,2)=PDR(I,2)
                PATT(NATT,3)=PDR(I,3)
                PATT(NATT,4)=PDR(I,4)
                EATT=EATT+PDR(I,4)
                GXAR(NATT) = rtdr(I,1)
                GYAR(NATT) = rtdr(I,2)
                GZAR(NATT) = 0.
                FTAR(NATT) = 0.
                ITYPAR(NATT) =KATT(NATT,1) 
                PXAR(NATT) = PATT(NATT,1)
                PYAR(NATT) = PATT(NATT,2)
                PZAR(NATT) = PATT(NATT,3)
                PEAR(NATT) = PATT(NATT,4)
                XMAR(NATT) = PDR(I,5)
 650        CONTINUE
         endif
         call embedHighPt
        CALL HJANA1
        call htop
        nsp=0
        nst=0
        nsg=natt
        NSI=NSG
        if(ioscar.eq.3) WRITE (95, *) IAEVT, mul
        CALL ZPCMN
        WRITE (14, 395) IAEVT, MISS, MUL, bimp, NELP,NINP,NELT,NINTHJ
        itest=itest+1
        DO 1016 I = 1, MUL
           if(dmax1(abs(GX5(I)),abs(GY5(I)),abs(GZ5(I)),abs(FT5(I)))
     1          .lt.9999) then
              write(14,210) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           else
              write(14,211) ITYP5(I), PX5(I), PY5(I), PZ5(I), XMASS5(I),
     1             GX5(I), GY5(I), GZ5(I), FT5(I)
           endif
 1016   CONTINUE
        DO 1018 I = 1, MAXSTR
           DO 1017 J = 1, 3
              K1SGS(I, J) = 0
              K2SGS(I, J) = 0
              PXSGS(I, J) = 0d0
              PYSGS(I, J) = 0d0
              PZSGS(I, J) = 0d0
              PESGS(I, J) = 0d0
              PMSGS(I, J) = 0d0
              GXSGS(I, J) = 0d0
              GYSGS(I, J) = 0d0
              GZSGS(I, J) = 0d0
              FTSGS(I, J) = 0d0
 1017      CONTINUE
 1018   CONTINUE
        DO 1019 I = 1, MUL
           IITYP=ITYP5(I)
           NSTRG = LSTRG1(I)
           NPART = LPART1(I)
           K2SGS(NSTRG, NPART) = ITYP5(I)
           PXSGS(NSTRG, NPART) = PX5(I)
           PYSGS(NSTRG, NPART) = PY5(I)
           PZSGS(NSTRG, NPART) = PZ5(I)
           PMSGS(NSTRG, NPART) = XMASS5(I)
           E5(I)=dsqrt(PX5(I)**2+PY5(I)**2+PZ5(I)**2+XMASS5(I)**2)
           PESGS(NSTRG, NPART) = E5(I)
           GXSGS(NSTRG, NPART) = GX5(I)
           GYSGS(NSTRG, NPART) = GY5(I)
           GZSGS(NSTRG, NPART) = GZ5(I)
           FTSGS(NSTRG, NPART) = FT5(I)
 1019   CONTINUE
        CALL HJANA2
        endif
        if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
           NATT=0
           EATT=0.
           call ptoh
           do 1006 I=1,nnozpc
              NATT=NATT+1
              KATT(NATT,1)=ITYPN(I)
              PATT(NATT,1)=PXN(I)
              PATT(NATT,2)=PYN(I)
              PATT(NATT,3)=PZN(I)
              PATT(NATT,4)=EEN(I)
              EATT=EATT+EEN(I)
              GXAR(NATT)=GXN(I)
              GYAR(NATT)=GYN(I)
              GZAR(NATT)=GZN(I)
              FTAR(NATT)=FTN(I)
              ITYPAR(NATT)=ITYPN(I)
              PXAR(NATT)=PXN(I)
              PYAR(NATT)=PYN(I)
              PZAR(NATT)=PZN(I)
              PEAR(NATT)=EEN(I)
              XMAR(NATT)=XMN(I)
 1006      continue
           goto 565
        endif
        IF(IHPR2(20).NE.0) THEN
           DO 360 ISG=1,NSG
                CALL HIJFRG(ISG,3,IERROR)
                IF(MSTU(24).NE.0 .OR.IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
                      WRITE(6,*) 'error occured ISG, repeat the event'
                  write(6,*) ISG
                   ENDIF
                   GO TO 50
                ENDIF
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
351                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  351
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
                nsbstR=0
                DO 360 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 360
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=20
                   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   LSG = NSP + NST + ISG
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
360           CONTINUE
           JTP(1)=IHNT2(1)
           JTP(2)=IHNT2(3)
           DO 400 NTP=1,2
           DO 400 jjtp=1,JTP(NTP)
                CALL HIJFRG(jjtp,NTP,IERROR)
                IF(MSTU(24).NE.0 .OR. IERROR.GT.0) THEN
                   MSTU(24)=0
                   MSTU(28)=0
                   IF(IHPR2(10).NE.0) THEN
                  WRITE(6,*) 'error occured P&T, repeat the event'
                  WRITE(6,*) NTP,jjtp
                   ENDIF
                   GO TO 50
                ENDIF
                nsbst=1
                IDSTR=92
                IF(IHPR2(21).EQ.0) THEN
                   CALL LUEDIT(2)
                ELSE
381                   nsbst=nsbst+1
                   IF(K(nsbst,2).LT.91.OR.K(nsbst,2).GT.93) GO TO  381
                   IDSTR=K(nsbst,2)
                   nsbst=nsbst+1
                ENDIF
                IF(FRAME.EQ.'LAB') THEN
                        CALL HBOOST
                ENDIF
                NFTP=NFP(jjtp,5)
                IF(NTP.EQ.2) NFTP=10+NFT(jjtp,5)
                nsbstR=0
                DO 390 I=nsbst,N
                   IF(K(I,2).EQ.IDSTR) THEN
                      nsbstR=nsbstR+1
                      GO TO 390
                   ENDIF
                   K(I,4)=nsbstR
                   NATT=NATT+1
                   KATT(NATT,1)=K(I,2)
                   KATT(NATT,2)=NFTP
                   KATT(NATT,4)=K(I,1)
                   IF(K(I,3).EQ.0) THEN
                      KATT(NATT,3)=0
                   ELSEIF(K(I,3).ne.0.and.K(K(I,3),2).EQ.IDSTR) THEN
                      KATT(NATT,3)=0
                   ELSE
                      KATT(NATT,3)=NATT-I+K(I,3)+nsbstR-K(K(I,3),4)
                   ENDIF
                   PATT(NATT,1)=P(I,1)
                   PATT(NATT,2)=P(I,2)
                   PATT(NATT,3)=P(I,3)
                   PATT(NATT,4)=P(I,4)
                   EATT=EATT+P(I,4)
                   IF (NTP .EQ. 1) THEN
                      LSG = jjtp
                   ELSE
                      LSG = jjtp + NSP
                   END IF
                   GXAR(NATT) = sngl(ZT1(LSG))
                   GYAR(NATT) = sngl(ZT2(LSG))
                   GZAR(NATT) = sngl(ZT3(LSG))
                   FTAR(NATT) = sngl(ATAUI(LSG))
                   ITYPAR(NATT) = K(I, 2)
                   PXAR(NATT) = P(I, 1)
                   PYAR(NATT) = P(I, 2)
                   PZAR(NATT) = P(I, 3)
                   PEAR(NATT) = P(I, 4)
                   XMAR(NATT) = P(I, 5)
390                CONTINUE 
400           CONTINUE
        ENDIF
        DO 450 I=1,NDR
           NATT=NATT+1
           KATT(NATT,1)=KFDR(I)
           KATT(NATT,2)=40
           KATT(NATT,3)=0
           PATT(NATT,1)=PDR(I,1)
           PATT(NATT,2)=PDR(I,2)
           PATT(NATT,3)=PDR(I,3)
           PATT(NATT,4)=PDR(I,4)
           EATT=EATT+PDR(I,4)
           GXAR(NATT) = rtdr(I,1)
           GYAR(NATT) = rtdr(I,2)
           GZAR(NATT) = 0.
           FTAR(NATT) = 0.
           ITYPAR(NATT) =KATT(NATT,1) 
           PXAR(NATT) = PATT(NATT,1)
           PYAR(NATT) = PATT(NATT,2)
           PZAR(NATT) = PATT(NATT,3)
           PEAR(NATT) = PATT(NATT,4)
           XMAR(NATT) = PDR(I,5)
 450    CONTINUE
 565    continue
        DENGY=EATT/(IHNT2(1)*HINT1(6)+IHNT2(3)*HINT1(7))-1.0
        IF(ABS(DENGY).GT.HIPR1(43).AND.IHPR2(20).NE.0
     &     .AND.IHPR2(21).EQ.0) THEN
         IF(IHPR2(10).NE.0) 
     &        WRITE(6,*) 'Energy not conserved, repeat the event'
         write(6,*) 'violated:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
         GO TO 50
        ENDIF
        write(6,*) 'satisfied:EATT(GeV),NATT,B(fm)=',EATT,NATT,bimp
        write(6,*) ' '
        write(94,*) IAEVT,MISS,IHNT2(1),IHNT2(3),bimp
        DO JP=1,IHNT2(1)
           write(94,243) YP(1,JP)+0.5*BB*cos(phiRP), 
     1 YP(2,JP)+0.5*BB*sin(phiRP),JP, NFP(JP,5),yp(3,jp),
     2 NFP(JP,3),NFP(JP,4)
        ENDDO
        DO JT=1,IHNT2(3)
           write(94,243) YT(1,JT)-0.5*BB*cos(phiRP), 
     1 YT(2,JT)-0.5*BB*sin(phiRP), -JT, NFT(JT,5),yt(3,jt),
     2 NFT(JT,3),NFT(JT,4)
        ENDDO
 243    format(f10.3,1x,f10.3,2(1x,I5),1x,f10.3,2(1x,I5))
        RETURN
        END
