       subroutine flow(nt)
       PARAMETER ( PI=3.1415926,APion=0.13957,aka=0.498)
        PARAMETER   (MAXSTR=150001,MAXR=1,AMU= 0.9383,etaM=0.5475)
       DIMENSION ypion(-80:80),ypr(-80:80),ykaon(-80:80)
       dimension pxpion(-80:80),pxpro(-80:80),pxkaon(-80:80)
      COMMON  /AA/      R(3,MAXSTR)
      COMMON  /BB/      P(3,MAXSTR)
      COMMON  /CC/      E(MAXSTR)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      COMMON  /RR/      MASSR(0:MAXR)
      COMMON  /RUN/     NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      SAVE   
       ycut1=-2.6
       ycut2=2.6
       DY=0.2
       LY=NINT((YCUT2-YCUT1)/DY)
       do 11 kk=-80,80
       pxpion(kk)=0
       pxpro(kk)=0
       pxkaon(kk)=0
11       continue
       DO 701 J=-LY,LY
       ypion(j)=0
       ykaon(j)=0
       ypr(j)=0
  701   CONTINUE
       nkaon=0
       npr=0
       npion=0
          IS=0
          DO 20 NRUN=1,NUM
          IS=IS+MASSR(NRUN-1)
          DO 20 J=1,MASSR(NRUN)
          I=J+IS
       e00=sqrt(P(1,I)**2+P(2,i)**2+P(3,i)**2+e(I)**2)
       y00=0.5*alog((e00+p(3,i))/(e00-p(3,i)))
       if(abs(y00).ge.ycut2)go to 20
       iy=nint(y00/DY)
       if(abs(iy).ge.80)go to 20
       if(e(i).eq.0)go to 20
       if(lb(i).ge.25)go to 20
       if((lb(i).le.5).and.(lb(i).ge.3))go to 50
       if(lb(i).eq.1.or.lb(i).eq.2)go to 200
       if(lb(i).ge.6.and.lb(i).le.17)go to 200
       if(lb(i).eq.23)go to 400
       go to 20
50       npion=npion+1
        ypion(iy)=ypion(iy)+1
       pxpion(iy)=pxpion(iy)+p(1,i)/e(I)
       go TO 20
200      npr=npr+1  
                pxpro(iy)=pxpro(iy)+p(1,I)/E(I)
                 ypr(iy)=ypr(iy)+1.
        go to 20
400     nkaon=nkaon+1  
                 ykaon(iy)=ykaon(iy)+1.
                pxkaon(iy)=pxkaon(iy)+p(1,i)/E(i)
20      CONTINUE
       do 3 npt=-10,10
       IF(ypr(npt).eq.0) go to 101
       pxpro(NPT)=-Pxpro(NPT)/ypr(NPT)
       DNUC=Pxpro(NPT)/SQRT(ypr(NPT))
101       IF(ypion(npt).eq.0) go to 102
       pxpion(NPT)=-pxpion(NPT)/ypion(NPT)
       DNUCp=pxpion(NPT)/SQRT(ypion(NPT))
102       IF(ykaon(npt).eq.0) go to 3
       pxkaon(NPT)=-pxkaon(NPT)/ykaon(NPT)
       DNUCk=pxkaon(NPT)/SQRT(ykaon(NPT))
3       CONTINUE
       DO 1001 M=-LY,LY
       DYPR=0
       IF(YPR(M).NE.0)DYPR=SQRT(YPR(M))/FLOAT(NRUN)/DY
       YPR(M)=YPR(M)/FLOAT(NRUN)/DY
       DYPION=0
       IF(YPION(M).NE.0)DYPION=SQRT(YPION(M))/FLOAT(NRUN)/DY
       YPION(M)=YPION(M)/FLOAT(NRUN)/DY
       DYKAON=0
       IF(YKAON(M).NE.0)DYKAON=SQRT(YKAON(M))/FLOAT(NRUN)/DY
       YKAON(M)=YKAON(M)/FLOAT(NRUN)/DY
 1001 CONTINUE
       return
       end
