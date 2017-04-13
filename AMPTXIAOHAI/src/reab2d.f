      real function reab2d(i1,i2,srt)
      PARAMETER      (MAXSTR=150001,MAXR=1,PI=3.1415926)
      parameter      (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
      PARAMETER      (AKA=0.498,ALA=1.1157,ASA=1.1974)
      parameter      (amn=0.938,ap1=0.14,arho=0.77,aomega=0.782)
       parameter       (maxx=20,maxz=24)
      COMMON   /AA/  R(3,MAXSTR)
      COMMON   /BB/  P(3,MAXSTR)
      COMMON   /CC/  E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      COMMON  /EE/      ID(MAXSTR),LB(MAXSTR)
      SAVE   
       reab2d=0
       LB1=iabs(LB(I1))
       LB2=iabs(LB(I2))
       ed1=e(i1)       
       ed2=e(i2)       
       pin2=(srt/2.)**2-amn**2
       pout2=((srt**2+ed1**2-ed2**2)/(2.*srt))**2-ed1**2
       if(pout2.le.0)return
       xpro=x2pi(srt)
       factor=1/4.
       if((lb1.ge.10.and.lb1.le.13).and.
     &    (lb2.ge.10.and.lb2.le.13))factor=1.
       if((lb1.ge.6.and.lb1.le.9).and.
     &    (lb2.gt.10.and.lb2.le.13))factor=1/2.
       if((lb2.ge.6.and.lb2.le.9).and.
     &    (lb1.gt.10.and.lb1.le.13))factor=1/2.
       reab2d=factor*pin2/pout2*xpro
       return
       end
