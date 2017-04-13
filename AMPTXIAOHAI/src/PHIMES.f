      SUBROUTINE PHIMES(I1, I2, SRT, XSK1, XSK2, XSK3, XSK4, XSK5,
     1     XSK6, XSK7, SIGPHI)
        PARAMETER (MAXSTR=150001,MAXR=1,AMN=0.939457,
     1  AMP=0.93828,AP1=0.13496,
     2  AP2=0.13957,AM0=1.232,PI=3.1415926,CUTOFF=1.8966,AVMASS=0.9383)
        PARAMETER  (AKA=0.498, AKS=0.895, AOMEGA=0.7819,
     3               ARHO=0.77, APHI=1.02)
        parameter     (MX=4,MY=4,MZ=8,MPX=4,MPY=4,mpz=10,mpzp=10)
        PARAMETER  (MAXX=20,  MAXZ=24)
        COMMON /AA/ R(3,MAXSTR)
        COMMON /BB/ P(3,MAXSTR)
        COMMON /CC/ E(MAXSTR)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
        COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      SAVE   
        S = SRT ** 2
       SIGPHI = 1.E-08
        XSK1 = 0.0
        XSK2 = 0.0
        XSK3 = 0.0
        XSK4 = 0.0
        XSK5 = 0.0
        XSK6 = 0.0
        XSK7 = 0.0
         em1 = E(i1)
         em2 = E(i2)
         LB1 = LB(i1)
         LB2 = LB(i2)
         akap = aka
         XSK1 = 5.0
         scheck=(S-(em1+em2)**2)*(S-(em1-em2)**2)
         if(scheck.le.0) then
            write(99,*) 'scheck48: ', scheck
            stop
         endif
         pii=sqrt(scheck)
       if( lb1.eq.23.or.lb2.eq.23 .or. lb1.eq.21.or.lb2.eq.21 )then
          if(srt .gt. (ap1+akap))then
           pff = sqrt((S-(ap1+akap)**2)*(S-(ap1-akap)**2))
           XSK2 = 195.639*pff/pii/32./pi/S 
          endif
          if(srt .gt. (arho+akap))then
           pff = sqrt((S-(arho+akap)**2)*(S-(arho-akap)**2))
           XSK3 = 526.702*pff/pii/32./pi/S 
          endif
          if(srt .gt. (aomega+akap))then
           pff = sqrt((S-(aomega+akap)**2)*(S-(aomega-akap)**2))
           XSK4 = 355.429*pff/pii/32./pi/S 
          endif
          if(srt .gt. (ap1+aks))then
           pff = sqrt((S-(ap1+aks)**2)*(S-(ap1-aks)**2))
           XSK5 = 2047.042*pff/pii/32./pi/S 
          endif
          if(srt .gt. (arho+aks))then
           pff = sqrt((S-(arho+aks)**2)*(S-(arho-aks)**2))
           XSK6 = 1371.257*pff/pii/32./pi/S 
          endif
          if(srt .gt. (aomega+aks))then
           pff = sqrt((S-(aomega+aks)**2)*(S-(aomega-aks)**2))
           XSK7 = 482.292*pff/pii/32./pi/S 
          endif
       elseif( iabs(lb1).eq.30.or.iabs(lb2).eq.30 )then
          if(srt .gt. (ap1+akap))then
           pff = sqrt((S-(ap1+akap)**2)*(S-(ap1-akap)**2))
           XSK2 = 372.378*pff/pii/32./pi/S 
          endif
          if(srt .gt. (arho+akap))then
           pff = sqrt((S-(arho+akap)**2)*(S-(arho-akap)**2))
           XSK3 = 1313.960*pff/pii/32./pi/S 
          endif
          if(srt .gt. (aomega+akap))then
           pff = sqrt((S-(aomega+akap)**2)*(S-(aomega-akap)**2))
           XSK4 = 440.558*pff/pii/32./pi/S 
          endif
          if(srt .gt. (ap1+aks))then
           pff = sqrt((S-(ap1+aks)**2)*(S-(ap1-aks)**2))
           XSK5 = 1496.692*pff/pii/32./pi/S 
          endif
          if(srt .gt. (arho+aks))then
           pff = sqrt((S-(arho+aks)**2)*(S-(arho-aks)**2))
           XSK6 = 6999.840*pff/pii/32./pi/S 
          endif
          if(srt .gt. (aomega+aks))then
           pff = sqrt((S-(aomega+aks)**2)*(S-(aomega-aks)**2))
           XSK7 = 1698.903*pff/pii/32./pi/S 
          endif
       else
           srr1 = em1+em2
         if(srt .gt. (akap+akap))then
          srrt = srt - srr1
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          XSK2 = 1.69/(srrt**0.141 - 0.407)
          else
          XSK2 = 3.74 + 0.008*srrt**1.9
          endif                 
         endif
         if(srt .gt. (akap+aks))then
          srr2 = akap+aks
          srr = amax1(srr1,srr2)
          srrt = srt - srr
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          XSK3 = 1.69/(srrt**0.141 - 0.407)
          else
          XSK3 = 3.74 + 0.008*srrt**1.9
          endif
         endif
         if(srt .gt. (aks+aks))then
          srr2 = aks+aks
          srr = amax1(srr1,srr2)
          srrt = srt - srr
          if(srrt .lt. 0.3 .and. srrt .gt. 0.01)then
          XSK4 = 1.69/(srrt**0.141 - 0.407)
          else
          XSK4 = 3.74 + 0.008*srrt**1.9
          endif
         endif
      endif
        SIGPHI = XSK1 + XSK2 + XSK3 + XSK4 + XSK5 + XSK6 + XSK7
       RETURN
       END
