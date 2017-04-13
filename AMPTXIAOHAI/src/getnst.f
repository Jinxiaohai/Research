      subroutine getnst(srt)
      parameter (pimass=0.140,pi=3.1415926)
      COMMON/ppbmas/niso(15),nstate,ppbm(15,2),thresh(15),weight(15)
      common/ppb1/ene,factr2(6),fsum,ppinnb,s,wtot
      common/ppmm/pprr,ppee,pppe,rpre,xopoe,rree
      SAVE   
      s=srt**2
      nstate=0
      wtot=0.
      if(srt.le.thresh(1)) return
      do 1001 i=1,15
         weight(i)=0.
         if(srt.gt.thresh(i)) nstate=i
 1001 continue
      do 1002 i=1,nstate
         pf2=(s-(ppbm(i,1)+ppbm(i,2))**2)
     1        *(s-(ppbm(i,1)-ppbm(i,2))**2)/4/s
         weight(i)=pf2*niso(i)
         wtot=wtot+weight(i)
 1002 continue
      ene=(srt/pimass)**3/(6.*pi**2)
      fsum=factr2(2)+factr2(3)*ene+factr2(4)*ene**2
     1     +factr2(5)*ene**3+factr2(6)*ene**4
      return
      END
