      subroutine pbarfs(srt,npion,iseed)
       parameter (pimass=0.140,pi=3.1415926) 
       Dimension factor(6),pnpi(6) 
      COMMON/RNDF77/NSEED
      SAVE   
       factor(2)=1.
       factor(3)=1.17e-01
       factor(4)=3.27e-03
       factor(5)=3.58e-05
       factor(6)=1.93e-07
       ene=(srt/pimass)**3/(6.*pi**2)
       do 1001 n=2,6 
           pnpi(n)=ene**n*factor(n)
 1001   continue
       pmax=max(pnpi(2),pnpi(3),pnpi(4),pnpi(5),pnpi(6))
       ntry=0
 10    npion=2+int(5*RANART(NSEED))
       if(npion.gt.6) goto 10
       thisp=pnpi(npion)/pmax  
       ntry=ntry+1 
       if((thisp.lt.RANART(NSEED)).and.(ntry.le.20)) go to 10
       return
       END
