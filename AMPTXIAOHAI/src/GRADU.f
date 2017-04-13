      SUBROUTINE GRADU(IOPT,IX,IY,IZ,GRADX,GRADY,GRADZ)
      PARAMETER         (MAXX =    20,  MAXZ =   24)
      PARAMETER         (RHO0 = 0.167)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      common  /ss/      inout(20)
      common  /tt/  PEL(-maxx:maxx,-maxx:maxx,-maxz:maxz)
     &,rxy(-maxx:maxx,-maxx:maxx,-maxz:maxz)
      SAVE   
      RXPLUS   = RHO(IX+1,IY,  IZ  ) / RHO0
      RXMINS   = RHO(IX-1,IY,  IZ  ) / RHO0
      RYPLUS   = RHO(IX,  IY+1,IZ  ) / RHO0
      RYMINS   = RHO(IX,  IY-1,IZ  ) / RHO0
      RZPLUS   = RHO(IX,  IY,  IZ+1) / RHO0
      RZMINS   = RHO(IX,  IY,  IZ-1) / RHO0
      den0     = RHO(IX,  IY,  IZ) / RHO0
      ene0     = pel(IX,  IY,  IZ) 
      GOTO (1,2,3,4,5) IOPT
       if(iopt.eq.6)go to 6
       if(iopt.eq.7)go to 7
    1 CONTINUE
           GRADX  = -0.062 * (RXPLUS - RXMINS) + 0.03525 * (RXPLUS**2 -
     &                                                      RXMINS**2)
           GRADY  = -0.062 * (RYPLUS - RYMINS) + 0.03525 * (RYPLUS**2 -
     &                                                      RYMINS**2)
           GRADZ  = -0.062 * (RZPLUS - RZMINS) + 0.03525 * (RZPLUS**2 -
     &                                                      RZMINS**2)
           RETURN
    2 CONTINUE
           EXPNT = 1.3333333
           GRADX = -0.109 * (RXPLUS - RXMINS) 
     &     + 0.082 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -0.109 * (RYPLUS - RYMINS) 
     &     + 0.082 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -0.109 * (RZPLUS - RZMINS) 
     &     + 0.082 * (RZPLUS**EXPNT-RZMINS**EXPNT)
           RETURN
    3 CONTINUE
           EXPNT = 1.1666667
          acoef = 0.178
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
                 RETURN
    4   CONTINUE
       eh=4.
       eqgp=7.
           acoef=0.178
           EXPNT = 1.1666667
       denr=rho(ix,iy,iz)/rho0
       if(denr.le.eh.or.denr.ge.eqgp)then
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
       else
          acoef1=0.178
          acoef2=0.0
          expnt2=2./3.
           GRADX =-acoef1* (RXPLUS**EXPNT-RXMINS**EXPNT)
     &                 -acoef2* (RXPLUS**expnt2 - RXMINS**expnt2) 
           GRADy =-acoef1* (RyPLUS**EXPNT-RyMINS**EXPNT)
     &                 -acoef2* (RyPLUS**expnt2 - RyMINS**expnt2) 
           GRADz =-acoef1* (RzPLUS**EXPNT-RzMINS**EXPNT)
     &                 -acoef2* (RzPLUS**expnt2 - RzMINS**expnt2) 
       endif
       return
    5   CONTINUE
           EXPNT = 2.77
           GRADX = -0.0516 * (RXPLUS - RXMINS) 
     &     + 0.02498 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -0.0516 * (RYPLUS - RYMINS) 
     &     + 0.02498 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -0.0516 * (RZPLUS - RZMINS) 
     &     + 0.02498 * (RZPLUS**EXPNT-RZMINS**EXPNT)
           RETURN
    6 CONTINUE
       if(ene0.le.0.5)then
           GRADX  = -0.062 * (RXPLUS - RXMINS) + 0.03525 * (RXPLUS**2 -
     &                                                      RXMINS**2)
           GRADY  = -0.062 * (RYPLUS - RYMINS) + 0.03525 * (RYPLUS**2 -
     &                                                      RYMINS**2)
           GRADZ  = -0.062 * (RZPLUS - RZMINS) + 0.03525 * (RZPLUS**2 -
     &                                                      RZMINS**2)
           RETURN
       endif
       if(ene0.gt.0.5.and.ene0.le.1.5)then
       ef=36./1000.
           GRADX  = -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
           GRADy  = -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
           GRADz  = -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
       if(ene0.gt.1.5)then
       ef=36./1000.
       cf0=0.8
        GRADX  =0.5*cf0*(rxplus**0.333-rxmins**0.333) 
     &         -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
        GRADy  =0.5*cf0*(ryplus**0.333-rymins**0.333) 
     &         -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
        GRADz  =0.5*cf0*(rzplus**0.333-rzmins**0.333) 
     &         -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
    7 CONTINUE
       if(den0.le.4.5)then
           EXPNT = 1.1666667
          acoef = 0.178
           GRADX = -acoef * (RXPLUS - RXMINS) 
     &     + 0.1515 * (RXPLUS**EXPNT-RXMINS**EXPNT)
           GRADY = -acoef * (RYPLUS - RYMINS) 
     &     + 0.1515 * (RYPLUS**EXPNT-RYMINS**EXPNT)
           GRADZ = -acoef * (RZPLUS - RZMINS) 
     &     + 0.1515 * (RZPLUS**EXPNT-RZMINS**EXPNT)
       return
       endif
       if(den0.gt.4.5.and.den0.le.5.1)then
       ef=36./1000.
           GRADX  = -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
           GRADy  = -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
           GRADz  = -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
       if(den0.gt.5.1)then
       ef=36./1000.
       cf0=0.8
        GRADX  =0.5*cf0*(rxplus**0.333-rxmins**0.333) 
     &         -0.5*ef* (RXPLUS**0.67-RXMINS**0.67)
        GRADy  =0.5*cf0*(ryplus**0.333-rymins**0.333) 
     &         -0.5*ef* (RyPLUS**0.67-RyMINS**0.67)
        GRADz  =0.5*cf0*(rzplus**0.333-rzmins**0.333) 
     &         -0.5*ef* (RzPLUS**0.67-RzMINS**0.67)
           RETURN
       endif
        END
