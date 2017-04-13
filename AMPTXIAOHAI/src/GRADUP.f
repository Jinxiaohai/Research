      SUBROUTINE GRADUP(IX,IY,IZ,GRADXP,GRADYP,GRADZP)
      PARAMETER         (MAXX =    20,  MAXZ =   24)
      PARAMETER         (RHO0 = 0.168)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                     RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      common  /ss/      inout(20)
      SAVE   
      RXPLUS   = RHOP(IX+1,IY,  IZ  ) / RHO0
      RXMINS   = RHOP(IX-1,IY,  IZ  ) / RHO0
      RYPLUS   = RHOP(IX,  IY+1,IZ  ) / RHO0
      RYMINS   = RHOP(IX,  IY-1,IZ  ) / RHO0
      RZPLUS   = RHOP(IX,  IY,  IZ+1) / RHO0
      RZMINS   = RHOP(IX,  IY,  IZ-1) / RHO0
           GRADXP  = (RXPLUS - RXMINS)/2. 
           GRADYP  = (RYPLUS - RYMINS)/2.
           GRADZP  = (RZPLUS - RZMINS)/2.
           RETURN
      END
