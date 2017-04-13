      SUBROUTINE GRADUK(IX,IY,IZ,GRADXk,GRADYk,GRADZk)
      PARAMETER         (MAXX =    20,  MAXZ =   24)
      PARAMETER         (RHO0 = 0.168)
      COMMON  /DD/      RHO(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHOP(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ),
     &                  RHON(-MAXX:MAXX,-MAXX:MAXX,-MAXZ:MAXZ)
      common  /ss/      inout(20)
      SAVE   
      RXPLUS   = RHO(IX+1,IY,  IZ  ) 
      RXMINS   = RHO(IX-1,IY,  IZ  ) 
      RYPLUS   = RHO(IX,  IY+1,IZ  ) 
      RYMINS   = RHO(IX,  IY-1,IZ  ) 
      RZPLUS   = RHO(IX,  IY,  IZ+1) 
      RZMINS   = RHO(IX,  IY,  IZ-1) 
           GRADXk  = (RXPLUS - RXMINS)/2. 
           GRADYk  = (RYPLUS - RYMINS)/2.
           GRADZk  = (RZPLUS - RZMINS)/2.
           RETURN
           END
