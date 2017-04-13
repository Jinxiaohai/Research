      subroutine coulin(masspr,massta,NUM)
        integer  zta,zpr
        PARAMETER (MAXSTR=150001)
        common  /EE/ ID(MAXSTR),LB(MAXSTR)
        COMMON  /ZZ/ ZTA,ZPR
      SAVE   
        MASS=MASSTA+MASSPR
        DO 500 IRUN=1,NUM
        do 100 i = 1+(IRUN-1)*MASS,zta+(IRUN-1)*MASS
        LB(i) = 1
  100   continue
        do 200 i = zta+1+(IRUN-1)*MASS,massta+(IRUN-1)*MASS
        LB(i) = 2
  200   continue
        do 300 i = massta+1+(IRUN-1)*MASS,massta+zpr+(IRUN-1)*MASS
        LB(i) = 1
  300   continue
        do 400 i = massta+zpr+1+(IRUN-1)*MASS,
     1  massta+masspr+(IRUN-1)*MASS
        LB(i) = 2
  400   continue
  500   CONTINUE
        return
        end
