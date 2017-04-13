            SUBROUTINE CMS(I1,I2,PX1CM,PY1CM,PZ1CM,SRT)
            PARAMETER (MAXSTR=150001)
            double precision px1,py1,pz1,px2,py2,pz2,em1,em2,e1,e2,
     1      s,ETOTAL,P1BETA,TRANSF,dBETAX,dBETAY,dBETAZ,dGAMMA,scheck
            COMMON   /BB/  P(3,MAXSTR)
            COMMON   /CC/  E(MAXSTR)
            COMMON   /BG/  BETAX,BETAY,BETAZ,GAMMA
            SAVE   
            PX1=dble(P(1,I1))
            PY1=dble(P(2,I1))
            PZ1=dble(P(3,I1))
            PX2=dble(P(1,I2))
            PY2=dble(P(2,I2))
            PZ2=dble(P(3,I2))
            EM1=dble(E(I1))
            EM2=dble(E(I2))
            E1=dSQRT(EM1**2+PX1**2+PY1**2+PZ1**2)
            E2=dSQRT(EM2**2+PX2**2+PY2**2+PZ2**2)
            S=(E1+E2)**2-(PX1+PX2)**2-(PY1+PY2)**2-(PZ1+PZ2)**2
            IF(S.LE.0) S=0d0
            SRT=sngl(dSQRT(S))
            ETOTAL = E1 + E2
            dBETAX  = (PX1+PX2) / ETOTAL
            dBETAY  = (PY1+PY2) / ETOTAL
            dBETAZ  = (PZ1+PZ2) / ETOTAL
            scheck=1.d0-dBETAX**2-dBETAY**2-dBETAZ**2
            if(scheck.le.0d0) then
               write(99,*) 'scheck1: ', scheck
               stop
            endif
            dGAMMA=1.d0/dSQRT(scheck)
            P1BETA = PX1*dBETAX + PY1*dBETAY + PZ1 * dBETAZ
            TRANSF = dGAMMA * ( dGAMMA * P1BETA / (dGAMMA + 1d0) - E1 )
            PX1CM  = sngl(dBETAX * TRANSF + PX1)
            PY1CM  = sngl(dBETAY * TRANSF + PY1)
            PZ1CM  = sngl(dBETAZ * TRANSF + PZ1)
            BETAX  = sngl(dBETAX)
            BETAY  = sngl(dBETAY)
            BETAZ  = sngl(dBETAZ)
            GAMMA  = sngl(dGAMMA)
            RETURN
            END
