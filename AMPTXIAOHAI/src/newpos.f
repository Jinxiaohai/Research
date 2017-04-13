        subroutine newpos(t, i)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        dt1 = ct(i) - ft(i)
        gx(i) = gx(i) + vx(i) * dt1
        gy(i) = gy(i) + vy(i) * dt1
        gz(i) = gz(i) + vz(i) * dt1
        ft(i) = ct(i)
        if (iconfg .le. 3) then
           if (ft(i) .le. abs(gz(i))) then
              eta(i) = 1000000.d0
           else
              eta(i) = 0.5d0 * log((ft(i) + gz(i)) / (ft(i) - gz(i)))
           end if
           if(eta(i).lt.1000000.d0) then
              tau(i) = ft(i) / cosh(eta(i))
           else
              tau(i) = 1d-10
           endif
        end if
        return
        end
