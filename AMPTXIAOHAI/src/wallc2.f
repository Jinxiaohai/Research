        subroutine wallc2(i, i1, i2, i3, t, tmin)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        x1p = gx(i)
        x2p = gy(i)
        x3p = gz(i)
        tf = ft(i)
        v1p = vx(i)
        v2p = vy(i)
        v3p = vz(i)
        if (v1p .gt. 0d0) then
           t1 = (5d0 * size1 - x1p) / v1p + tf
        else if (v1p .lt. 0d0) then
           t1 = (-5d0 * size1 - x1p) / v1p + tf
        else
           t1 = tlarge
        end if
        if (v2p .gt. 0d0) then
           t2 = (5d0 * size2 - x2p) / v2p + tf
        else if (v2p .lt. 0d0) then
           t2 = (- 5d0 * size2 - x2p) / v2p +tf
        else
           t2 = tlarge
        end if
        if (iconfg .eq. 5) then
           if (v3p .gt. 0d0) then
              t3 = (5d0 * size3 - x3p) / v3p + tf
           else if (v3p .lt. 0d0) then
              t3 = (- 5d0 * size3 - x3p) / v3p +tf
           else
              t3 = tlarge
           end if
        else
           t3 = tlarge
        end if
        tmin = min(t1, t2, t3)
        if (tmin .eq. t1) then
           if (v1p .gt. 0d0) then
              icsta(i) = 101
           else
              icsta(i) = 102
           end if
        end if
        if (tmin .eq. t2) then
           if (v2p .gt. 0d0) then
              icsta(i) = 103
           else
              icsta(i) = 104
           end if
        end if
        if (tmin .eq. t3) then
           if (v3p .gt. 0d0) then
              icsta(i) = 105
           else
              icsta(i) = 106
           end if
        end if
        return
        end
