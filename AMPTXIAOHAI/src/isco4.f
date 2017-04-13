        subroutine isco4(i, j, allok, tm, t1, t2)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
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
        logical allok
        allok = last(i) .ne. j .or. last(j) .ne. i
        icels1 = icels(i)
        ii1 = icels1 / 10000
        jj1 = (icels1 - ii1 * 10000) / 100
        kk1 = icels1 - ii1 * 10000 - jj1 * 100
        icels2 = icels(j)
        ii2 = icels2 / 10000
        jj2 = (icels2 - ii2 * 10000) / 100
        kk2 = icels2 - ii2 * 10000 - jj2 * 100
        i1 = i
        i2 = j
        p4 = ft(i2) - ft(i1)
        p1 = gx(i2) - gx(i1)
        p2 = gy(i2) - gy(i1)
        p3 = gz(i2) - gz(i1)
        if (ii1 - ii2 .gt. 5) then
           p1 = p1 + 10d0 * size1
        else if (ii1 - ii2 .lt. -5) then
           p1 = p1 - 10d0 * size1
        end if
        if (jj1 - jj2 .gt. 5) then
           p2 = p2 + 10d0 * size2
        else if (jj1 - jj2 .lt. -5) then
           p2 = p2 - 10d0 * size2
        end if
        if (kk1 - kk2 .gt. 5) then
           p3 = p3 + 10d0 * size3
        else if (kk1 - kk2 .lt. -5) then
           p3 = p3 - 10d0 * size3
        end if
        q4 = e(i1)
        q1 = px(i1)
        q2 = py(i1)
        q3 = pz(i1)
        r4 = e(i2)
        r1 = px(i2)
        r2 = py(i2)
        r3 = pz(i2)
        a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3
        b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3
        c = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3
        d = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3
        ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3
        f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3
        h = a + b
        if (h .gt. 0d0) then
           g = a
           a = -b
           b = -g
           g = c
           c = d
           d = g
           i1 = j
           i2 = i
        end if
        if (allok) then
           vp = a * d - b * ee
           allok = allok .and. vp .lt. 0d0
        end if
         if (allok) then
           dm2 = - f - (a ** 2 * d + b ** 2 * c - 2d0 * a * b * ee) /
     &           (ee ** 2 - c * d)
           allok = allok .and. dm2 .lt. cutof2
        end if
        if (allok) then
           tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (ee ** 2 - c * d)
           tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (ee ** 2 - c * d)
           tm = 0.5d0 * (tc1 + tc2)
           allok = allok .and. tm .gt. ft(i) .and. tm .gt. ft(j)
        end if
        if (allok) then
           rts2 = (q4 + r4) ** 2 - (q1 + r1) ** 2
     &          - (q2 + r2) ** 2 - (q3 + r3) ** 2
           allok = allok .and. rts2 .gt. rscut2
        end if
        if (.not. allok) then
           tm = tlarge
           t1 = tlarge
           t2 = tlarge
        else if (h .gt. 0d0) then
           t1 = tm
           t2 = tm
        else
           t1 = tm
           t2 = tm
        end if
        return
        end
