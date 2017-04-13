        subroutine isco11(i, j, allok, tm, t1, t2)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para5/ iconfg, iordsc
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /aurec1/ jxa, jya, jza
        common /aurec2/ dgxa(MAXPTN), dgya(MAXPTN), dgza(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        logical allok, allokp
        allok = last(i) .ne. j .or. last(j) .ne. i
        tm = tlarge
        if (allok) then
           do 1000 ii = - 1, 1
              do 2000 jj = - 1, 1
                 do 3000 kk = - 1, 1
                 allokp = .true.
                 i1 = i
                 i2 = j
                 p4 = ft(j) - ft(i)
                 p1 = gx(j) - gx(i)
                 p2 = gy(j) - gy(i)
                 p3 = gz(j) - gz(i)
                 p1 = p1 + ii * 10d0 * size1
                 p2 = p2 + jj * 10d0 * size2
                 p3 = p3 + kk * 10d0 * size3
                 q4 = e(i)
                 q1 = px(i)
                 q2 = py(i)
                 q3 = pz(i)
                 r4 = e(j)
                 r1 = px(j)
                 r2 = py(j)
                 r3 = pz(j)
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
                 if (allokp) then
                    vp = a * d - b * ee
                    allokp = allokp .and. vp .lt. 0d0
                 end if
                 if (allokp) then
           dm2 = - f - (a ** 2 * d + b ** 2 * c - 2d0 * a * b * ee) /
     &            (ee ** 2 - c * d)
                    allokp = allokp .and. dm2 .lt. cutof2
                 end if
                 if (allokp) then
           tc1 = ft(i1) - e(i1) * (a * d - b * ee) / (ee ** 2 - c * d)
           tc2 = ft(i2) + e(i2) * (b * c - a * ee) / (ee ** 2 - c * d)
                    if (iordsc .eq. 20) then
                       tmp = min(tc1, tc2)
                    else if (iordsc .eq. 21) then
                       tmp = 0.5d0 * (tc1 + tc2)
                    else
                       tmp = max(tc1, tc2)
                    end if
           allokp = allokp .and. tmp .gt. ft(i) .and. tmp .gt. ft(j)
                 end if
                 if (allokp .and. tmp .lt. tm) then
                    tm = tmp
                    jxa = ii
                    jya = jj
                    jza = kk
                    ha = h
                    tc1a = tc1
                    tc2a = tc2
                 end if
 3000                 continue
 2000              continue
 1000           continue
           if (tm .eq. tlarge) then
              allok = .false.
           end if
        end if
        if (allok) then
           q4 = e(i1)
           q1 = px(i1)
           q2 = py(i1)
           q3 = pz(i1)
           r4 = e(i2)
           r1 = px(i2)
           r2 = py(i2)
           r3 = pz(i2)
           rts2 = (q4 + r4) ** 2 - (q1 + r1) ** 2
     &          - (q2 + r2) ** 2 - (q3 + r3) ** 2
           allok = allok .and. rts2 .gt. rscut2
        end if
        if (.not. allok) then
           tm = tlarge
           t1 = tlarge
           t2 = tlarge
        else if (ha .gt. 0d0) then
           t1 = tc2a
           t2 = tc1a
        else
           t1 = tc1a
           t2 = tc2a
        end if
        return
        end
