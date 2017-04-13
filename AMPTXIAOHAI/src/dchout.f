        subroutine dchout(l, ii, t)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        external integ
        tt = ft(l)
        td = t - size
        x1 = gx(l) + vx(l) * (t - tt)
        x2 = gy(l) + vy(l) * (t - tt)
        x3 = gz(l) + vz(l) * (t - tt)
        if (td .le. 0d0) then
           i1 = integ(x1 / size1) + 6
           i2 = integ(x2 / size2 ) + 6
           i3 = integ(x3 / size3 ) + 6
           if (integ(x1 / size1) .eq. x1 / size1 .and. vx(l) .lt. 0d0)
     &        i1 = i1 - 1
           if (integ(x2 / size2) .eq. x2 / size2 .and. vy(l) .lt. 0d0)
     &        i2 = i2 - 1
           if (integ(x3 / size3) .eq. x3 / size3 .and. vz(l) .lt. 0d0)
     &        i3 = i3 - 1
        else
           i1 = integ(x1 / (size1 + v1 * td)) + 6
           i2 = integ(x2 / (size2 + v2 * td)) + 6
           i3 = integ(x3 / (size3 + v3 * td)) + 6
           if (integ(x1 / (size1 + v1 * td)) .eq. 
     &        x1 / (size1 +v1 * td) .and. 
     &        vx(l) .lt. (i1 - 6) * v1) i1 = i1 - 1
           if (integ(x2 / (size2 + v2 * td)) .eq.
     &        x2 / (size2 + v2 * td) .and.
     &        vy(l) .lt. (i2 - 6) * v2) i2 = i2 - 1
           if (integ(x3 / (size3 + v3 * td)) .eq. 
     &        x3 / (size3 + v3 * td) .and.
     &        vz(l) .lt. (i3 - 6) * v3) i3 = i3 - 1
        end if
        if (ii .eq. 1) then
           i = 9
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k = i3 - 1, i3 + 1
                 if (j .ge. 1 .and. j .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1001         continue
 1002      continue
        end if
        if (ii .eq. 2) then
           i = 2
           do 1004 j = i2 - 1, i2 + 1
              do 1003 k = i3 - 1, i3 + 1
                 if (j .ge. 1 .and. j .le. 10 .and. k .ge. 1 .and. 
     &              k .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1003         continue
 1004      continue
        end if
        if (ii .eq. 3) then
           j = 9
           do 1006 i = i1 - 1, i1 + 1
              do 1005 k = i3 - 1, i3 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1005         continue
 1006      continue
        end if
        if (ii .eq. 4) then
           j = 2
           do 1008 i = i1 - 1, i1 + 1
              do 1007 k = i3 - 1, i3 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1007         continue
 1008      continue
        end if
        if (ii .eq. 5) then
           k = 9
           do 1010 i = i1 - 1, i1 + 1
              do 1009 j = i2 - 1, i2 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. j .ge. 1 .and.
     &              j .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1009         continue
 1010      continue
        end if
        if (ii .eq. 6) then
           k = 2
           do 1012 i = i1 - 1, i1 + 1
              do 1011 j = i2 - 1, i2 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. j .ge. 1 .and.
     &              j .le. 10) then
                    call dchcel(l, i, j, k, t)
                 end if
 1011         continue
 1012      continue
        end if
        return
        end
