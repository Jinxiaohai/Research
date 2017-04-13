        subroutine celasn
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /para5/ iconfg, iordsc
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist2/ icell, icel(10,10,10)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        SAVE   
        external integ
        i = ifmpt
        tt = ft(i)
        td = tt - size
        if (iconfg .eq. 1 .and. (size1 .eq. 0d0 .or.
     &     size2 .eq. 0d0 .or. size3 .eq. 0d0)) then
           i1 = 11
           i2 = 11
           i3 = 11
        else if (iconfg .eq. 4 .or. td .le. 0d0) then
           i1 = integ(gx(i) / size1) + 6
           i2 = integ(gy(i) / size2) + 6
           i3 = integ(gz(i) / size3) + 6
           if (integ(gx(i) / size1) .eq. gx(i) / size1 .and. 
     &        vx(i) .lt. 0d0)
     &        i1 = i1 - 1
           if (integ(gy(i) / size2) .eq. gy(i) / size2 .and. 
     &        vy(i) .lt. 0d0)
     &        i2 = i2 - 1
           if (integ(gz(i) / size3) .eq. gz(i) / size3 .and. 
     &        vz(i) .lt. 0d0)
     &        i3 = i3 - 1
        else
           i1 = integ(gx(i) / (size1 + v1 * td)) + 6
           i2 = integ(gy(i) / (size2 + v2 * td)) + 6
           i3 = integ(gz(i) / (size3 + v3 * td)) + 6
           if (integ(gx(i) / (size1 + v1 * td)) .eq. gx(i) / 
     &        (size1 + v1 * td) .and. vx(i) .lt. (i1 - 6) * v1)
     &        i1 = i1 - 1
           if (integ(gy(i) / (size2 + v2 * td)) .eq. gy(i)/
     &        (size2 + v2 * td) .and. vy(i) .lt. (i2 - 6) * v2)
     &        i2 = i2 - 1
           if (integ(gz(i) / (size3 + v3 * td)) .eq. gz(i)/
     &        (size3 + v3 * td) .and. vz(i) .lt. (i3 - 6) * v3)
     &        i3 = i3 - 1
        end if
        if (i1 .le. 0 .or. i1 .ge. 11 .or. i2 .le. 0 .or.
     &     i2 .ge. 11 .or. i3 .le. 0 .or. i3 .ge. 11) then
           i1 = 11
           i2 = 11
           i3 = 11
        end if
        if (i1 .eq. 11) then
           j = icell
           call newcre(i, j)
           icell = j
           icels(i) = 111111
        else
           j = icel(i1, i2, i3)
           call newcre(i, j)
           icel(i1, i2, i3) = j
           icels(i) = i1 * 10000 + i2 * 100 + i3
        end if
        return
        end
