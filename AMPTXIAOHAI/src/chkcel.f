        subroutine chkcel(il, i1, i2, i3, t, tmin, nc)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para5/ iconfg, iordsc
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist2/ icell, icel(10, 10, 10)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        SAVE   
        if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           jj = ichkpt
           do 1001 j = 1, jj
              call ck(j, ick)
                            jud2=j
              if (ick .eq. 1) call ud2(jud2, il, t, tmin, nc)
 1001      continue
           return
        end if
        if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
           l = icell
        else
           l = icel(i1, i2, i3)
        end if
        if (l .eq. 0) then
           return
        end if
        j = nic(l)
        if (j .eq. 0) then
           call ck(l, ick)
           if (ick .eq. 1) call ud2(l, il, t, tmin, nc)
        else
           call ck(l, ick)
           if (ick .eq. 1) call ud2(l, il, t, tmin, nc)
           do 10 while(j .ne. l)
              call ck(j, ick)
              if (ick .eq. 1) call ud2(j, il, t, tmin, nc)
              j = nic(j)
 10           continue
        end if
        return
        end
