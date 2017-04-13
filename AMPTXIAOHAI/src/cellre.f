        subroutine cellre(i, t)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
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
        common /ilist2/ icell, icel(10,10,10)
        common /ilist3/ size1, size2, size3, v1, v2, v3, size
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        logical good
        external integ
        t0 = t
 1000        continue
        if (iconfg .eq. 3 .or. iconfg .eq. 5) then
           k = mod(icsta(i), 10)
           if (k .eq. 1) then
              gx(i) = gx(i) - 10d0 * size1
              dgxa(i) = dgxa(i) + 10d0 * size1
              do 1001 ii = 1, ichkpt
                 if (next(ii) .eq. i) then
                    dgxa(ii) = dgxa(ii) - 10d0 * size1
                 end if
 1001         continue
           end if
           if (k .eq. 2) then
              gx(i) = gx(i) + 10d0 * size1
              dgxa(i) = dgxa(i) - 10d0 * size1
              do 1002 ii = 1, ichkpt
                 if (next(ii) .eq. i) then
                    dgxa(ii) = dgxa(ii) + 10d0 * size1
                 end if
 1002         continue
           end if
           if (k .eq. 3) then
              gy(i) = gy(i) - 10d0 * size2
              dgya(i) = dgya(i) + 10d0 * size2
              do 1003 ii = 1, ichkpt
                 if (next(ii) .eq. i) then
                    dgya(ii) = dgya(ii) - 10d0 * size2
                 end if
 1003         continue
           end if
           if (k .eq. 4) then
              gy(i) = gy(i) + 10d0 * size2
              dgya(i) = dgya(i) - 10d0 * size2
              do 1004 ii = 1, ichkpt
                 if (next(ii) .eq. i) then
                    dgya(ii) = dgya(ii) + 10d0 * size2
                 end if
 1004         continue
           end if
           if (iconfg .eq. 5) then
              if (k .eq. 5) then
                 gz(i) = gz(i) - 10d0 * size3
                 dgza(i) = dgza(i) + 10d0 * size3
                 do 1005 ii = 1, ichkpt
                    if (next(ii) .eq. i) then
                       dgza(ii) = dgza(ii) - 10d0 * size3
                    end if
 1005            continue
              end if
              if (k .eq. 6) then
                 gz(i) = gz(i) + 10d0 * size3
                 dgza(i) = dgza(i) - 10d0 * size3
                 do 1006 ii = 1, ichkpt
                    if (next(ii) .eq. i) then
                       dgza(ii) = dgza(ii) + 10d0 * size3
                    end if
 1006               continue
              end if
           end if
        else
           icels0 = icels(i)
           i1 = icels0 / 10000
           i2 = (icels0 - i1 * 10000) / 100
           i3 = icels0 - i1 * 10000 - i2 * 100
           if (i1 .ge. 1 .and. i1 .le. 10
     &        .and. i2 .ge. 1 .and. i2 .le. 10
     &        .and. i3 .ge. 1 .and. i3 .le. 10) then
              if (icel(i1, i2, i3) .eq. i) icel(i1, i2, i3) = nic(i)
              call oldcre(i)
              k = mod(icsta(i), 10)
              if (iconfg .eq. 1) then
                 good = (i1 .eq. 1 .and. k .eq. 2)
     &              .or. (i1 .eq. 10 .and. k .eq. 1)
     &              .or. (i2 .eq. 1 .and. k .eq. 4)
     &              .or. (i2 .eq. 10 .and. k .eq. 3)
     &              .or. (i3 .eq. 1 .and. k .eq. 6)
     &              .or. (i3 .eq. 10 .and. k .eq. 5)
              end if
              if (iconfg .eq. 2) then
                 good = (i3 .eq. 1 .and. k .eq. 6)
     &              .or. (i3 .eq. 10 .and. k .eq. 5)
              end if
              if (good) then
                 call newcre(i, icell)
                 icels(i) = 111111
              else
                 if (k .eq. 1) i1 = i1 + 1
                 if (k .eq. 2) i1 = i1 - 1
                 if (k .eq. 3) i2 = i2 + 1
                 if (k .eq. 4) i2 = i2 - 1
                 if (k .eq. 5) i3 = i3 + 1
                 if (k .eq. 6) i3 = i3 - 1
                 if (iconfg .eq. 2 .or. iconfg .eq. 4) then
                    if (i1 .eq. 0) then
                       i1 = 10
                       gx(i) = gx(i) + 10d0 * size1
                    end if
                    if (i1 .eq. 11) then
                       i1 = 1
                       gx(i) = gx(i) - 10d0 * size1
                    end if
                    if (i2 .eq. 0) then
                       i2 = 10
                       gy(i) = gy(i) + 10d0 * size2
                    end if
                    if (i2 .eq. 11) then
                       i2 = 1
                       gy(i) = gy(i) - 10d0 * size2
                    end if
                    if (iconfg .eq. 4) then
                       if (i3 .eq. 0) then
                          i3 = 10
                          gz(i) = gz(i) + 10d0 * size3
                       end if
                       if (i3 .eq. 11) then
                          i3 = 1
                          gz(i) = gz(i) - 10d0 * size3
                       end if
                    end if
                 end if
                 j = icel(i1, i2, i3)
                 call newcre(i, j)
                 icel(i1 ,i2, i3) = j
                 icels(i) = i1 * 10000 + i2 * 100 + i3
              end if
           else
              if (icell .eq. i) icell = nic(i)
              call oldcre(i)
              k = mod(icsta(i), 10)
              ddt = t - ft(i)
              dtt = t - size
              if (dtt .le. 0d0) then
                 i1 = integ((gx(i) + vx(i) * ddt) / size1) + 6
                 i2 = integ((gy(i) + vy(i) * ddt) / size2) + 6
                 i3 = integ((gz(i) + vz(i) * ddt) / size3) + 6
              else
                 i1 = integ((gx(i) + vx(i) * ddt) / 
     &               (size1 + v1 * dtt)) + 6
                 i2 = integ((gy(i) + vy(i) * ddt) /
     &               (size2 + v2 * dtt)) + 6
                 i3 = integ((gz(i) + vz(i) * ddt) /
     &               (size3 + v3 * dtt)) + 6
              end if 
              if (k .eq. 1) i1 = 1
              if (k .eq. 2) i1 = 10
              if (k .eq. 3) i2 = 1
              if (k .eq. 4) i2 = 10
              if (k .eq. 5) i3 = 1
              if (k .eq. 6) i3 = 10
              j = icel(i1, i2, i3)
              call newcre(i, j)
              icel(i1, i2, i3) = j
              icels(i) = i1 * 10000 + i2 * 100 + i3
           end if
        end if
        if (next(i) .ne. 0) then
           otmp = ot(next(i))
           ctmp = ct(next(i))
        end if
        if (i1 .eq. 11 .and. i2 .eq. 11 .and. i3 .eq. 11) then
           call dchout(i, k, t)
        else
           if (iconfg .eq. 1) then
              call dchin1(i, k, i1, i2, i3, t)
           else if (iconfg .eq. 2) then
              call dchin2(i, k, i1, i2, i3, t)
           else if (iconfg .eq. 4) then
              call dchin3(i, k, i1, i2, i3, t)              
           end if
        end if
        if (icsta(i) / 10 .eq. 11) then
           ot(next(i)) = otmp
           ct(next(i)) = ctmp
           next(next(i)) = i
           call wallc(i, i1, i2, i3, t0, tmin1)
           if (tmin1 .lt. ct(i)) then
              icsta(i) = icsta(i) + 10
              t0 = tmin1
              goto 1000
           end if
        end if
        return
        end
