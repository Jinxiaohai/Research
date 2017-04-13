        subroutine dchin2(l, ii, i1, i2, i3, t)
        implicit double precision (a-h, o-z)
        SAVE   
        if (ii .eq. 1) then
           i = i1 - 2
           if (i .le. 0) i = i + 10
           ia = i
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k = i3 - 1, i3 + 1
                 ib = j
                 ic = k
                 if (j .eq. 0) ib = 10
                 if (j .eq. 11) ib = 1
                 if (k .ge. 1 .and. k .le. 10) then
                    call dchcel(l, ia, ib, ic, t)
                 end if
 1001         continue
 1002      continue
        end if
        if (ii .eq. 2) then
           i = i1 + 2
           if (i .ge. 11) i = i - 10
           ia = i
           do 1004 j = i2 - 1, i2 + 1
              do 1003 k = i3 - 1, i3 + 1
                 ib = j
                 ic = k
                 if (j .eq. 0) ib = 10
                 if (j .eq. 11) ib = 1
                 if (k .ge. 1 .and. k .le. 10) then
                    call dchcel(l, ia, ib, ic, t)
                 end if
 1003         continue
 1004      continue
        end if
        if (ii .eq. 3) then
           j = i2 - 2
           if (j .le. 0) j = j + 10
           ib = j
           do 1006 i = i1 - 1, i1 + 1
              do 1005 k = i3 - 1, i3 + 1
                 ia = i
                 ic = k
                 if (i .eq. 0) ia = 10
                 if (i .eq. 11) ia = 1
                 if (k .ge. 1 .and. k .le. 10) then
                    call dchcel(l, ia, ib, ic, t)
                 end if
 1005         continue
 1006      continue
        end if
        if (ii .eq. 4) then
           j = i2 + 2
           if (j .ge. 11) j = j - 10
           ib = j
           do 1008 i = i1 - 1, i1 + 1
              do 1007 k = i3 - 1, i3 + 1
                 ia = i
                 ic = k
                 if (i .eq. 0) ia = 10
                 if (i .eq. 11) ia = 1
                 if (k .ge. 1 .and. k .le. 10) then
                    call dchcel(l, ia, ib, ic, t)
                 end if
 1007         continue
 1008      continue
        end if
        if (ii .eq. 5) then
           if (i3 .eq. 2) goto 100
           k = i3 - 2
           ic = k
           do 1010 i = i1 - 1, i1 + 1
              do 1009 j = i2 - 1, i2 + 1
                 ia = i
                 ib = j
                 if (i .eq. 0) ia = 10
                 if (i .eq. 11) ia = 1
                 if (j .eq. 0) ib = 10
                 if (j .eq. 11) ib = 1
                     call dchcel(l, ia, ib, ic, t)
 1009         continue
 1010      continue
        end if
        if (ii .eq. 6) then
           if (i3 .eq. 9) goto 100
           k = i3 + 2
           ic = k
           do 1012 i = i1 - 1, i1 + 1
              do 1011 j = i2 - 1, i2 + 1
                 ia = i
                 ib = j
                 if (i .eq. 0) ia = 10
                 if (i .eq. 11) ia = 1
                 if (j .eq. 0) ib = 10
                 if (j .eq. 11) ib = 1
                     call dchcel(l, ia, ib, ic, t)
 1011         continue
 1012      continue
        end if
 100        continue
        return
        end
