        subroutine dchin1(l, ii, i1, i2, i3, t)
        implicit double precision (a-h, o-z)
        SAVE   
        itest = 0
        if (ii .eq. 1) then
           if (i1 .eq. 1) goto 100
           if (i1 .eq. 2) then
              if (i2 .ge. 2 .and. i2 .le. 9 .and. i3 .ge. 2 .and.
     &           i3 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           i = i1 - 2
           do 1002 j = i2 - 1, i2 + 1
              do 1001 k = i3 - 1, i3 + 1
                 if (j .ge. 1 .and. j .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10)
     &                    call dchcel(l, i, j, k, t)
 1001         continue
 1002      continue
        end if
        if (ii .eq. 2) then
           if (i1 .eq. 10) goto 100
           if (i1 .eq. 9) then
              if (i2 .ge. 2 .and. i2 .le. 9 .and. i3 .ge. 2 .and.
     &           i3 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           i = i1 + 2
           do 1004 j = i2 - 1, i2 + 1
              do 1003 k = i3 - 1, i3 + 1
                 if (j .ge. 1 .and. j .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10)
     &                    call dchcel(l, i, j, k, t)
 1003         continue
 1004      continue
        end if
        if (ii .eq. 3) then
           if (i2 .eq. 1) goto 100
           if (i2 .eq. 2) then
              if (i1 .ge. 2 .and. i1 .le. 9 .and. i3 .ge. 2 .and.
     &           i3 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           j = i2 - 2
           do 1006 i = i1 - 1, i1 + 1
              do 1005 k = i3 - 1, i3 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. k .ge. 1 .and.
     &              k .le. 10)
     &              call dchcel(l, i, j, k, t)
 1005         continue
 1006      continue
        end if
        if (ii .eq. 4) then
           if (i2 .eq. 10) goto 100
           if (i2 .eq. 9) then
              if (i1 .ge. 2 .and. i1 .le. 9 .and. i3 .ge. 2 .and.
     &           i3 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           j = i2 + 2
           do 1008 i = i1 - 1, i1 + 1
              do 1007 k = i3 - 1, i3 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. k .ge. 1 .and.
     &           k .le. 10)
     &                 call dchcel(l, i, j, k, t)
 1007         continue
 1008      continue
        end if
        if (ii .eq. 5) then
           if (i3 .eq. 1) goto 100
           if (i3 .eq. 2) then
              if (i1 .ge. 2 .and. i1 .le. 9 .and. i2 .ge. 2 .and.
     &           i2 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           k = i3 - 2
           do 1010 i = i1 - 1, i1 + 1
              do 1009 j = i2 - 1, i2 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. j .ge. 1 .and.
     &           j .le. 10)
     &                 call dchcel(l, i, j, k, t)
 1009         continue
 1010      continue
        end if
        if (ii .eq. 6) then
           if (i3 .eq. 10) goto 100
           if (i3 .eq. 9) then
              if (i1 .ge. 2 .and. i1 .le. 9 .and. i2 .ge. 2 .and.
     &           i2 .le. 9) then
                 i = 11
                 j = 11
                 k = 11
                 call dchcel(l, i, j, k, t)
              end if
              goto 100
           end if
           k = i3 + 2
           do 1012 i = i1 - 1, i1 + 1
              do 1011 j = i2 - 1, i2 + 1
                 if (i .ge. 1 .and. i .le. 10 .and. j .ge. 1 .and.
     &           j .le. 10)
     &                 call dchcel(l, i, j, k, t)
 1011         continue
 1012      continue
        end if
 100        continue
        return
        end
