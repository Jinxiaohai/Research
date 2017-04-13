        subroutine inian1
        implicit double precision (a-h,o-z)
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        common /ana1/ ts(12)
        SAVE   
        if (ibstfg .ne. 0) then
           a = cosh(6d0)
           do 1001 i = 1, 12
              ts(i) = ts(i) * a
 1001      continue
        end if
        return
        end
