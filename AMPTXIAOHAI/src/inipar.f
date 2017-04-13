        subroutine inipar
        implicit double precision (a-h,o-z)
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        common /para6/ centy
        SAVE   
        if (ibstfg .ne. 0) then
           centy = -6d0
        end if
        return
        end
