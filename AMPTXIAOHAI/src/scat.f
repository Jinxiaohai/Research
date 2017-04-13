        subroutine scat(t, iscat, jscat)
        implicit double precision (a-h, o-z)
        SAVE   
        call newpos(t, iscat)
        call newpos(t, jscat)
        call newmom(t)
        return
        end
