        subroutine ck(l, ick)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        SAVE   
        ick = 1
        if (ictype .eq. 1) then
           if (l .eq. ifmpt) ick = 0
        else if (ictype .eq. 0 .or. ictype .eq. 3 .or. 
     &     ictype .eq. 4) then
           if (l .eq. iscat .or. l .eq. jscat) ick = 0
        else
           if (l .eq. iscat .or. l .eq. jscat .or.
     &         l .eq. ifmpt) ick = 0
        end if
        return
        end
