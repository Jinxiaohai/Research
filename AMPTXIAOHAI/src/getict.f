        subroutine getict(t1)
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &       px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &       xmass(MAXPTN), ityp(MAXPTN)
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        SAVE   
        t1 = tlarge
        iscat = 0
        jscat = 0
        do 1001 i = 1, ichkpt
           if (ot(i) .lt. t1) then
              t1 = ot(i)
              iscat = i
           end if
 1001   continue
        if (iscat .ne. 0) jscat = next(iscat)
        if (iscat .ne. 0 .and. jscat .ne. 0) then
           if (icsta(iscat) .eq. 0 .and. icsta(jscat) .eq. 0) then
              ictype = 0
           else
              ictype = 4
           end if
        else if (iscat .ne. 0 .or. jscat .ne. 0) then
           ictype = 3
        end if
        if (ifmpt .le. mul) then
           if (ft(ifmpt) .lt. t1) then
              ictype = 1
              t1 = ft(ifmpt)
           else if (ft(ifmpt) .eq. t1) then
              if (ictype .eq. 0) ictype = 2
              if (ictype .eq. 3) ictype = 5
              if (ictype .eq. 4) ictype = 6
           end if
        end if
        return
        end
