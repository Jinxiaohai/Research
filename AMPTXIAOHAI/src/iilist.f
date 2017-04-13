        subroutine iilist
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /ilist1/
     &     iscat, jscat, next(MAXPTN), last(MAXPTN),
     &     ictype, icsta(MAXPTN),
     &     nic(MAXPTN), icels(MAXPTN)
        common /ilist2/ icell, icel(10,10,10)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        common /ilist6/ t, iopern, icolln
        SAVE   
        iscat = MAXPTN
        jscat = MAXPTN
        do 1001 i = 1, mul
           next(i) = 0
           last(i) = 0
           icsta(i) = 0
           nic(i) = 0
           icels(i) = 0
 1001   continue
        icell = 0
        do 1004 i1 = 1, 10
           do 1003 i2 = 1, 10
              do 1002 i3 = 1, 10
                 icel(i1, i2, i3) = 0
 1002         continue
 1003      continue
 1004   continue
        ichkpt = 0
        ifmpt = 1
        do 1005 i = 1, mul
           ct(i) = tlarge
           ot(i) = tlarge
 1005   continue
        iopern = 0
        icolln = 0
        t = 0.d0
        return
        end
