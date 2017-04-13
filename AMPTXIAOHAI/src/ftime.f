        subroutine ftime
        implicit double precision (a-h, o-z)
        external ftime1
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para4/ iftflg, ireflg, igeflg, ibstfg
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        common /prec4/ vx(MAXPTN), vy(MAXPTN), vz(MAXPTN)
        common /ilist4/ ifmpt, ichkpt, indx(MAXPTN)
        common /ilist5/ ct(MAXPTN), ot(MAXPTN), tlarge
        common /par1/ formt
        common/anim/nevent,isoft,isflag,izpc
        common /rndm3/ iseedp
        SAVE   
        iseed=iseedp
        do 1001 i = 1, MAXPTN
           ct(i)=0d0
           ot(i)=0d0
 1001   continue
        tlarge=1000000.d0
        if (iftflg .eq. 0) then
           if(isoft.eq.3.or.isoft.eq.4.or.isoft.eq.5) then
              do 1002 i = 1, mul
                 if (ft0(i) .gt. tlarge) ft0(i) = tlarge
 1002         continue
              goto 150
           else
           do 1003 i = 1, MAXPTN
              ft0(i) = tlarge
 1003      continue
           do 1004 i = 1, mul
              xmt2 = px0(i) ** 2 + py0(i) ** 2 + xmp ** 2
              formt = xmt2 / e0(i)           
              ft0(i) = ftime1(iseed)
              if (ft0(i) .gt. tlarge) ft0(i) = tlarge
 1004      continue
        endif
        end if
 150        continue
        if (mul .gt. 1) then
           call index1(MAXPTN, mul, ft0, indx)
        else
           indx(1)=1
        end if
        return
        end
