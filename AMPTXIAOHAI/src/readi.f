        subroutine readi
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        double precision field(9)
        common /para1/ mul
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        SAVE   
        do 1001 i = 1, MAXPTN
           if (ievt .ne. 1 .and. i .eq. 1) then
              ityp0(i) = ntyp
              gx0(1) = field(1)
              gy0(1) = field(2)
              gz0(1) = field(3)
              ft0(1) = field(4)
              px0(1) = field(5)
              py0(1) = field(6)
              pz0(1) = field(7)
              e0(1) = field(8)
              xmass0(i) = field(9)
              mul = 1
           else
 900              read (27, *, end = 1000) neve, ntyp, field
              if (neve .lt. nsevt) goto 900
              if (neve .gt.
     &           nsevt + ievt - 1) goto 1000
              ityp0(i) = ntyp
              gx0(i) = field(1)
              gy0(i) = field(2)
              gz0(i) = field(3)
              ft0(i) = field(4)
              px0(i) = field(5)
              py0(i) = field(6)
              pz0(i) = field(7)
              e0(i) = field(8)
              xmass0(i) = field(9)
              mul = mul + 1
           end if
 1001   continue
 1000        continue
        return
        end
