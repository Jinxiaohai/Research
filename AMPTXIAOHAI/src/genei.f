        subroutine genei
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /para2/ xmp, xmu, alpha, rscut2, cutof2
        common /para3/ nsevt, nevnt, nsbrun, ievt, isbrun
        common /para5/ iconfg, iordsc
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        common /prec5/ eta(MAXPTN), rap(MAXPTN), tau(MAXPTN)
        common /lor/ enenew, pxnew, pynew, pznew
        common /rndm3/ iseedp
        SAVE   
        external ran1
        iseed=iseedp
        incmul = 4000
        temp = 0.5d0
        etamin = -5d0        
        etamax = 5d0
        r0 = 5d0
        tau0 = 0.1d0
        deta = etamax - etamin
        do 1001 i = mul + 1, mul + incmul
           ityp0(i) = 21
           xmass0(i) = xmp
           call energy(e, temp)
           call momntm(px, py, pz, e)
           e = dsqrt(e ** 2 + xmp ** 2)
           if (iconfg .le. 3) then
              eta(i) = etamin + deta * ran1(iseed)
              bex = 0d0
              bey = 0d0
              bez = -tanh(eta(i))
              call lorenz(e, px, py, pz, bex, bey, bez)
              px0(i) = pxnew
              py0(i) = pynew
              pz0(i) = pznew
              e0(i) = enenew
           else
              px0(i) = px
              py0(i) = py
              pz0(i) = pz
              e0(i) = e
           end if
 1001   continue
        do 1002 i = mul + 1, mul + incmul
           if (iconfg .le. 3) then
              gz0(i) = tau0 * sinh(eta(i))
              ft0(i) = tau0 * cosh(eta(i))
              if (iconfg .eq. 1) then
                 call posit1(x, y, r0)
                 gx0(i) = x + px0(i) * ft0(i)/e0(i)
                 gy0(i) = y + py0(i) * ft0(i)/e0(i)
              else if (iconfg .eq. 2 .or. iconfg .eq. 3) then
                 call posit2(x, y)
                 gx0(i) = x
                 gy0(i) = y
              end if
           else
              ft0(i) = 0d0
              call posit3(x, y, z)
              gx0(i) = x
              gy0(i) = y
              gz0(i) = z
           end if
 1002   continue
        mul = mul + incmul
            if (mul .ge. MAXPTN .or. mul .eq. 0) then
           print *, 'event',ievt,'has',mul,'number of gluon',
     &          'adjusting counting is necessary'
           stop 'adarr'
        end if
        return
        end
