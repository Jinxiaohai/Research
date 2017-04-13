        subroutine boosti
        implicit double precision (a-h, o-z)
        parameter (MAXPTN=400001)
        common /para1/ mul
        common /para6/ centy
        COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &       PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &       XMASS0(MAXPTN), ITYP0(MAXPTN)
        common /lor/ enenew, pxnew, pynew, pznew
        SAVE   
        external lorenz
        bex = 0d0 
        bey = 0d0
        bez = - tanh(centy)
        do 1001 i = 1, mul
           px1 = gx0(i)
           py1 = gy0(i)
           pz1 = gz0(i)
           e1 = ft0(i)
           call lorenz(e1, px1, py1, pz1, bex, bey, bez)
           gx0(i) = pxnew
           gy0(i) = pynew
           gz0(i) = pznew
           ft0(i) = enenew
           px1 = px0(i)
           py1 = py0(i)
           pz1 = pz0(i)
           e1 = e0(i)
           call lorenz(e1, px1, py1, pz1, bex, bey, bez)
           px0(i) = pxnew
           py0(i) = pynew
           pz0(i) = pznew
           e0(i) = enenew
 1001   continue
        return
        end
