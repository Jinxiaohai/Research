      subroutine lorenz(energy, px, py, pz, bex, bey, bez)
      implicit double precision (a-h, o-z)
      common /lor/ enenew, pxnew, pynew, pznew
      SAVE   
      beta2 = bex ** 2 + bey ** 2 + bez ** 2
      if (beta2 .eq. 0d0) then
         enenew = energy
         pxnew = px
         pynew = py
         pznew = pz
      else
         if (beta2 .gt. 0.999999999999999d0) then
            beta2 = 0.999999999999999d0
            print *,'beta2=0.999999999999999'
         end if
         gam = 1.d0 / dsqrt(1.d0 - beta2)
         enenew = gam * (energy - bex * px - bey * py - bez * pz)
         pxnew = - gam * bex * energy + (1.d0 
     &        + (gam - 1.d0) * bex ** 2 / beta2) * px
     &        + (gam - 1.d0) * bex * bey/beta2 * py
     &        + (gam - 1.d0) * bex * bez/beta2 * pz     
         pynew = - gam * bey * energy 
     &        + (gam - 1.d0) * bex * bey / beta2 * px
     &        + (1.d0 + (gam - 1.d0) * bey ** 2 / beta2) * py
     &        + (gam - 1.d0) * bey * bez / beta2 * pz         
         pznew = - gam * bez * energy
     &        +  (gam - 1.d0) * bex * bez / beta2 * px
     &        + (gam - 1.d0) * bey * bez / beta2 * py
     &        + (1.d0 + (gam - 1.d0) * bez ** 2 / beta2) * pz    
      endif
      return
      end
