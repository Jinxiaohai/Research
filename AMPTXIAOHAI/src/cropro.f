      subroutine cropro(vx1, vy1, vz1, vx2, vy2, vz2)
      implicit double precision (a-h, o-z)
      common/cprod/ vx3, vy3, vz3
      SAVE   
      vx3 = vy1 * vz2 - vz1 * vy2
      vy3 = vz1 * vx2 - vx1 * vz2
      vz3 = vx1 * vy2 - vy1 * vx2
      return
      end
