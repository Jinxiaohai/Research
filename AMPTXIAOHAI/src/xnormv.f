      subroutine xnormv(vx, vy, vz)
      implicit double precision (a-h, o-z)
      SAVE   
      vv = dsqrt(vx ** 2 + vy ** 2 + vz ** 2)
      vx = vx / vv
      vy = vy / vv
      vz = vz / vv
      return
      end
