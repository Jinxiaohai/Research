      FUNCTION asinh(x)
      SAVE
      if(x.gt.0) then
         ASINH=alog(x+sqrt(x**2+1.))
      else
         ASINH=-alog(-x+sqrt(x**2+1.))
      endif
      return
      end
