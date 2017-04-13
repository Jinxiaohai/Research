      SUBROUTINE exchge(isg,ipi,jsg,ipj)
      implicit double precision  (a-h, o-z)
      PARAMETER (MAXSTR=150001)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      SAVE   
      k1=K1SGS(isg,ipi)
      k2=K2SGS(isg,ipi)
      px=PXSGS(isg,ipi)
      py=PYSGS(isg,ipi)
      pz=PZSGS(isg,ipi)
      pe=PESGS(isg,ipi)
      pm=PMSGS(isg,ipi)
      gx=GXSGS(isg,ipi)
      gy=GYSGS(isg,ipi)
      gz=GZSGS(isg,ipi)
      ft=FTSGS(isg,ipi)
      K1SGS(isg,ipi)=K1SGS(jsg,ipj)
      K2SGS(isg,ipi)=K2SGS(jsg,ipj)
      PXSGS(isg,ipi)=PXSGS(jsg,ipj)
      PYSGS(isg,ipi)=PYSGS(jsg,ipj)
      PZSGS(isg,ipi)=PZSGS(jsg,ipj)
      PESGS(isg,ipi)=PESGS(jsg,ipj)
      PMSGS(isg,ipi)=PMSGS(jsg,ipj)
      GXSGS(isg,ipi)=GXSGS(jsg,ipj)
      GYSGS(isg,ipi)=GYSGS(jsg,ipj)
      GZSGS(isg,ipi)=GZSGS(jsg,ipj)
      FTSGS(isg,ipi)=FTSGS(jsg,ipj)
      K1SGS(jsg,ipj)=k1
      K2SGS(jsg,ipj)=k2
      PXSGS(jsg,ipj)=px
      PYSGS(jsg,ipj)=py
      PZSGS(jsg,ipj)=pz
      PESGS(jsg,ipj)=pe
      PMSGS(jsg,ipj)=pm
      GXSGS(jsg,ipj)=gx
      GYSGS(jsg,ipj)=gy
      GZSGS(jsg,ipj)=gz
      FTSGS(jsg,ipj)=ft
      RETURN
      END
