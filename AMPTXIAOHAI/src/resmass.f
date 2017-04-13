      FUNCTION resmass(kf)
      PARAMETER  (arho=0.775,aomega=0.783,aeta=0.548,aks=0.894,
     1     aphi=1.019,adelta=1.232)
      PARAMETER  (wrho=0.149,womega=0.00849,weta=1.30E-6,wks=0.0498,
     1     wphi=0.00426,wdelta=0.118)
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      SAVE   
      if(kf.eq.113.or.abs(kf).eq.213) then
         amass=arho
         wid=wrho
      elseif(kf.eq.221) then
         amass=aeta
         wid=weta
      elseif(kf.eq.223) then
         amass=aomega
         wid=womega
      elseif(abs(kf).eq.313.or.abs(kf).eq.323) then
         amass=aks
         wid=wks
      elseif(kf.eq.333) then
         amass=aphi
         wid=wphi
      elseif(abs(kf).eq.1114.or.abs(kf).eq.2114
     1        .or.abs(kf).eq.2214.or.abs(kf).eq.2224) then
         amass=adelta
         wid=wdelta
      endif
      dmin=amass-2*wid
      dmax=amass+2*wid
      if(amass.eq.adelta) dmin=1.078
      FM=1.
      NTRY1=0
 10   DM = RANART(NSEED) * (DMAX-DMIN) + DMIN
      NTRY1=NTRY1+1
      fmass=(amass*wid)**2/((DM**2-amass**2)**2+(amass*wid)**2)
      IF((RANART(NSEED) .GT. FMASS/FM).AND. (NTRY1.LE.10)) GOTO 10
      resmass=DM
      RETURN
      END
