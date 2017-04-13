      FUNCTION PYI3AU(BE,EPS,IREIM) 
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      IF(EPS.LT.1.) GA=0.5*(1.+SQRT(1.-EPS))    
      IF(EPS.LT.0.) THEN    
        F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-PYSPEN(GA/(GA+BE-1.),0.,1)+    
     &  PYSPEN((BE-GA)/BE,0.,1)-PYSPEN((BE-GA)/(BE-1.),0.,1)+   
     &  (LOG(BE)**2-LOG(BE-1.)**2)/2.+LOG(GA)*LOG((GA+BE-1.)/BE)+   
     &  LOG(GA-1.)*LOG((BE-1.)/(GA+BE-1.))  
        F3IM=0. 
      ELSEIF(EPS.LT.1.) THEN    
        F3RE=PYSPEN((GA-1.)/(GA+BE-1.),0.,1)-PYSPEN(GA/(GA+BE-1.),0.,1)+    
     &  PYSPEN(GA/(GA-BE),0.,1)-PYSPEN((GA-1.)/(GA-BE),0.,1)+   
     &  LOG(GA/(1.-GA))*LOG((GA+BE-1.)/(BE-GA)) 
        F3IM=-PARU(1)*LOG((GA+BE-1.)/(BE-GA))   
      ELSE  
        RSQ=EPS/(EPS-1.+(2.*BE-1.)**2)  
        RCTHE=RSQ*(1.-2.*BE/EPS)    
        RSTHE=SQRT(RSQ-RCTHE**2)    
        RCPHI=RSQ*(1.+2.*(BE-1.)/EPS)   
        RSPHI=SQRT(RSQ-RCPHI**2)    
        R=SQRT(RSQ) 
        THE=ACOS(RCTHE/R)   
        PHI=ACOS(RCPHI/R)   
        F3RE=PYSPEN(RCTHE,RSTHE,1)+PYSPEN(RCTHE,-RSTHE,1)-  
     &  PYSPEN(RCPHI,RSPHI,1)-PYSPEN(RCPHI,-RSPHI,1)+   
     &  (PHI-THE)*(PHI+THE-PARU(1)) 
        F3IM=PYSPEN(RCTHE,RSTHE,2)+PYSPEN(RCTHE,-RSTHE,2)-  
     &  PYSPEN(RCPHI,RSPHI,2)-PYSPEN(RCPHI,-RSPHI,2)    
      ENDIF 
      IF(IREIM.EQ.1) PYI3AU=2./(2.*BE-1.)*F3RE  
      IF(IREIM.EQ.2) PYI3AU=2./(2.*BE-1.)*F3IM  
      RETURN    
      END   
