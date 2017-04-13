      SUBROUTINE PYSTFE(KF,X,Q2,XPQ)    
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /LUDAT1/ 
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)    
      SAVE /LUDAT2/ 
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      SAVE /PYPARS/ 
      DIMENSION XPQ(-6:6),XFDFLM(9) 
      CHARACTER CHDFLM(9)*5,HEADER*40   
      DATA CHDFLM/'UPVAL','DOVAL','GLUON','QBAR ','UBAR ','SBAR ',  
     &'CBAR ','BBAR ','TBAR '/  
      DATA HEADER/'Tung evolution package has been invoked'/    
      DATA INIT/0/  
      IF(MSTP(51).GE.11.AND.MSTP(51).LE.13.AND.MSTP(52).LE.1) THEN  
        XDFLM=MAX(0.51E-4,X)    
        Q2DFLM=MAX(10.,MIN(1E8,Q2)) 
        IF(MSTP(52).EQ.0) Q2DFLM=10.    
        DO 100 J=1,9    
        IF(MSTP(52).EQ.1.AND.J.EQ.9) THEN   
          Q2DFLM=Q2DFLM*(40./PMAS(6,1))**2  
          Q2DFLM=MAX(10.,MIN(1E8,Q2))   
        ENDIF   
        XFDFLM(J)=0.    
  100   CONTINUE    
        IF(X.LT.0.51E-4.AND.ABS(PARP(51)-1.).GT.0.01) THEN  
          CXS=(0.51E-4/X)**(PARP(51)-1.)    
          DO 110 J=1,7  
  110     XFDFLM(J)=XFDFLM(J)*CXS   
        ENDIF   
        XPQ(0)=XFDFLM(3)    
        XPQ(1)=XFDFLM(2)+XFDFLM(5)  
        XPQ(2)=XFDFLM(1)+XFDFLM(5)  
        XPQ(3)=XFDFLM(6)    
        XPQ(4)=XFDFLM(7)    
        XPQ(5)=XFDFLM(8)    
        XPQ(6)=XFDFLM(9)    
        XPQ(-1)=XFDFLM(5)   
        XPQ(-2)=XFDFLM(5)   
        XPQ(-3)=XFDFLM(6)   
        XPQ(-4)=XFDFLM(7)   
        XPQ(-5)=XFDFLM(8)   
        XPQ(-6)=XFDFLM(9)   
      ELSE  
        IF(INIT.EQ.0) THEN  
          I1=0  
          IF(MSTP(52).EQ.4) I1=1    
          IHDRN=1   
          NU=MSTP(53)   
          I2=MSTP(51)   
          IF(MSTP(51).GE.11) I2=MSTP(51)-3  
          I3=0  
          IF(MSTP(52).EQ.3) I3=1    
          ALAM=0.75*PARP(1) 
          TPMS=PMAS(6,1)    
          QINI=PARP(52) 
          QMAX=PARP(53) 
          XMIN=PARP(54) 
          INIT=1    
        ENDIF   
        Q=SQRT(Q2)  
        DO 200 I=-6,6   
        FIXQ=0. 
  200   XPQ(I)=X*FIXQ   
        XPS=XPQ(1)  
        XPQ(1)=XPQ(2)   
        XPQ(2)=XPS  
        XPS=XPQ(-1) 
        XPQ(-1)=XPQ(-2) 
        XPQ(-2)=XPS 
      ENDIF 
      RETURN    
      END   
