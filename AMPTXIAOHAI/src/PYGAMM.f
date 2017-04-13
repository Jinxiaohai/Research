      FUNCTION PYGAMM(X)    
      DIMENSION B(8)    
      DATA B/-0.57719165,0.98820589,-0.89705694,0.91820686, 
     &-0.75670408,0.48219939,-0.19352782,0.03586834/    
      NX=INT(X) 
      DX=X-NX   
      PYGAMM=1. 
      DO 100 I=1,8  
  100 PYGAMM=PYGAMM+B(I)*DX**I  
      IF(X.LT.1.) THEN  
        PYGAMM=PYGAMM/X 
      ELSE  
        DO 110 IX=1,NX-1    
  110   PYGAMM=(X-IX)*PYGAMM    
      ENDIF 
      RETURN    
      END   
