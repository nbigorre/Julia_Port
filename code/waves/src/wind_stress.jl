function wind_stress(udif,vdif,step) 
#  !     ---------------------------------------------                     
#  USE header
#
#  implicit none 
#  integer i,j,k,m,step
#                                                                        
#  REAL(kind=rc_kind) :: udif(NI,NJ,NK),vdif(NI,NJ,NK)
#
#  REAL(kind=rc_kind) :: Kdudzt,Kdvdzt,fact,fac,rhoinv
#
#  INTEGER :: yw, ymid, yset
#
#  real :: stressmax, yyy

  udif .= 0e0
  vdif .= 0e0


#!*********************************
#! COMPUTATION OF THE WIND STRESS

  local stress_top_x=0e0
  local stress_top_y=0e0 
  local stressmax= 0.26e0  
  
  local yw= 50e0
  local ymid= 0.5* rc_kind(NJ)
  local yset =30e0

  for j in 1:NJ/2
     yyy= rc_kind(j)-0.5
     stress_top_y[1,j] = stressmax* 0.5* (1e0 +tanh((( yyy-yset)/yw)*PI ))
  end

  yset= NJ-30
  for j in NJ/2+1:NJ
    yyy= rc_kind(j)-0.5
    stress_top_y[1,j] = stressmax*0.5* (1e0 +tanh((-( yyy-yset)/yw)*PI ))
  end

   for j in 1:NJ
     if (abs(stress_top_y[1,j]-stressmax) < 0.001) 
      stress_top_y[1,j]=stressmax
     end
     if (abs(stress_top_y[1,j]) < 0.001) 
      stress_top_y[1,j]=0e0
     end
   end

   for i in 1:NI
   for j in 1:NJ
     stress_top_y[i,j]=stress_top_y[1,j]
   end
   end

   if (step>240) 
    stress_top_y=0e0
   end

  #************************************
# COMPUTATION OF THE SOURCE TERM


  fac= 1e0/(UL*DL*delta) 
  fact = DL/UL 

  for j in 1:NJ 
  for i in 1:NI 
      stress_top[i,j] = sqrt(stress_top_x[i,j]*stress_top_x[i,j]+ stress_top_y[i,j]*stress_top_y[i,j])
      #stress_top[i,j]=0
      rhoinv = 1e0/rho[i,j,NK]
      #rhoinv = 1e0/R0 
      Kdudzt= stress_top_x[i,j]*rhoinv*fact 
      Kdvdzt= stress_top_y[i,j]*rhoinv*fact 

      udif[i,j,NK]= fac*Jac[i,j,NK]*wz[i,j,NK]*Kdudzt   
      vdif[i,j,NK]= fac*Jac[i,j,NK]*wz[i,j,NK]*Kdvdzt   
    
    end
  end
                                                                        
end                                                         
                                                                        
