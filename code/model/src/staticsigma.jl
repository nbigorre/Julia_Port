function staticsigma()

#  USE header
#  !     ------------------------                                          
#  !     FOR THE REGION BELOW THE TOPMOST LAYER                            
#  !     modified a for periodicew bc                                      
#  !     This subroutine updates all those quantities that are a function  
#  !     of time and  (time and depth). It is called  every time step.     
#  !     The surface bc  wfbct is updated here.                            
#!  implicit REAL(kind=rc_kind) :: (a-h,o-z) 
#  integer i,j,k 
#  REAL(kind=rc_kind) :: hpd,hpdinv,hu,hv,hx,hy,hxpdx,hypdy,z,temp,       &
#       be2,wxk,wyk,d2,dnk,dnkm1,sig,ep,epm1,dnkm1p,zpd,           &
#       g13(0:NI+1,0:NJ+1,0:NK+1),g23(0:NI+1,0:NJ+1,0:NK+1)          
#  !     Ddx and Ddy are known from init                                   
#  !                                                                       
#  !     pfac is defined in findzall     
  local g13 = OffsetArrays.zeros(rc_kind, 0:NI+1,0:NJ+1,0:NK+1)
  local g23 = OffsetArrays.zeros(rc_kind, 0:NI+1,0:NJ+1,0:NK+1)                                  
  local dnk= rc_kind(NK) 

  local dnkm1 = 0e0
@static if (cppdefs.fixed_bottom_thickness)
  dnkm1= rc_kind(NK-1-1) 
else
  dnkm1= rc_kind(NK-1) 
end

  local dnkm1p= dnkm1/pfac 
  local be2= beta*EPS*EPS 
  local epm1= exp(pfac) -1e0 
  local ep = pfac*exp(1e0)                                                                      
  for j in 0:NJ+1 
     for i in 0:NI+1 
        #     All these variables are functions of time                         
            local hpd= dztop +D[i,j] 
            local hpdinv= 1e0/hpd 
#          local hu = 0e0
#            if (i == 0)                                          
#               hu = HDL*(h[i+1,j]-h[NI-1,j])                           
#            elseif  (i == NI+1)                                 
#               hu= HDL*(h[2,j]-h[i-1,j])                              
#            else                                                      
#               hu= 0.5e0*HDL*(h[i+1,j]-h[i-1,j])                      
#            end 
#             local hv = 0e0                      
#            if (j == 0)                                          
#               hv= HDL*(h[i,j+1]-h[i,j])                              
#            elseif (j == NJ+1)                                  
#               hv= HDL*(h[i,j]-h[i,j-1])                              
#            else                                                      
#               hv= 0.5e0*HDL*(h[i,j+1]-h[i,j-1])                      
#            end                       
#            local hx= hu*ux[i,j] + hv*vx[i,j]                               
#            local hy= hu*uy[i,j] + hv*vy[i,j]                               
#            local hxpdx= hx + Ddx[i,j]                               
#            local hypdy= hy + Ddy[i,j]                               
#     wz is not a function of depth when the sigma lines are equally spa
#     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k
#     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)
#                                                                       
                                                                       
            for k in 0:NK-1 
               local z = zc[i,j,k] 
               local zpd= z +dztop 
               wz[i,j,k]= -epm1*dnkm1p/(epm1*zpd +hpd)   
               Jac[i,j,k]= J2d[i,j]/wz[i,j,k] 
#     All these variables are functions of time and depth               
#     z is  dimensionaless and is equal to  z*/DL                       
#               z= ((rc_kind(k)-0.5e0)/rc_kind(NK)) *hpd - D[i,j]        
#               local temp= ( z+D[i,j])*hpdinv                               
#               wt(i,j,k)= -temp*HDL*hdt[i,j]*rc_kind(NK)*hpdinv         
#     now hdt computed in vhydro already contains HDL                   
               local sig= rc_kind(k)-0.5e0 
#              wt is at cell faces                                      
#               wt[i,j,k]= -sig*hdt[i,j]*hpdinv                        
#               wt[i,j,k]= -temp*hdt[i,j]*rc_kind(NK)*hpdinv            
#               wx[i,j,k]= (Ddx[i,j] - temp*hxpdx)*rc_kind(NK)*hpdinv                                 
#               wy[i,j,k]= (Ddy[i,j] - temp*hypdy)* rc_kind(NK)*hpdinv                                 
               wx[i,j,k]= wz[i,j,k]*zpd*hpdinv*Ddx[i,j] 
               wy[i,j,k]= wz[i,j,k]*zpd*hpdinv*Ddy[i,j] 
               g13[i,j,k]= ux[i,j]*wx[i,j,k] +uy[i,j]*wy[i,j,k] 
               g23[i,j,k]= vx[i,j]*wx[i,j,k] +vy[i,j]*wy[i,j,k] 
#               g33[i,j,k]= wx[i,j,k]*wx[i,j,k] +wy[i,j,k]*wy[i,j,k] + be2*wz[i,j]*wz[i,j]                               
            end 
            for k in 0:NK-1 
               local sig= rc_kind(k) 
               wt[i,j,k]= 0e0 
               local z= zc[i,j,k] 
               local zpd= z +dztop 
               wzk[i,j,k] = -epm1*dnkm1p/(epm1*zpd +hpd) 
            end
            wzk[i,j,NK-1]= 0.5*(wz[i,j,NK] + wz[i,j,NK-1]) 
#     temporary fix ****** for k=0  - use linear extrapolation          
#            wx[i,j,0]= 2e0*wx[i,j,1] -wx[i,j,2]                      
#            wz[i,j,0]= 2e0*wz[i,j,1] -wz[i,j,2]                      
#            wy[i,j,0]= 2e0*wy[i,j,1] -wy[i,j,2]                      
#            Jac[i,j,0]= 2e0*Jac[i,j,1] -Jac[i,j,2]                   
#            g13[i,j,0]= 2e0*g13[i,j,1] -g13[i,j,2]                   
#            g23[i,j,0]= 2e0*g23[i,j,1] -g23[i,j,2]                   
#                                                                       
#     We evaluate gk at i=0,NI+1,j=0,NJ+1 only because these are used   
#     at k=0,NK to fill the horizontal edges in mgpfill or hnbc.        
            for k in 0:NK 
#               local z= (rc_kind(k)/rc_kind(NK))*hpd - D[i,j]                 
               local z= zf[i,j,k] 
               local zpd= z +dztop 
               local sig= rc_kind(k) 
#               local temp= ( z+D[i,j])*hpdinv                               
#              local wtk= -temp*HDL*hdt[i,j]*rc_kind(NK)*hpdinv               
#    now hdt computed in vhydro already contains HDL                   
#                wtk= -temp*hdt[i,j]*rc_kind(NK)*hpdinv                  
#                wxk= (Ddx[i,j] - temp*hxpdx)*rc_kind(NK)*hpdinv                                 
#                wyk= (Ddy[i,j] - temp*hypdy)*rc_kind(NK)*hpdinv                                 
#                wtk= -sig*hdt[i,j]*hpdinv                              
                                                                        
               local temp= (epm1*dnkm1p/(epm1*zpd +hpd))*zpd*hpdinv 
               local wxk= Ddx[i,j]*temp 
               local wyk= Ddy[i,j]*temp 
               gqk[i,j,k,1]= qpr*Jac[i,j,k]*(ux[i,j]*wxk +uy[i,j]*wyk) 
               gqk[i,j,k,2]= qpr*Jac[i,j,k]*(vx[i,j]*wxk +vy[i,j]*wyk) 
               gqk[i,j,k,3]= Jac[i,j,k]*(qpr*(wxk*wxk +wyk*wyk) + be2*wz[i,j,k]*wz[i,j,k])                            
            end
         end
      end
                                                                
      for k in 1:NK 
      for i in 0:NI 
         for j in 1:NJ 
            Jifc[i,j,k]= 0.5e0*(Jac[i,j,k]+ Jac[i+1,j,k]) 
            gi[i,j,k,1]= 0.5e0*(g11[i,j] +g11[i+1,j])*Jifc[i,j,k] 
            gi[i,j,k,2]= 0.5e0*(g12[i,j] +g12[i+1,j])*Jifc[i,j,k] 
            gqi[i,j,k,1]= qpr*gi[i,j,k,1]
            gqi[i,j,k,2]= qpr*gi[i,j,k,2]
         end 
      end 
   end 
      for k in 1:NK 
      for i in 1:NI 
         for j in 0:NJ 
            Jjfc[i,j,k]= 0.5e0*(Jac[i,j,k]+ Jac[i,j+1,k]) 
            gj[i,j,k,1]= 0.5e0*(g12[i,j] +g12[i,j+1])*Jjfc[i,j,k] 
            gj[i,j,k,2]= 0.5e0*(g22[i,j] +g22[i,j+1])*Jjfc[i,j,k] 
            gqj[i,j,k,1]= qpr*gj[i,j,k,1]
            gqj[i,j,k,2]= qpr*gj[i,j,k,2]
                                                                        
         end
      end
   end
                                                                      
      for j in 1:NJ 
         for i in 0:NI 
            for k in 1:NK 
               gi3[i,j,k]= 0.5e0*(g13[i,j,k] +g13[i+1,j,k])*Jifc[i,j,k] 
               gqi3[i,j,k]= qpr*gi3[i,j,k]
            end
         end
      end

      for j in 0:NJ 
         for i in 1:NI 
            for k in 1:NK 
               gj3[i,j,k]= 0.5e0*(g23[i,j,k] +g23[i,j+1,k])*Jjfc[i,j,k] 
               gqj3[i,j,k]= qpr*gj3[i,j,k]

            end
         end
      end

#                                                                       
#     vent recomputed using vfbcs                                       
#     v= [(-vx/d2)*uf +(ux/d2)*vf]/Jac , but uf=0                       
#      local j=0                                                              
#      for i in ient1:ient2                                              
#         local d2= ux[i,j]*vy[i,j] -vx[i,j]*uy[i,j]                          
#        for k in 1:NK                                                   
#            vent[i,k]= (ux[i,j]/d2)*vfbcs[i,k]/Jac[i,j]                
#        end
#     end

@static if (cppdefs.fixed_bottom_thickness)
# ==================================================================================
#     Bottom boundary layers - evenly spaced, no stretching
#     -----------------------------------------------------
      for j in 0:NJ+1
         for i in 0:NI+1

            for k in 0:1        #     at cell centers
               local z= zc[i,j,k]
               wz[i,j,k]= 1e0/dzbot
               Jac[i,j,k]= J2d[i,j]/wz[i,j,k] #  z is  dimensionaless and is equal to  z*/DL
               local sig= rc_kind(k)-0.5e0
               wx[i,j,k]= -wz[i,j,k]*Ddx[i,j]
               wy[i,j,k]= -wz[i,j,k]*Ddy[i,j]
               g13[i,j,k]= ux[i,j]*wx[i,j,k] +uy[i,j]*wy[i,j,k]
               g23[i,j,k]= vx[i,j]*wx[i,j,k] +vy[i,j]*wy[i,j,k]
               # g33[i,j,k]= wx[i,j,k]*wx[i,j,k] +wy[i,j,k]*wy[i,j,k] +
            end

            for k in 0:1-1      # at cell faces
               local sig= rc_kind(k)
               wt[i,j,k]= 0e0
               local z= zc[i,j,k]
               local zpd= z +dztop
               wzk[i,j,k]= 1e0/dzbot
            end

            wzk[i,j,1]= 0.5*(wz[i,j,1] + wz[i,j,1+1])

            for k in 0:1
               local z= zf[i,j,k]
               local sig= rc_kind(k)
               local temp= 1e0/dzbot
               local wxk= -Ddx[i,j]*temp
               local wyk= -Ddy[i,j]*temp
               gqk[i,j,k,1]= qpr*Jac[i,j,k]*(ux[i,j]*wxk +uy[i,j]*wyk)
               gqk[i,j,k,2]= qpr*Jac[i,j,k]*(vx[i,j]*wxk +vy[i,j]*wyk)
               gqk[i,j,k,3]= Jac[i,j,k]*(qpr*(wxk*wxk +wyk*wyk) + be2*wzk[i,j,k]*wzk[i,j,k])
            end

         end
      end

      # interpolation of Jac, g11, g12 and gi on the face x. 
      for k in 1:1
         for i in 0:NI
            for j in 1:NJ
               Jifc[i,j,k]= 0.5e0*(Jac[i,j,k]+ Jac[i+1,j,k])
               gi[i,j,k,1]= 0.5e0*(g11[i,j] +g11[i+1,j])*Jifc[i,j,k]
               gi[i,j,k,2]= 0.5e0*(g12[i,j] +g12[i+1,j])*Jifc[i,j,k]
               gqi[i,j,k,1]= qpr*gi[i,j,k,1]
               gqi[i,j,k,2]= qpr*gi[i,j,k,2]
            end
         end
      end

      # interpolation of Jac, g11, g12 and gi on the face y. 
      for k in 1:1
         for i in 1:NI
            for j in 0:NJ
               Jjfc[i,j,k]= 0.5e0*(Jac[i,j,k]+ Jac[i,j+1,k])
               gj[i,j,k,1]= 0.5e0*(g12[i,j] +g12[i,j+1])*Jjfc[i,j,k]
               gj[i,j,k,2]= 0.5e0*(g22[i,j] +g22[i,j+1])*Jjfc[i,j,k]
               gqj[i,j,k,1]= qpr*gj[i,j,k,1]
               gqj[i,j,k,2]= qpr*gj[i,j,k,2]
            end
         end
      end


   end



 for i in 0:NI
   for j in 0:NJ
     for k in 0:NK
       JacInv[i,j,k]=1e0/Jac[i,j,k]
     end
   end
 end

end                                          
