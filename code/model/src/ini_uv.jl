function ini_uv(n)
   #     ---------------------------------------------                     
   # USE header
   #  USE rpgrads
   #     modified for periodicew bc                                        
   #     ---------------------------                                       
   #     Sets up the initial velocities so that they are in geostrophic bal
   #     with the initial density field.                                   
   #      IMPLICIT NONE                                           
   #     integer i,j,k,n,iter,imax,jmax,kmax,m 
   #     REAL(kind=rc_kind) :: uxi,vyj,hxi,heta,hx,hy,px,py,ujfc,pyjfc 
   #     REAL(kind=rc_kind) :: res,resmax,fiter,fzero 
   #     REAL(kind=rc_kind) :: ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
   #    &     psig,pgrad,con,pz                                            
   #     REAL(kind=rc_kind) ::  kaph1
   #     ! pyjfc added on 2012/02/14                                                                      


   local imax = 0
   local jmax = 0
   local kmax = 0

   if (use_shchepetkin == 0) # calculates baroclinic press gradients
      rpevalgrad_Song(n)
      #@ccall "./PSOM_LIB.so".rpevalgrad_song_(Ref(n)::Ref{Int})::Cvoid
   else
      rpevalgrad_Sche(n)
      #@ccall "./PSOM_LIB.so".rpevalgrad_sche_(Ref(n)::Ref{Int})::Cvoid

      drpx[1:NI, 1:NJ, 1:NK] = ru4_Sche[1:NI, 1:NJ, 1:NK]
      drpy[1:NI, 1:NJ, 1:NK] = rv4_Sche[1:NI, 1:NJ, 1:NK]
      grpifc[0:NI, 1:NJ, 1:NK] = ru2_Sche[0:NI, 1:NJ, 1:NK]
      grpjfc[1:NI, 0:NJ, 1:NK] = rv2_Sche[1:NI, 0:NJ, 1:NK]
   end

   local kaph1 = 1e0 - kappah
   local fzero = 0e0
   local fiter = 0e0
   local iter = 1
   local con = 1e0 - qpr
   @label iniuvg101
   local resmax = 1e-20
   for j in 1:NJ
      for i in 1:NI
         local hxi = 0.5e0 * (h[i+1, j] - h[i-1, j])
         local heta = 0.5e0 * (h[i, j+1] - h[i, j-1])
         local hx = ux[i, j] * hxi + vx[i, j] * heta
         local hy = uy[i, j] * hxi + vy[i, j] * heta
         for k in 1:NK
            local pxi = 0.5e0 * (p[i+1, j, k] - p[i-1, j, k])
            local peta = 0.5e0 * (p[i, j+1, k] - p[i, j-1, k])
            local psig = 0.5e0 * (p[i, j, k+1] - p[i, j, k-1])
            local px = (ux[i, j] * pxi + vx[i, j] * peta + wx[i, j, k] * psig)
            local py = (uy[i, j] * pxi + vy[i, j] * peta + wy[i, j, k] * psig)

            #               res = u[i,j,k,n] + (qpr*py +drpy[i,j,k]+gpr*hy)/ffc[i,j]
            #c               res = ffc[i,j]*u[i,j,k,n]*con + (qpr*py +drpy[i,j,k]*con+gpr*hy)                  
            local res = ffc[i, j] * u[i, j, k, n] + (qpr * py + drpy[i, j, k] + gpr * hy)
            if (abs(res) > abs(resmax))
               resmax = res
               imax = i
               jmax = j
               kmax = k
            end

            #               u[i,j,k,n] = - (qpr*py +con*drpy[i,j,k]+gpr*hy)/ (ffc[i,j]*con)                                     
            #               v[i,j,k,n] =   (qpr*px +con*drpx[i,j,k]+gpr*hx)/ (ffc[i,j]*con)                                     
            u[i, j, k, n] = -(qpr * py + drpy[i, j, k] + gpr * hy) / (ffc[i, j])

            v[i, j, k, n] = (qpr * px + drpx[i, j, k] + gpr * hx) / (ffc[i, j])

         end
      end
      for k in 1:NK
         u[0, j, k, n] = u[NI, j, k, n]
         v[0, j, k, n] = v[NI, j, k, n]
         u[NI+1, j, k, n] = u[1, j, k, n]
         v[NI+1, j, k, n] = v[1, j, k, n]
      end
   end

   #     COMPUTE NH PRESSURE                                               
   if (fnhhy > 0.9)
      coriolis(n)
      #@ccall "./PSOM_LIB.so".coriolis_(Ref(n)::Ref{Int})::Cvoid
      
      srcface_nopy(n)
      #@ccall "./PSOM_LIB.so".srcface_nopy_(Ref(n)::Ref{Int})::Cvoid

      calcfkfc()
      #@ccall "./PSOM_LIB.so".calcskfc_()::Cvoid

      local fac2 = EPS * lambda
      local be2 = beta * EPS * EPS
      local ainv = 1e0 / apr
      for j in 1:NJ
         for i in 1:NI
            for k in NK:-1:0
               #     pz= (p[i,j,k+1] -p[i,j,k])*gqk[i,j,k,3]                           
               #     pz = -skfc                                                        
               local pgrad = fiter * (
                  0.25 * (p[i+1, j, k+1]
                          +
                          p[i+1, j, k] - p[i-1, j, k+1] - p[i-1, j, k]) * gqk[i, j, k, 1]
                  +
                  0.25 * (p[i, j+1, k+1] + p[i, j+1, k] - p[i, j-1, k+1]
                          -
                          p[i, j-1, k]) * gqk[i, j, k, 2])

               if (k == NK)
                  #     Use the BC that p(NK)+p(NK+1) = 0.0  (p=0 at the free-surf)       
                  p[i, j, k] = 0.5 * (pgrad + skfc[i, j, k]) / gqk[i, j, k, 3]
               else
                  p[i, j, k] = (p[i, j, k+1] * gqk[i, j, k, 3] + pgrad + skfc[i, j, k]) / gqk[i, j, k, 3]
               end
            end
         end
      end
      mgpfill(fzero, p)
      #@ccall "./PSOM_LIB.so".mgpfill_(Ref(fzero)::Ref{rc_kind}, pointer(p)::Ptr{rc_kind})::Cvoid
      #     filter                                                            
      for m in 1:8
         for k in 0:NK+1
            for j in 1:NJ
               for i in 0:NI+1
                  p[i, j, k] = 0.5 * p[i, j, k] + 0.25 * (p[i, j-1, k] + p[i, j+1, k])
               end
            end
         end
      end
      @goto iniuvg301
      #     -------------- p at boundaries (mgpfill not suitable)-----        
      #     faces j=0,j=NJ                                                    
      #     --------------                                                    
      local j = 0
      for k in 1:NK
         for i in 1:NI
            local vyj = 0.5 * (vy[i, j] + vy[i, j+1])
            #               ujfc=(15.0*u[i,1,k,n]-10.*u[i,2,k,n]+3.0*u[i,3,k,n])/8.0                               
            #               ujfc= u[i,1,k,n]    
            u[i, 0, k, n] = 2e0 * u[i, 1, k, n] - u[i, 2, k, n]
            ujfc = 0.5 * (u[i, 1, k, n] + u[i, 0, k, n])
            pyjfc = -(gpr * hyn[i, j, k] + grpjfc[i, j, k] + ujfc * Jjfc[i, j, k] * vyj * ffc[i, j])
            p[i, j, k] = p[i, j+1, k] - pyjfc / gqj[i, j, k, 2]
            p[i, 0, k] = 2e0 * p[i, 1, k] - p[i, 2, k]
         end
      end
      j = NJ
      for k in 1:NK
         for i in 1:NI
            vyj = 0.5 * (vy[i, j] + vy[i, j+1])
            #ujfc = (15.0 * u[i, NJ, k, n] - 10.0 * u[i, NJ-1, k, n] + 3.0 * u[i, NJ-2, k, n]) / 8.0
            #ujfc = u[i, NJ, k, n]
            u[i, NJ+1, k, n] = 2.0 * u[i, NJ, k, n] - u[i, NJ-1, k, n]
            ujfc = 0.5 * (u[i, NJ, k, n] + u[i, NJ+1, k, n])
            pyjfc = -(gpr * hyn[i, j, k] + grpjfc[i, j, k] + ujfc * Jjfc[i, j, k] * vyj * ffc[i, j])
            p[i, j+1, k] = p[i, j, k] + pyjfc / gqj[i, j, k, 2]
            p[i, NJ+1, k] = 2.0 * p[i, NJ, k] - p[i, NJ-1, k]
         end
      end
      local k = NK
      for j in 0:NJ+1
         for i in 1:NI
            p[i, j, k+1] = -p[i, j, k]
         end
      end
      #     periodic e-w boundaries                                           
      for k in 0:NK+1
         for j in 1:NJ
            p[NI+1, j, k] = p[1, j, k]
            p[0, j, k] = p[NI, j, k]
         end
      end
      @label iniuvg301
      #     --------------------------------                                  
      iter = iter + 1
      fiter = 1.0
   end
   #     end computation of NH pressure                                    
   #!====                                                                   
   #!      if (iter.eq.24) goto 202                                         
   #!==                                                                     
   if ((iter < 50) && (fnhhy >= 0.001) && (abs(resmax) >= 1e-15))
      @goto iniuvg101
   end
   #                                                                       
   #     -------                                                           
   #     filter                                                            
   #      for m in 1:3                                                         
   #      for k in 0:NK+1                                                      
   #         for j in 1:NJ                                                     
   #            for i in 0:NI+1                                                
   #               p[i,j,k]= 0.5*p[i,j,k]+0.25*(p[i,j-1,k]+p[i,j+1,k])     
   #            end                                                     
   #         end                                                        
   #      end                                                           
   #      end                                                                                                                   
   for j in 1:NJ
      for i in 0:NI
         #            vyj= 0.5*(vy[i,j+1]+vy[i,j])                               
         for k in 1:NK
            uf[i, j, k] = (ux[i+1, j] * u[i+1, j, k, n] + ux[i, j] * u[i, j, k, n] + uy[i+1, j] * v[i+1, j, k, n] + uy[i, j] * v[i, j, k, n]) * 0.5 * Jifc[i, j, k]
            #               uf[i,j,k]= grpjfc[i,j,k]/vyj*ffj[i,j]*Jjfc[i,j,k]       
         end
      end
   end
   for j in 0:NJ
      for i in 1:NI
         #            hxi= 0.25*(h[i+1,j+1] +h[i+1,j] -h[i-1,j+1] -h[i-1,j])     
         #            heta= h[i,j+1] -h[i,j]                                     
         #            uxi= 0.5*(ux[1,j]+ux[NI,j])                                
         for k in 1:NK
            #               hy= gj[i,j,k,1]*hxi +gj[i,j,k,2]*heta                   
            #               gradh= gpr*(kappah*hy + kaph1*hyn[i,j,k])               
            #               vf[i,j,k]= -gradh -grpjfc[i,j,k]/vxj*ffi[i,j]*Jifc[i,j,k]
            vf[i, j, k] = 0.5 * Jjfc[i, j, k] * (vx[i, j+1] * u[i, j+1, k, n] + vx[i, j] * u[i, j, k, n] + vy[i, j+1] * v[i, j+1, k, n] + vy[i, j] * v[i, j, k, n])
         end
      end
   end


end
