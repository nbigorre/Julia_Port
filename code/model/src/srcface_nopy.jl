function srcface_nopy(n)
   #!----------------------------------------------------                   
   #  USE header
   #  !USE rpgrads
   #!     FOR PERIODICEW boundaries                                         
   #!     --------------------------                                        
   #!     We interpolate the source terms onto the cell faces               
   #!     It is important that u,v,rp have been filled outside the boundarie
   #      implicit none
   #      integer i,j,k,n,in 
   #      REAL(kind=rc_kind) ::  uint,vint,wint,uxi(0:NI,NJ),uyi(0:NI,NJ),        &
   #     &     vxj(NI,0:NJ),vyj(NI,0:NJ),rpxi,rpeta,rpsig,fac,fac2,         &
   #     &     ainv,be2,vcif,vcjf, Jack                                     
   #      REAL(kind=rc_kind) ::  hxi, heta, gradh, hy, py,px 
   #      REAL(kind=rc_kind) ::  wxsk,wysk,wzsk
   #                                                                        
   #!     We are using the values at the ghost points.                      
   local uxi = OffsetArrays.zeros(rc_kind, 0:NI, 1:NJ)                  
   local uyi = OffsetArrays.zeros(rc_kind, 0:NI, 1:NJ)
   local vxj = OffsetArrays.zeros(rc_kind, 1:NI, 0:NJ)                  
   local vyj = OffsetArrays.zeros(rc_kind, 1:NI, 0:NJ)
   local fac = EPS * delta
   local ainv = 1e0 / apr
   local fac2 = EPS * lambda
   local be2 = beta * EPS * EPS
   #c      for i in 0:NI                                                    
   for j in 1:NJ
      for i in 1:NI-1
         uxi[i, j] = 0.5 * (ux[i+1, j] + ux[i, j])
         uyi[i, j] = 0.5 * (uy[i+1, j] + uy[i, j])
      end
      #     periodic e-w boundaries                                           
      uxi[NI, j] = 0.5 * (ux[1, j] + ux[NI, j])
      uyi[NI, j] = 0.5 * (uy[1, j] + uy[NI, j])
      uxi[0, j] = uxi[NI, j]
      uyi[0, j] = uyi[NI, j]
   end

   for k in 1:NK
      for j in 1:NJ
         for i in 1:NI-1
            local vint = 0.5 * (v[i+1, j, k, n] + v[i, j, k, n])
            local uint = 0.5 * (u[i+1, j, k, n] + u[i, j, k, n])
            local wint = 0.5 * (w[i+1, j, k, n] + w[i, j, k, n])
            local vcif = 0.5 * EPS * (uvis[i+1, j, k] + uvis[i, j, k])
            local vcjf = 0.5 * EPS * (vvis[i+1, j, k] + vvis[i, j, k])
            sifc[i, j, k] = ((uxi[i, j] * (-ffi[i, j] * vint + fac * bbi[i, j] * wint - vcif) + uyi[i, j] * (ffi[i, j] * uint - vcjf)) * Jifc[i, j, k] + grpifc[i, j, k])
         end
      end
   end

   #                                                                       
   #     periodic-ew boundaries                                            
   local i = NI
   for j in 1:NJ
      for k in 1:NK
         local vint = 0.5 * (v[1, j, k, n] + v[NI, j, k, n])
         local uint = 0.5 * (u[1, j, k, n] + u[NI, j, k, n])
         local wint = 0.5 * (w[1, j, k, n] + w[NI, j, k, n])
         local vcif = 0.5 * EPS * (uvis[1, j, k] + uvis[NI, j, k])
         local vcjf = 0.5 * EPS * (vvis[1, j, k] + vvis[NI, j, k])
         sifc[i, j, k] = ((uxi[i, j] * (-ffi[i, j] * vint + fac * bbi[i, j] * wint - vcif) + uyi[i, j] * (ffi[i, j] * uint - vcjf)) * Jifc[i, j, k] + grpifc[i, j, k])
      end
   end

   for k in 1:NK
      for j in 1:NJ
         sifc[0, j, k] = sifc[NI, j, k]
      end
   end


   # sifc is computed.
   # ------------------


   for j in 0:NJ
      for i in 1:NI
         vxj[i, j] = 0.5 * (vx[i, j+1] + vx[i, j])
         vyj[i, j] = 0.5 * (vy[i, j+1] + vy[i, j])
      end
   end

   for k in 1:NK
      for i in 1:NI
         for j in 1:NJ-1
            local vint = 0.5 * (v[i, j+1, k, n] + v[i, j, k, n])
            local uint = 0.5 * (u[i, j+1, k, n] + u[i, j, k, n])
            local wint = 0.5 * (w[i, j+1, k, n] + w[i, j, k, n])
            local vcif = 0.5 * EPS * (uvis[i, j+1, k] + uvis[i, j, k])
            local vcjf = 0.5 * EPS * (vvis[i, j+1, k] + vvis[i, j, k])
            sjfc[i, j, k] = ((vxj[i, j] * (-ffj[i, j] * vint + fac * bbj[i, j] * wint - vcif) + vyj[i, j] * (ffc[i, j] * uint - vcjf)) * Jjfc[i, j, k] + grpjfc[i, j, k])
         end

         for j in 0:NJ:NJ
            local in = 0
            if (j == 0)
               in = 1
            end
            if (j == NJ)
               in = NJ
            end
            local vint = v[i, in, k, n]
            local uint = u[i, in, k, n]
            local wint = w[i, in, k, n]
            local vcif = EPS * uvis[i, in, k]
            local vcjf = EPS * vvis[i, in, k]
            sjfc[i, j, k] = ((vxj[i, j] * (-ffj[i, j] * vint + fac * bbj[i, j] * wint - vcif) + vyj[i, j] * (ffc[i, j] * uint - vcjf)) * Jjfc[i, j, k] + grpjfc[i, j, k])

         end
      end
   end


   # sjfc is computed.
   # ------------------

end