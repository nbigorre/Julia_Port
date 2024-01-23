function uv_geostrophic(ugeo, vgeo)
   #     ---------------------------------------------                     
   # !   Finds the geostrophic velocities ugeo,vgeo
   # USE header
   #  USE rpgrads
   #     modified for periodicew bc                                        
   #     ---------------------------                                       
   #     Sets up the initial velocities so that they are in geostrophic bal
   #     with the initial density field.                                   
   #      implicit logical (a-z)                                           
   #     implicit none 
   #     integer i,j,k,n,imax,jmax,kmax,m 
   #     real(kind=rc_kind):: uxi,vyj,hxi,heta,hx,hy,px,py,ujfc 
   #     real(kind=rc_kind):: res,resmax,kaph1
   #     real(kind=rc_kind):: ainv,be2,fac2,wzsk,wxsk,wysk,Jack,pxi,peta,      &
   #    &     psig,pgrad,con,pz                                            
   #     real(kind=rc_kind), dimension (0:NI+1,0:NJ+1,0:NK+1) :: ugeo, vgeo

   rpevalgrad_Song(0)
   #@ccall "./PSOM_LIB.so".rpevalgrad_song_(Ref(0)::Ref{Int})::Cvoid

   local kaph1 = 1e0 - kappah
   local con = 1e0 - qpr

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

            ugeo[i, j, k] = -(qpr * py + drpy[i, j, k] + gpr * hy) / (ffc[i, j])
            vgeo[i, j, k] = (qpr * px + drpx[i, j, k] + gpr * hx) / (ffc[i, j])
         end
      end
      for k in 1:NK
         ugeo[0, j, k] = ugeo[NI, j, k]
         vgeo[0, j, k] = vgeo[NI, j, k]
         ugeo[NI+1, j, k] = ugeo[1, j, k]
         vgeo[NI+1, j, k] = vgeo[1, j, k]
      end
   end

end



