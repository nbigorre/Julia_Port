function hsave()
   #----------------------------------------------------                   
   # USE header
   #     saves the hydrostatic pressure gradients at the initial time step 
   #     integer i,j,k 
   #     REAL(kind=rc_kind) :: kap1,hxi,heta 
   #                                                                       
   local kap1 = 1e0 - kappah
   #     hxn,hyn,hzn need to be evaluated at the boundaries too, because   
   #     they are used in mghfill.                                         
   #                                                                       
   #     x-direction                                                       
   #     -----------                                                       
   for j in 1:NJ
      for i in 0:NI
         local hxi = h[i+1, j] - h[i, j]
         local heta = 0.25 * (h[i+1, j+1] + h[i, j+1] - h[i+1, j-1] - h[i, j-1])
         for k in 1:NK
            hxn[i, j, k] = hxi * gi[i, j, k, 1] + heta * gi[i, j, k, 2]
         end
      end
   end
   #                                                                       
   #     y-direction                                                       
   #     -----------                                                       
   for j in 0:NJ
      for i in 1:NI
         local heta = h[i, j+1] - h[i, j]
         local hxi = 0.25 * (h[i+1, j+1] + h[i+1, j] - h[i-1, j+1] - h[i-1, j])
         for k in 1:NK
            hyn[i, j, k] = heta * gj[i, j, k, 2] + hxi * gj[i, j, k, 1]
         end

      end
   end


   #     Now save grad h as is done in uvchy.f for later time steps.       
   for j in 1:NJ
      for i in 1:NI
         local hxi = 0.5e0 * (h[i+1, j] - h[i-1, j])
         local heta = 0.5e0 * (h[i, j+1] - h[i, j-1])
         gradhn[i, j, 1] = ux[i, j] * hxi + vx[i, j] * heta
         gradhn[i, j, 2] = uy[i, j] * hxi + vy[i, j] * heta
      end
   end

end
