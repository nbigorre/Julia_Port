

function intpol()
   if (NI < 4 || NJ < 4)
      println("prob in intpol,insuff i pts")
      exit(1)
   end
   local xa::rc_kind = 9e0 / 16e0
   local xb::rc_kind = -1e0 / 16e0

   for i in 0:NI
      local im1::Int = 0
      local ip1::Int = 0 
      local ip2::Int = 0
      local i0::Int = 0
      if (i == 0)
         i0 = NI 
         im1 = NI-1 
         ip1 = i+1 
         ip2 = i+2 
      elseif (i == 1)
         i0 = i 
         im1 = NI 
         ip1 = i+1 
         ip2 = i+2 
      elseif (i == NI)
         i0 = i 
         im1 = i-1 
         ip1 = 1 
         ip2 = 2 
      elseif  (i == NI-1)
         i0= i 
         im1 = i-1 
         ip1 = i+1 
         ip2 = 1 
      else
         i0 = i 
         im1 = i-1 
         ip1 = i+1 
         ip2 = i+2 
      end

      for k in 1:NK
         for j in 1:NJ
            cxf[i+1,j,k] = ( xa * (ux[i0+1, j+1] * cx[i0+1, j+1, k+1] + ux[ip1+1, j+1] * cx[ip1+1, j+1, k+1])
                           + xb * (ux[im1+1, j+1] * cx[im1+1, j+1, k+1] + ux[ip2+1, j+1] * cx[ip2+1, j+1, k+1])
                           + 0.5e0 * (uy[i0+1, j+1] * cy[i0+1, j+1, k+1] + uy[ip1+1, j+1] * cy[ip1+1, j+1, k+1])) * Jifc[i+1, j, k]
         end
      end
   end

   if (!periodicew)
      for k in 1:NK
         for j in 1:NJ
            cxf[1, j, k] = ufbcw[j, k]
            cxf[NI+1, j, k] = ufbce[j, k]
         end
      end
   end

   for j in 1:NJ-1
      local j0::Int = 0
      local jm1::Int = 0
      local jp1::Int = 0
      local jp2::Int = 0
      if (j == 1)
         j0= j 
         jm1 = NJ 
         jp1 = j+1 
         jp2 = j+2 
         xa = 0.5 
         xb = 0.0 
      elseif (j == NJ-1)
         j0 = j 
         jm1 = j-1 
         jp1 = NJ 
         jp2 = 1 
         xa = 0.5 
         xb = 0.0 
      else
         j0 = j 
         jm1 = j-1 
         jp1 = j+1 
         jp2 = j+2 
         xa = 9e0/16e0 
         xb = -1e0/16e0 
      end
      for k in 1:NK
         for i in 1:NI
            cyf[i, j+1, k] = Jjfc[i, j+1, k] * (
               0.5e0 * (vx[i+1, j0+1] * cx[i+1, j0+1, k+1] + vx[i+1, jp1+1] * cx[i+1, jp1+1, k+1])
                + xa * (vy[i+1, j0+1] * cy[i+1, j0+1, k+1] + vy[i+1, jp1+1] * cy[i+1, jp1+1, k+1])
                + xb * (vy[i+1, jm1+1] * cy[i+1, jm1+1, k+1] + vy[i+1, jp2+1] * cy[i+1, jp2+1, k+1])
            )
         end
      end
   end

   for k in 1:NK
      for i in 1:NI
         cyf[i, NJ+1, k] = vfbcn[i,k]
         cyf[i, 1, k] = vfbcs[i,k]
      end
   end

   for k in 1:NK-1
      for j in 1:NJ
         for i in 1:NI
            local Jack::rc_kind = 0.5e0 * (Jac[i+1, j+1, k+1] + Jac[i+1, j+1, k+2]) 
            czf[i, j, k+1] = 0.5 * Jack * ( wx[i+1, j+1, k+1] * cx[i+1, j+1, k+1] + wx[i+1, j+1, k+2] * cx[i+1, j+1, k+2]
                                          + wy[i+1, j+1, k+1] * cy[i+1, j+1, k+1] + wy[i+1, j+1, k+2] * cy[i+1, j+1, k+2]
                                    + EPS *(wz[i+1, j+1, k+1] * cz[i+1, j+1, k+1] + wz[i+1, j+1, k+2] * cz[i+1, j+1, k+2]))
         end
      end
   end

   for j in 1:NJ
      for i in 1:NI
         czf[i, j, 1] = wfbcb[i,j]
         czf[i, j, NK+1] = Jac[i+1, j+1, NK+1] * (wx[i+1, j+1, NK+1] * cx[i+1, j+1, NK+1] + wy[i+1, j+1, NK+1] * cy[i+1, j+1, NK+1] + EPS * wz[i+1, j+1, NK+1] * cz[i+1, j+1, NK+1])
      end
   end

end