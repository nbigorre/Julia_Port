function rpevalgrad_Song(n)
   local Jy = OffsetArray(zeros(rc_kind, (NJ + 1, NK + 1)), 0:NJ, 0:NK)
   local dzeta = OffsetArray(zeros(rc_kind, (NJ, NK + 1)), 1:NJ, 0:NK)
   local Jy2 = OffsetArray(zeros(rc_kind, (NJ, NK + 1)), 1:NJ, 0:NK)
   local dzsig = OffsetArray(zeros(rc_kind, (NJ, NK + 1)), 1:NJ, 0:NK)
   local Jx = OffsetArray(zeros(rc_kind, (NI + 1, NK + 1)), 0:NI, 0:NK)
   local dzxi = OffsetArray(zeros(rc_kind, (NI + 1, NK + 1)), 0:NI, 0:NK)
   local Jx2 = OffsetArray(zeros(rc_kind, (NI, NK + 1)), 1:NI, 0:NK)
   local dzsig_x = OffsetArray(zeros(rc_kind, (NI + 1, NK + 1)), 0:NI, 0:NK)
   if (!rect)
      println("modify grpifc,grpjfc")
      exit(1)
   end
   #      Evaluate Jacobian at the edge-centers
   #      multiply press and press gradient by const
   local vconst = 10e0 * gpr / @fortGet("p1", rc_kind)

   #@ccall "./PSOM_LIB.so".findzall_()::Cvoid
   findzall()
   
   #@ccall "./PSOM_LIB.so".evalrho_(pointer(rho)::Ptr{rc_kind}, Ref(n)::Ref{Int})::Cvoid
   evalrho(rho, n)
   # ===================================================================================



   for i in 0:NI+1
      for j in 0:NJ+1
         local hpd = h[i+1, j+1] * @fortGet("hdl", rc_kind) + D[i, j]
         local dep = D[i, j]
         local ht = h[i+1, j+1] * @fortGet("hdl", rc_kind)
         #     for k=NK (firt for the toppmost layer of cells)                   
         #     z is the dimensional z and is therefore multiplied by DL 
         local k = NK
         local z = DL * zc[i, j, k]
         local zt = DL * zf[i, j, k]
         local zb = DL * zf[i, j, k-1]
         local rg = (rho[i, j, k] - R0) * vconst
         rp[i, j, k] = rg * (zt - z)
         local temp = rg * (z - zb)
         #     for the rest of the column
         #     remember that findz will not work for the k=0 cell where sig is -v
         #     We will simply fill rp at k=0 using extrapolation - since this    
         #     is what we use to fill s and T in any case.  
         for k in NK-1:-1:1
            z = DL * zc[i, j, k]
            zt = zb
            zb = DL * zf[i, j, k]
            rg = (rho[i, j, k] - R0) * vconst
            rp[i, j, k] = rg * (zt - z) + temp + rp[i, j, k+1]
            temp = rg * (z - zb)
         end
      end
   end

   
   
   #     Boundary conditions
   for k in 1:NK
      for i in 1:NI
         rp[i, 0, k] = rp[i, 1, k]
         rp[i, NJ+1, k] = rp[i, NJ, k]
      end
      for j in 1:NJ
         rp[NI+1, j, k] = rp[1, j, k]
         rp[0, j, k] = rp[NI, j, k]
      end
   end
   
   
   # ===================================================================================
   # c     y-direction
   # c     -----------
   # c     sloping sigma surfaces
   # c     solid boundaries
   
   for i in 1:NI
      
      for j in 1:NJ-1
         local k = 0
         
         Jy[j, k] = 0.5e0 * (rho[i, j+1, k+1] + rho[i, j+1, k+2] - (rho[i, j, k+1] + rho[i, j, k+2]))
         dzeta[j, k] = 0.5e0 * (zc[i, j+1, k+1] + zc[i, j+1, k+2] - (zc[i, j, k+1] + zc[i, j, k+2])) * DL
         for k in 1:NK-1
            Jy[j, k] = 0.5e0 * (rho[i, j+1, k] + rho[i, j+1, k+1] - (rho[i, j, k] + rho[i, j, k+1]))
            dzeta[j, k] = 0.5e0 * (zc[i, j+1, k] + zc[i, j+1, k+1] - (zc[i, j, k] + zc[i, j, k+1])) * DL
         end
         local k = NK
         Jy[j, k] = 0.5e0 * (rho[i,j+1, k] + rho[i, j+1, k-1] - (rho[i, j, k] + rho[i, j, k-1]))
         dzeta[j, k] = 0.5e0 * (zc[i,j+1, k] + zc[i, j+1, k-1] - (zc[i, j, k] + zc[i, j, k-1])) * DL
      end

      
      for j in 1:NJ-1
         local k = 0
         Jy2[j, 0] = 0.5e0 * (rho[i, j, k+2] + rho[i, j+1, k+2] - (rho[i, j, k+1] + rho[i, j+1, k+1]))
         dzsig[j, 0] = 0.5e0 * (zc[i, j, 2] + zc[i, j+1, 2] - (zc[i, j, 1] + zc[i, j+1, 1])) * DL

         for k in 1:NK-1
            Jy2[j, k+1] = 0.5e0 * (rho[i, j, k+1] + rho[i, j+1, k+1] - (rho[i, j, k] + rho[i, j+1, k]))
            dzsig[j,k] = 0.5e0 * (zc[i, j, k+1] + zc[i, j+1, k+1] - (zc[i, j, k] + zc[i, j+1, k])) * DL
         end
         Jy2[j, NK] = 0.5e0 * (rho[i, j, NK] + rho[i, j+1, NK] - (rho[i, j, NK-1] + rho[i, j+1, NK-1]))
         dzsig[j, NK] = 0.5e0 * (zc[i, j, NK] + zc[i, j+1, NK] - (zc[i, j, NK-1] + zc[i, j+1, NK-1])) * DL
      end
      
      
      
      for k in 0:NK
         for j in 1:NJ-1
            Jy[j, k] *= dzsig[j, k]
            Jy2[j, k] *= dzeta[j, k]
            Jy[j, k] -= Jy2[j, k]
         end
      end
      
      for k in 0:NK
         Jy[0, k] = 0e0
         Jy[NJ, k] = 0e0
      end
      
      for j in 1:NJ-1
         local k = NK
         local Jsum = Jy[j, k] * 0.5e0
         grpjfc[i, j+1, NK] = Jsum * vconst * gj[i, j, k, 2]
         for k in NK-1:-1:1
            Jsum += Jy[j, k]
            grpjfc[i, j+1, k] = Jsum * vconst * gj[i, j, k, 2]
         end
      end
      
      for k in 1:NK
         grpjfc[i,1,k] = 0e0
         grpjfc[i, NJ+1, k] = 0e0
      end
      for j in 1:NJ
         local k = NK
         drpy[i,j,k] = vconst * vy[i, j] * 0.5e0 * (Jy[j, NK] + Jy[j-1, NK]) * 0.5e0
         for k in NK-1:-1:1
            drpy[i,j,k] = drpy[i,j,k+1] + vconst * vy[i, j] * 0.5e0 * (Jy[j, k] + Jy[j-1, k])
         end
      end
  end

   
   #   ===================================================================================
   # c     x- direction
   # c     ------------------------
   # c     sloping sigma surfaces
   # c     periodic boundaries
   
   for j in 1:NJ

      for i in 1:NI-1
         local k = 0
         Jx[i, k] = rho[i+1, j, k+1] - rho[i, j, k+1]
         dzxi[i, k] = (zc[i+1, j, k+1] - zc[i, j, k+1]) * DL
         for k in 1:NK-1
            Jx[i, k] = 0.5e0 * (rho[i+1, j, k] + rho[i+1, j, k+1] - (rho[i, j, k] + rho[i, j, k+1]))
            dzxi[i, k] = 0.5e0 * (zc[i+1, j, k] + zc[i+1, j, k+1] - (zc[i, j, k] + zc[i, j, k+1])) * DL
         end
         local k = NK
         Jx[i, k] = rho[i+1, j, k] - rho[i, j, k]
         dzxi[i, k] = (zc[i+1, j, k] - zc[i, j, k]) * DL
      end
      local i = NI
      local k = 0
      Jx[i, k] = rho[1, j, k+1] - rho[i, j, k+1]
      dzxi[i, k] = (zc[1, j, k+1] - zc[i, j, k+1]) * DL
      for k in 1:NK-1
         Jx[i, k] = 0.5e0 * (rho[1, j, k] + rho[1, j, k+1] - (rho[i, j, k] + rho[i, j, k+1]))
         dzxi[i, k] = 0.5e0 * (zc[1, j, k] + zc[1, j, k+1] - (zc[i, j, k] + zc[i, j, k+1])) * DL
      end
      local k = NK
      Jx[i, k] = rho[1, j, k] - rho[i, j, k]
      dzxi[i, k] = (zc[1, j, k] - zc[i, j, k]) * DL
      
      
      for i in 1:NI-1
         local k = 0
         Jx2[j, 0] = 0.5e0 * (rho[i, j, k+2] + rho[i+1, j, k+2] - (rho[i, j, k+1] + rho[i+1, j+1, k+1]))            ##############################
         dzsig_x[i, 1] = 0.5e0 * (zc[i, j, k+2] + zc[i+1, j, k+2] - (zc[i, j, k+1] + zc[i+1, j+1, k+1])) * DL
         for k in 1:NK-1
            Jx2[i, k] = 0.5e0 * (rho[i, j, k+1] + rho[i+1, j, k+1] - (rho[i, j, k] + rho[i+1, j, k]))
            dzsig_x[i, k] = 0.5e0 * (zc[i, j, k+1] + zc[i+1, j, k+1] - (zc[i, j, k] + zc[i+1, j, k])) * DL
         end
         Jx2[i, NK] = 0.5e0 * (rho[i, j, NK] + rho[i+1, j, NK] - (rho[i, j, NK-1] + rho[i+1, j, NK-1]))
         dzsig_x[i, NK] = 0.5e0 * (zc[i, j, NK] + zc[i+1, j, NK] - (zc[i, j, NK-1] + zc[i+1, j, NK-1])) * DL
      end
      local i = NI
      local k = 0
      Jx2[j, 0] = 0.5e0 * (rho[i, j, k+2] + rho[1, j, k+2] - (rho[i, j, k+1] + rho[1, j, k+1]))                    ###############################
      dzsig_x[i, 1] = 0.5e0 * (zc[i, j, k+2] + zc[1, j, k+2] - (zc[i, j, k+1] + zc[1, j, k+1])) * DL
      for k in 1:NK-1
         Jx2[i, k] = 0.5e0 * (rho[i, j, k+1] + rho[1, j, k+1] - (rho[i, j, k] + rho[1, j, k]))
         dzsig_x[i, k] = 0.5e0 * (zc[i, j, k+1] + zc[1, j, k+1] - (zc[i, j, k] + zc[1, j, k])) * DL
      end
      Jx2[i, NK] = 0.5e0 * (rho[i, j, NK] + rho[1, j, NK] - (rho[i, j, NK-1] + rho[1, j, NK-1]))
      dzsig_x[i, NK] = 0.5e0 * (zc[i, j, NK] + zc[1, j, NK] - (zc[i, j, NK-1] + zc[1, j, NK-1])) * DL
      
      for k in 0:NK
         for i in 1:NI
            Jx[i, k] *= dzsig_x[i, k]
            Jx2[i, k] *= dzxi[i, k]
            Jx[i, k] -= Jx2[i, k]
         end
      end

      
      for k in 0:NK
         Jx[0, k] = Jx[NI, k]
      end

      for i in 1:NI
         local k = NK
         local Jsum = Jx[i, k] * 0.5e0
         grpifc[i+1, j, NK] = Jsum * vconst * gi[i, j, k, 1]
         for k in NK-1:-1:1
            Jsum += Jx[i, k]
            grpifc[i+1, j, k] = Jsum * vconst * gi[i, j, k, 1]
         end
      end
      for k in 1:NK
         grpifc[1, j, k] = grpifc[NI+1, j, k]
      end
      for i in 1:NI
         local k = NK
         drpx[i,j,k] = vconst * ux[i, j] * 0.5e0 * (Jx[i, NK] + Jx[i-1, NK]) * 0.5e0
         for k in NK-1:-1:1
            drpx[i,j,k] = drpx[i,j,k+1] + vconst * ux[i, j] * 0.5e0 * (Jx[i, k] + Jx[i-1, k])
         end
      end

   end

end