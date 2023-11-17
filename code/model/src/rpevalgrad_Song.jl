function rpevalgrad_Song(n)
   local Jy = zeros(rc_kind, (NJ + 1, NK + 1))
   local dzeta = zeros(rc_kind, (NJ, NK + 1))
   local Jy2 = zeros(rc_kind, (NJ, NK + 1))
   local dzsig = zeros(rc_kind, (NJ, NK + 1))
   local Jx = zeros(rc_kind, (NI + 1, NK + 1))
   local dzxi = zeros(rc_kind, (NI + 1, NK + 1))
   local Jx2 = zeros(rc_kind, (NI, NK + 1))
   local dzsig_x = zeros(rc_kind, (NI + 1, NK + 1))
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
         local hpd = h[i+1, j+1] * @fortGet("hdl", rc_kind) + D[i+1, j+1]
         local dep = D[i+1, j+1]
         local ht = h[i+1, j+1] * @fortGet("hdl", rc_kind)
         #     for k=NK (firt for the toppmost layer of cells)                   
         #     z is the dimensional z and is therefore multiplied by DL 
         local k = NK
         local z = DL * zc[i+1, j+1, k+1]
         local zt = DL * zf[i+1, j+1, k+2]
         local zb = DL * zf[i+1, j+1, k+1]
         local rg = (rho[i+1, j+1, k+1] - R0) * vconst
         rp[i+1, j+1, k+1] = rg * (zt - z)
         local temp = rg * (z - zb)
         #     for the rest of the column
         #     remember that findz will not work for the k=0 cell where sig is -v
         #     We will simply fill rp at k=0 using extrapolation - since this    
         #     is what we use to fill s and T in any case.  
         for k in NK-1:-1:1
            z = DL * zc[i+1, j+1, k+1]
            zt = zb
            zb = DL * zf[i+1, j+1, k+1]
            rg = (rho[i+1, j+1, k+1] - R0) * vconst
            rp[i+1, j+1, k+1] = rg * (zt - z) + temp + rp[i+1, j+1, k+2]
            temp = rg * (z - zb)
         end
      end
   end

   
   
   #     Boundary conditions
   for k in 1:NK
      for i in 1:NI
         rp[i+1, 1, k+1] = rp[i+1, 2, k+1]
         rp[i+1, NJ+2, k+1] = rp[i+1, NJ+1, k+1]
      end
      for j in 1:NJ
         rp[NI+2, j+1, k+1] = rp[2, j+1, k+1]
         rp[1, j+1, k+1] = rp[NI+1, j+1, k+1]
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
         
         Jy[j+1, k+1] = 0.5e0 * (rho[i+1, j+2, k+2] + rho[i+1, j+2, k+3] - (rho[i+1, j+1, k+2] + rho[i+1, j+1, k+3]))
         dzeta[j, k+1] = 0.5e0 * (zc[i+1, j+2, k+2] + zc[i+1, j+2, k+3] - (zc[i+1, j+1, k+2] + zc[i+1, j+1, k+3])) * DL
         for k in 1:NK-1
            Jy[j+1, k+1] = 0.5e0 * (rho[i+1, j+2, k+1] + rho[i+1, j+2, k+2] - (rho[i+1, j+1, k+1] + rho[i+1, j+1, k+2]))
            dzeta[j, k+1] = 0.5e0 * (zc[i+1, j+2, k+1] + zc[i+1, j+2, k+2] - (zc[i+1, j+1, k+1] + zc[i+1, j+1, k+2])) * DL
         end
         local k = NK
         Jy[j+1, k+1] = 0.5e0 * (rho[i+1,j+2, k+1] + rho[i+1, j+2, k] - (rho[i+1, j+1, k+1] + rho[i+1, j+1, k]))
         dzeta[j, k+1] = 0.5e0 * (zc[i+1,j+2, k+1] + zc[i+1, j+2, k] - (zc[i+1, j+1, k+1] + zc[i+1, j+1, k])) * DL
      end

      
      for j in 1:NJ-1
         local k = 0
         Jy2[j, 1] = 0.5e0 * (rho[i+1, j+1, k+3] + rho[i+1, j+2, k+3] - (rho[i+1, j+1, k+2] + rho[i+1, j+2, k+2]))
         dzsig[j, 2] = 0.5e0 * (zc[i+1, j+1, 3] + zc[i+1, j+2, 3] - (zc[i+1, j+1, 2] + zc[i+1, j+2, 2])) * DL

         for k in 1:NK-1
            Jy2[j, k+1] = 0.5e0 * (rho[i+1, j+1, k+2] + rho[i+1, j+2, k+2] - (rho[i+1, j+1, k+1] + rho[i+1, j+2, k+1]))
            dzsig[j,k+1] = 0.5e0 * (zc[i+1, j+1, k+2] + zc[i+1, j+2, k+2] - (zc[i+1, j+1, k+1] + zc[i+1, j+2, k+1])) * DL
         end
         Jy2[j, NK+1] = 0.5e0 * (rho[i+1, j+1, NK+1] + rho[i+1, j+2, NK+1] - (rho[i+1, j+1, NK] + rho[i+1, j+2, NK]))
         dzsig[j, NK+1] = 0.5e0 * (zc[i+1, j+1, NK+1] + zc[i+1, j+2, NK+1] - (zc[i+1, j+1, NK] + zc[i+1, j+2, NK])) * DL
      end
      
      
      
      for k in 0:NK
         for j in 1:NJ-1
            Jy[j+1, k+1] *= dzsig[j, k+1]
            Jy2[j, k+1] *= dzeta[j, k+1]
            Jy[j+1, k+1] -= Jy2[j, k+1]
         end
      end
      
      for k in 0:NK
         Jy[1, k+1] = 0e0
         Jy[NJ+1, k+1] = 0e0
      end
      
      for j in 1:NJ-1
         local k = NK
         local Jsum = Jy[j+1, k+1] * 0.5e0
         grpjfc[i, j+1, NK] = Jsum * vconst * gj[i, j+1, k, 2]
         for k in NK-1:-1:1
            Jsum += Jy[j+1, k+1]
            grpjfc[i, j+1, k] = Jsum * vconst * gj[i, j+1, k, 2]
         end
      end
      
      for k in 1:NK
         grpjfc[i,1,k] = 0e0
         grpjfc[i, NJ+1, k] = 0e0
      end
      for j in 1:NJ
         local k = NK
         drpy[i,j,k] = vconst * vy[i+1, j+1] * 0.5e0 * (Jy[j+1, NK+1] + Jy[j, NK+1]) * 0.5e0
         for k in NK-1:-1:1
            drpy[i,j,k] = drpy[i,j,k+1] + vconst * vy[i+1, j+1] * 0.5e0 * (Jy[j+1, k+1] + Jy[j, k+1])
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
         Jx[i+1, k+1] = rho[i+2, j+1, k+2] - rho[i+1, j+1, k+2]
         dzxi[i+1, k+1] = (zc[i+2, j+1, k+2] - zc[i+1, j+1, k+2]) * DL
         for k in 1:NK-1
            Jx[i+1, k+1] = 0.5e0 * (rho[i+2, j+1, k+1] + rho[i+2, j+1, k+2] - (rho[i+1, j+1, k+1] + rho[i+1, j+1, k+2]))
            dzxi[i+1, k+1] = 0.5e0 * (zc[i+2, j+1, k+1] + zc[i+2, j+1, k+2] - (zc[i+1, j+1, k+1] + zc[i+1, j+1, k+2])) * DL
         end
         local k = NK
         Jx[i+1, k+1] = rho[i+2, j+1, k+1] - rho[i+1, j+1, k+1]
         dzxi[i+1, k+1] = (zc[i+2, j+1, k+1] - zc[i+1, j+1, k+1]) * DL
      end
      local i = NI
      local k = 0
      Jx[i+1, k+1] = rho[2, j+1, k+2] - rho[i+1, j+1, k+2]
      dzxi[i+1, k+1] = (zc[2, j+1, k+2] - zc[i+1, j+1, k+2]) * DL
      for k in 1:NK-1
         Jx[i+1, k+1] = 0.5e0 * (rho[2, j+1, k+1] + rho[2, j+1, k+2] - (rho[i+1, j+1, k+1] + rho[i+1, j+1, k+2]))
         dzxi[i+1, k+1] = 0.5e0 * (zc[2, j+1, k+1] + zc[2, j+1, k+2] - (zc[i+1, j+1, k+1] + zc[i+1, j+1, k+2])) * DL
      end
      local k = NK
      Jx[i+1, k+1] = rho[2, j+1, k+1] - rho[i+1, j+1, k+1]
      dzxi[i+1, k+1] = (zc[2, j+1, k+1] - zc[i+1, j+1, k+1]) * DL
      
      
      for i in 1:NI-1
         local k = 0
         Jx2[j, 1] = 0.5e0 * (rho[i+1, j+1, k+3] + rho[i+2, j+1, k+3] - (rho[i+1, j+1, k+2] + rho[i+2, j+2, k+2]))            ##############################
         dzsig_x[i+1, 2] = 0.5e0 * (zc[i+1, j+1, k+3] + zc[i+2, j+1, k+3] - (zc[i+1, j+1, k+2] + zc[i+2, j+2, k+2])) * DL
         for k in 1:NK-1
            Jx2[i, k+1] = 0.5e0 * (rho[i+1, j+1, k+2] + rho[i+2, j+1, k+2] - (rho[i+1, j+1, k+1] + rho[i+2, j+1, k+1]))
            dzsig_x[i+1, k+1] = 0.5e0 * (zc[i+1, j+1, k+2] + zc[i+2, j+1, k+2] - (zc[i+1, j+1, k+1] + zc[i+2, j+1, k+1])) * DL
         end
         Jx2[i, NK+1] = 0.5e0 * (rho[i+1, j+1, NK+1] + rho[i+2, j+1, NK+1] - (rho[i+1, j+1, NK] + rho[i+2, j+1, NK]))
         dzsig_x[i+1, NK+1] = 0.5e0 * (zc[i+1, j+1, NK+1] + zc[i+2, j+1, NK+1] - (zc[i+1, j+1, NK] + zc[i+2, j+1, NK])) * DL
      end
      local i = NI
      local k = 0
      Jx2[j, 1] = 0.5e0 * (rho[i+1, j+1, k+3] + rho[2, j+1, k+3] - (rho[i+1, j+1, k+2] + rho[2, j+1, k+2]))                    ###############################
      dzsig_x[i+1, 2] = 0.5e0 * (zc[i+1, j+1, k+3] + zc[2, j+1, k+3] - (zc[i+1, j+1, k+2] + zc[2, j+1, k+2])) * DL
      for k in 1:NK-1
         Jx2[i, k+1] = 0.5e0 * (rho[i+1, j+1, k+2] + rho[2, j+1, k+2] - (rho[i+1, j+1, k+1] + rho[2, j+1, k+1]))
         dzsig_x[i+1, k+1] = 0.5e0 * (zc[i+1, j+1, k+2] + zc[2, j+1, k+2] - (zc[i+1, j+1, k+1] + zc[2, j+1, k+1])) * DL
      end
      Jx2[i, NK+1] = 0.5e0 * (rho[i+1, j+1, NK+1] + rho[2, j+1, NK+1] - (rho[i+1, j+1, NK] + rho[2, j+1, NK]))
      dzsig_x[i+1, NK+1] = 0.5e0 * (zc[i+1, j+1, NK+1] + zc[2, j+1, NK+1] - (zc[i+1, j+1, NK] + zc[2, j+1, NK])) * DL
      
      for k in 0:NK
         for i in 1:NI
            Jx[i+1, k+1] *= dzsig_x[i+1, k+1]
            Jx2[i, k+1] *= dzxi[i+1, k+1]
            Jx[i+1, k+1] -= Jx2[i, k+1]
         end
      end

      
      for k in 0:NK
         Jx[1, k+1] = Jx[NI+1, k+1]
      end

      for i in 1:NI
         local k = NK
         local Jsum = Jx[i+1, k+1] * 0.5e0
         grpifc[i+1, j, NK] = Jsum * vconst * gi[i+1, j, k, 1]
         for k in NK-1:-1:1
            Jsum += Jx[i+1, k+1]
            grpifc[i+1, j, k] = Jsum * vconst * gi[i+1, j, k, 1]
         end
      end
      for k in 1:NK
         grpifc[1, j, k] = grpifc[NI+1, j, k]
      end
      for i in 1:NI
         local k = NK
         drpx[i,j,k] = vconst * ux[i+1, j+1] * 0.5e0 * (Jx[i+1, NK+1] + Jx[i, NK+1]) * 0.5e0
         for k in NK-1:-1:1
            drpx[i,j,k] = drpx[i,j,k+1] + vconst * ux[i+1, j+1] * 0.5e0 * (Jx[i+1, k+1] + Jx[i, k+1])
         end
      end

   end

end