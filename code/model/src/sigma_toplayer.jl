global const sigma_g13 = OffsetArray(zeros(rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const sigma_g23 = OffsetArray(zeros(rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)

function sigma()
  local g13 = sigma_g13
  local g23 = sigma_g23

  local dnkm1 = rc_kind(NK-1)
  local dnk = rc_kind(NK)

  local be2 = beta * EPS * EPS

  for j in 0:NJ+1
    for i in 0:NI+1
      # All these variables are functions of time
      local hpd = h[i, j] * HDL + dztop
      local hpdinv = 1e0 / hpd

      # Computation of hu:
      local hu = 0e0
      if (i == 0)
        hu = HDL * (h[i+1, j] - h[NI-1, j])
      elseif (i == NI +1)
        hu = HDL * (h[2, j] - h[i-1, j])
      else
        hu = 0.5e0 * HDL * (h[i+1, j] - h[i-1, j])
      end

      # Computation of hv:
      local hv = 0e0
      if (j == 0)
        hv = HDL * (h[i, j+1] - h[i, j])
      elseif (j == NJ +1)
        hv = HDL * (h[i, j] - h[i, j-1])
      else
        hv = 0.5e0 * HDL * (h[i, j+1] - h[i, j-1])
      end

      # Computation of hx and hy:
      local hx = hu * ux[i, j] + hv * vx[i, j]
      local hy = hu * uy[i, j] + hv * vy[i, j]

      #     wz is not a function of depth when the sigma lines are equally spa
      #     Hence wz is wz(i,j,time). For a stretched grid wz would be w(i,j,k
      #     Then Jac which is now Jac(i,j,time) would  become  Jac(i,j,k,time)
      
      # Computation of wx, wy, g13, g23:

      for k in NK:NK+1
        wz[i,j,k] = hpdinv
        Jac[i, j, k] = J2d[i, j] / wz[i, j, k]
        #     All these variables are functions of time and depth               
        #     now hdt computed in vhydro already contains HDL  
        local sig = rc_kind(k) - 0.5e0
        local z = (sig - dnkm1) * hpd - dztop
        wx[i, j, k] = (dnkm1 - sig) * hx * hpdinv
        wy[i, j, k] = (dnkm1 - sig) * hy * hpdinv
        g13[i, j, k] = ux[i, j] * wx[i, j, k] + uy[i, j] * wy[i, j, k]
        g23[i, j, k] = vx[i, j] * wx[i, j, k] + vy[i, j] * wy[i, j, k]
      end

      for k in NK-1:NK
        local sig = rc_kind(k)
        wt[i, j, k] = (dnkm1 - sig) * hdt[i, j] * hpdinv
      end

      wzk[i, j, NK] = hpdinv
      wzk[i, j, NK-1] = 0.5e0 * (wz[i, j, NK] + wz[i, j, NK-1])

      for k in NK-1:NK
        local sig = rc_kind(k)
        local wxk = (dnkm1 - sig) * hx * hpdinv
        local wyk = (dnkm1 - sig) * hy * hpdinv
        gqk[i, j, k, 1] = qpr * Jac[i, j, k] * (ux[i, j] * wxk +  uy[i, j] * wyk)
        gqk[i, j, k, 2] = qpr * Jac[i, j, k] * (vx[i, j] * wxk +  vy[i, j] * wyk)
        gqk[i, j, k, 3] = Jac[i, j, k] * (qpr * (wxk*wxk + wyk*wyk) + be2 * wz[i, j, k] * wz[i, j, k])

      end

    end
  end


  
  for i in 0:NI
    for j in 1:NJ
      Jifc[i, j, NK] = 0.5e0 * (Jac[i, j, NK] + Jac[i+1, j, NK])
      gi[i, j, NK, 1] = 0.5e0 * (g11[i, j] + g11[i+1, j]) * Jifc[i, j, NK]
      gi[i, j, NK, 2] = 0.5e0 * (g12[i, j] + g12[i+1, j]) * Jifc[i, j, NK]
      gqi[i, j, NK, 1] = qpr * gi[i, j, NK, 1]
      gqi[i, j, NK, 2] = qpr * gi[i, j, NK, 2]
    end
  end
  
  for i in 1:NI
    for j in 0:NJ
      Jjfc[i, j, NK] = 0.5e0 * (Jac[i, j, NK] + Jac[i, j+1, NK])
      gj[i, j, NK, 1] = 0.5e0 * (g12[i, j] + g12[i, j+1]) * Jjfc[i, j, NK]
      gj[i, j, NK, 2] = 0.5e0 * (g22[i, j] + g22[i, j+1]) * Jjfc[i, j, NK]
      gqj[i, j, NK, 1] = qpr * gj[i, j, NK, 1]
      gqj[i, j, NK, 2] = qpr * gj[i, j, NK, 2]
    end
  end
  
  for j in 1:NJ
    for i in 0:NI
      gi3[i, j, NK] = 0.5e0 * (g13[i, j, NK] + g13[i+1, j, NK]) * Jifc[i, j, NK]
      gqi3[i, j, NK] = qpr * gi3[i, j, NK]
    end
  end

  for j in 0:NJ
    for i in 1:NI
      gj3[i, j, NK] = 0.5e0 * (g23[i, j, NK] + g23[i, j+1, NK]) * Jjfc[i, j, NK]
      gqj3[i, j, NK] = qpr * gj3[i, j, NK]
    end
  end
  
  for i in 0:NI
    for j in 0:NJ
      for k in NK:NK+1
        JacInv[i, j, k] = 1e0 / Jac[i, j, k]
      end
    end
  end

  
end