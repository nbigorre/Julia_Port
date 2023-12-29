function vface(pf, dtimel::rc_kind)
  pf = OffsetArrays.reshape(view(pf, 1:((NI+2)*(NJ+2)*(NK+2))), 0:NI+1, 0:NJ+1, 0:NK+1)
  local dte = dtimel / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)
  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
          local px = ((pf[i+1, j, k] - pf[i, j, k]) * gqi[i, j, k, 1]
          + 0.25e0 * (pf[i+1, j+1, k] + pf[i, j+1, k] - pf[i+1, j-1, k] - pf[i, j-1, k]) * gqi[i, j, k, 2]
           + 0.25e0 * (pf[i+1, j, k+1] + pf[i, j, k+1] - pf[i+1, j, k-1] - pf[i, j, k-1]) * gqi3[i, j, k])
           uf[i, j, k] = cxf[i, j, k] - dte * (px + sifc[i, j, k])
      end
    end
  end
  
  uf[0,:,:] .= uf[NI,:,:]

  for i in 1:NI
    for j in 0:NJ
      for k in 1:NK
        if (j == 0)
          vf[i, j, k] = vfbcs[i, k]
        elseif (j == NJ)
          vf[i, j, k] = vfbcn[i, k]
        else 
          local py = ((pf[i, j+1, k] - pf[i, j, k]) * gqj[i, j, k, 2]
            + 0.25e0 * (pf[i+1, j+1, k] + pf[i+1, j, k] - pf[i-1, j+1, k] - pf[i-1, j, k]) * gqj[i, j, k, 1]
            + 0.25e0 * (pf[i, j+1, k+1] + pf[i, j, k+1] - pf[i, j+1, k-1] - pf[i, j, k-1]) * gqj3[i, j, k])
          vf[i, j, k] = cyf[i, j, k] - dte * (py + sjfc[i, j, k])
        end
      end
    end
  end

  for j in 1:NJ
    for i in 1:NI
      wf[i, j, 0] = wfbcb[i, j]
      for k in 1:NK
        local pz = ((pf[i, j, k+1] - pf[i, j, k]) * gqk[i, j, k, 3]
          + 0.25e0 * (pf[i+1, j, k+1] + pf[i+1, j, k] - pf[i-1, j, k+1] - pf[i-1, j, k]) * gqk[i, j, k, 1]
          + 0.25e0 * (pf[i, j+1, k+1] + pf[i, j+1, k] - pf[i, j-1, k+1] - pf[i, j-1, k]) * gqk[i, j, k, 2])
        wf[i, j, k] = czf[i, j, k] - dte * (pz + skfc[i, j, k])
      end
    end
  end

end