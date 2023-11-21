function vface(pf, dtimel::rc_kind)
  pf = reshape(view(pf, 1:((NI+2)*(NJ+2)*(NK+2))), NI+2, NJ+2, NK+2)
  local dte = dtimel / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)
  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
          local px = ((pf[i+2, j+1, k+1] - pf[i+1, j+1, k+1]) * gqi[i+1, j, k, 1]
          + 0.25e0 * (pf[i+2, j+2, k+1] + pf[i+1, j+2, k+1] - pf[i+2, j, k+1] - pf[i+1, j, k+1]) * gqi[i+1, j, k, 2]
           + 0.25e0 * (pf[i+2, j+1, k+2] + pf[i+1, j+1, k+2] - pf[i+2, j+1, k] - pf[i+1, j+1, k]) * gqi3[i+1, j, k])
           uf[i+1, j, k] = cxf[i+1, j, k] - dte * (px + sifc[i+1, j, k])
      end
    end
  end
  
  uf[1,:,:] .= uf[NI+1,:,:]

  for i in 1:NI
    for j in 0:NJ
      for k in 1:NK
        if (j == 0)
          vf[i, j+1, k] = vfbcs[i, k]
        elseif (j == NJ)
          vf[i, j+1, k] = vfbcn[i, k]
        else 
          local py = ((pf[i+1, j+2, k+1] - pf[i+1, j+1, k+1]) * gqj[i, j+1, k, 2]
            + 0.25e0 * (pf[i+2, j+2, k+1] + pf[i+2, j+1, k+1] - pf[i, j+2, k+1] - pf[i, j+1, k+1]) * gqj[i, j+1, k, 1]
            + 0.25e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+1, k+2] - pf[i+1, j+2, k] - pf[i+1, j+1, k]) * gqj3[i, j+1, k])
          vf[i, j+1, k] = cyf[i, j+1, k] - dte * (py + sjfc[i, j+1, k])
        end
      end
    end
  end

  for j in 1:NJ
    for i in 1:NI
      wf[i, j, 1] = wfbcb[i, j]
      for k in 1:NK
        local pz = ((pf[i+1, j+1, k+2] - pf[i+1, j+1, k+1]) * gqk[i+1, j+1, k+1, 3]
          + 0.25e0 * (pf[i+2, j+1, k+2] + pf[i+2, j+1, k+1] - pf[i, j+1, k+2] - pf[i, j+1, k+1]) * gqk[i+1, j+1, k+1, 1]
          + 0.25e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+2, k+1] - pf[i+1, j, k+2] - pf[i+1, j, k+1]) * gqk[i+1, j+1, k+1, 2])
        wf[i, j, k+1] = czf[i, j, k+1] - dte * (pz + skfc[i+1, j+1, k+1])
      end
    end
  end

end