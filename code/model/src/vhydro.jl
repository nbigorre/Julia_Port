function vhydro(dtimel)
  local dte = dtimel / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

  for j in 1:NJ
    for i in 1:NI
      local hxi = h[i+1, j] - h[i, j]
      local heta = 0.25e0 * ( h[i+1, j+1] + h[i, j+1] - h[i+1, j-1] - h[i, j-1])
      for k in 1:NK
        local hx = gi[i, j, k, 1] * hxi + gi[i, j, k, 2] * heta
        local gradh = gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * hxn[i, j, k])
        cxf[i, j, k] -= dte * (gradh + sifc[i, j, k])
        hxn[i, j, k] = hx
      end
    end

    for k in 1:NK
      cxf[0, j, k] = cxf[NI, j, k]
      hxn[0, j, k] = hxn[NI, j, k]
    end
  end

  for i in 1:NI
    for j in 0:NJ
      local hxi = 0.25e0 * (h[i+1, j+1] + h[i+1, j] - h[i-1, j+1] - h[i-1, j])
      local heta = h[i, j+1] - h[i, j]
      for k in 1:NK
        local hy = gj[i, j, k, 1] * hxi + gj[i, j, k, 2] * heta
        local gradh = gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * hyn[i, j, k])
        cyf[i, j, k] -= dte * (gradh + sjfc[i, j, k])
        hyn[i, j, k] = hy
      end
    end
  end

  for i in 1:NI
    for k in 1:NK
      cyf[i, NJ, k] = vfbcn[i,k]
      cyf[i, 0, k] = vfbcs[i,k]
    end
  end

  for k in 0:NK
    for j in 1:NJ
      for i in 1:NI
        local pz = (p[i, j, k+1] - p[i, j, k]) * gqk[i, j, k, 3] + 0.25e0 * (p[i+1, j, k+1] + p[i+1, j, k] - p[i-1, j, k+1] - p[i-1, j, k]) * gqk[i, j, k, 1] + 0.25e0 * (p[i, j+1, k+1] + p[i, j+1, k] - p[i, j-1, k+1] - p[i, j-1, k]) * gqk[i, j, k, 2] 
        czf[i, j, k] -= dte * (pz + skfc[i, j, k])
      end
    end
  end

  uvchy(dtimel)
  #@ccall "./PSOM_LIB.so".uvchy_(Ref(dtimel)::Ref{rc_kind})::Cvoid
end