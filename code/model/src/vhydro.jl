function vhydro(dtimel)
  local dte = dtimel / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

  for j in 1:NJ
    for i in 1:NI
      local hxi = h[i+2, j+1] - h[i+1, j+1]
      local heta = 0.25e0 * ( h[i+2, j+2] + h[i+1, j+2] - h[i+2, j] - h[i+1, j])
      for k in 1:NK
        local hx = gi[i+1, j, k, 1] * hxi + gi[i+1, j, k, 2] * heta
        local gradh = gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * hxn[i+1, j, k])
        cxf[i+1, j, k] -= dte * (gradh + sifc[i+1, j, k])
        hxn[i+1, j, k] = hx
      end
    end

    for k in 1:NK
      cxf[1, j, k] = cxf[NI+1, j, k]
      hxn[1, j, k] = hxn[NI+1, j, k]
    end
  end

  for i in 1:NI
    for j in 0:NJ
      local hxi = 0.25e0 * (h[i+2, j+2] + h[i+2, j+1] - h[i, j+2] - h[i, j+1])
      local heta = h[i+1, j+2] - h[i+1, j+1]
      for k in 1:NK
        local hy = gj[i, j+1, k, 1] * hxi + gj[i, j+1, k, 2] * heta
        local gradh = gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * hyn[i, j+1, k])
        cyf[i, j+1, k] -= dte * (gradh + sjfc[i, j+1, k])
        hyn[i, j+1, k] = hy
      end
    end
  end

  for i in 1:NI
    for k in 1:NK
      cyf[i, NJ+1, k] = vfbcn[i,k]
      cyf[i, 1, k] = vfbcs[i,k]
    end
  end

  for k in 0:NK
    for j in 1:NJ
      for i in 1:NI
        local pz = (p[i+1, j+1, k+2] - p[i+1, j+1, k+1]) * gqk[i+1, j+1, k+1, 3] + 0.25e0 * (p[i+2, j+1, k+2] + p[i+2, j+1, k+1] - p[i, j+1, k+2] - p[i, j+1, k+1]) * gqk[i+1, j+1, k+1, 1] + 0.25e0 * (p[i+1, j+2, k+2] + p[i+1, j+2, k+1] - p[i+1, j, k+2] - p[i+1, j, k+1]) * gqk[i+1, j+1, k+1, 2] 
        czf[i, j, k+1] -= dte * (pz + skfc[i+1, j+1, k+1])
      end
    end
  end

  @ccall "./PSOM_LIB.so".uvchy_(Ref(dtimel)::Ref{rc_kind})::Cvoid
end