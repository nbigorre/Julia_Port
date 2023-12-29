function uvchy(dtimel)
  local dte = dtimel / EPS
  local epsinv = 1e0 / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

  for j in 1:NJ
    for i in 1:NI
      local hxi = 0.5e0 * (h[i+1, j] - h[i-1, j])
      local heta = 0.5e0 * (h[i, j+1] - h[i, j-1])
      local hx = ux[i, j] * hxi + vx[i, j] * heta
      local hy = uy[i, j] * hxi + vy[i, j] * heta
      for k in 1:NK
        local pxi = 0.5e0 * (p[i+1, j, k] - p[i-1, j, k])
        local peta = 0.5e0 * (p[i, j+1, k] - p[i, j-1, k])
        if (j == 1)
          peta = p[i, j+1, k] - p[i, j, k]
        elseif (j == NJ)
          peta = p[i, j, k] - p[i, j-1, k]
        end
        local psig = 0.5e0 * (p[i, j, k+1] - p[i, j, k-1])
        local px = ux[i, j] * pxi + vx[i, j] * peta + wx[i, j, k] * psig
        local py = uy[i, j] * pxi + vy[i, j] * peta + wy[i, j, k] * psig

        cx[i, j, k] = cx[i, j, k] - dte * (gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * gradhn[i, j, 1]) + si[i, j, k] + @fortGet("qpr", rc_kind) * px)
        cy[i, j, k] = cy[i, j, k] - dte * (gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * gradhn[i, j, 2]) + sj[i, j, k] + @fortGet("qpr", rc_kind) * py)
      end
      gradhn[i, j, 1] = hx
      gradhn[i, j, 2] = hy
    end
  end

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local wfk = 0.5e0 * (czf[i, j, k] + czf[i, j, k-1])
        cz[i, j, k] = (wfk / Jac[i, j, k] - cx[i, j, k] * wx[i, j, k] - cy[i, j, k] * wy[i+1, j+1, k+1]) / (EPS * wz[i, j, k])
      end
    end
  end

end