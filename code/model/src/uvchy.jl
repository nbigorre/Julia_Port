function uvchy(dtimel)
  local dte = dtimel / EPS
  local epsinv = 1e0 / EPS
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

  for j in 1:NJ
    for i in 1:NI
      local hxi = 0.5e0 * (h[i+2, j+1] - h[i, j+1])
      local heta = 0.5e0 * (h[i+1, j+2] - h[i+1, j])
      local hx = ux[i+1, j+1] * hxi + vx[i+1, j+1] * heta
      local hy = uy[i+1, j+1] * hxi + vy[i+1, j+1] * heta
      for k in 1:NK
        local pxi = 0.5e0 * (p[i+2, j+1, k+1] - p[i, j+1, k+1])
        local peta = 0.5e0 * (p[i+1, j+2, k+1] - p[i+1, j, k+1])
        if (j == 1)
          peta = p[i+1, j+2, k+1] - p[i+1, j+1, k+1]
        elseif (j == NJ)
          peta = p[i+1, j+1, k+1] - p[i+1, j, k+1]
        end
        local psig = 0.5e0 * (p[i+1, j+1, k+2] - p[i+1, j+1, k])
        local px = ux[i+1, j+1] * pxi + vx[i+1, j+1] * peta + wx[i+1, j+1, k+1] * psig
        local py = uy[i+1, j+1] * pxi + vy[i+1, j+1] * peta + wy[i+1, j+1, k+1] * psig

        cx[i+1, j+1, k+1] = cx[i+1, j+1, k+1] - dte * (gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * gradhn[i+1, j+1, 1]) + si[i+1, j+1, k+1] + @fortGet("qpr", rc_kind) * px)
        cy[i+1, j+1, k+1] = cy[i+1, j+1, k+1] - dte * (gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * gradhn[i+1, j+1, 2]) + sj[i+1, j+1, k+1] + @fortGet("qpr", rc_kind) * py)
      end
      gradhn[i+1, j+1, 1] = hx
      gradhn[i+1, j+1, 2] = hy
    end
  end

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local wfk = 0.5e0 * (czf[i, j, k+1] + czf[i, j, k])
        cz[i+1, j+1, k+1] = (wfk / Jac[i+1, j+1, k+1] - cx[i+1, j+1, k+1] * wx[i+1, j+1, k+1] - cy[i+1, j+1, k+1] * wy[i+1, j+1, k+1]) / (EPS * wz[i+1, j+1, k+1])
      end
    end
  end

end