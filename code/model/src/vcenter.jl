function vcenter(pf, dtimel, n)
  local dte = dtimel / EPS

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local pxi = 0.5e0 * (pf[i+2, j+1, k+1] - pf[i, j+1, k+1])
        local peta = 0.5e0 * (pf[i+1, j+2, k+1] - pf[i+1, j, k+1])
        local psig = 0.5e0 * (pf[i+1, j+1, k+2] - pf[i+1, j+1, k])
        local px = ux[i+1, j+1] * pxi + vx[i+1, j+1] * peta + wx[i+1, j+1, k+1] * psig
        local py = uy[i+1, j+1] * pxi + vy[i+1, j+1] * peta + wy[i+1, j+1, k+1] * psig
        local pz = wz[i+1, j+1, k+1] * psig

        u[i+1, j+1, k+1, n+1] = cx[i+1, j+1, k+1] - dte * (@fortGet("qpr", rc_kind) * px + si[i+1, j+1, k+1])
        v[i+1, j+1, k+1, n+1] = cy[i+1, j+1, k+1] - dte * (@fortGet("qpr", rc_kind) * py + sj[i+1, j+1, k+1])

      end
    end
  end

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local wfk = 0.5e0 * (wf[i, j, k+1] + wf[i, j, k])
        w[i+1, j+1, k+1, n+1] = (wfk / Jac[i+1, j+1, k+1] - u[i+1, j+1, k+1, n+1] * wx[i+1, j+1, k+1] - v[i+1, j+1, k+1, n+1] * wy[i+1, j+1, k+1]) / (EPS * wz[i+1, j+1, k+1])
      end
    end
  end
end