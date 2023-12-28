function vcenter(pf, dtimel, n)
  local dte = dtimel / EPS

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local pxi = 0.5e0 * (pf[i+2, j+1, k+1] - pf[i, j+1, k+1])
        local peta = 0.5e0 * (pf[i+1, j+2, k+1] - pf[i+1, j, k+1])
        local psig = 0.5e0 * (pf[i+1, j+1, k+2] - pf[i+1, j+1, k])
        local px = ux[i, j] * pxi + vx[i, j] * peta + wx[i+1, j+1, k+1] * psig
        local py = uy[i, j] * pxi + vy[i, j] * peta + wy[i+1, j+1, k+1] * psig
        local pz = wz[i+1, j+1, k+1] * psig

        u[i+1, j+1, k+1, n+1] = cx[i, j, k] - dte * (@fortGet("qpr", rc_kind) * px + si[i, j, k])
        v[i+1, j+1, k+1, n+1] = cy[i, j, k] - dte * (@fortGet("qpr", rc_kind) * py + sj[i, j, k])

      end
    end
  end

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local wfk = 0.5e0 * (wf[i, j, k] + wf[i, j, k-1])
        w[i+1, j+1, k+1, n+1] = (wfk / Jac[i, j, k] - u[i+1, j+1, k+1, n+1] * wx[i+1, j+1, k+1] - v[i+1, j+1, k+1, n+1] * wy[i+1, j+1, k+1]) / (EPS * wz[i+1, j+1, k+1])
      end
    end
  end
end