function vcenter(pf, dtimel, n)
  local dte = dtimel / EPS

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local pxi = 0.5e0 * (pf[i+1, j, k] - pf[i-1, j, k])
        local peta = 0.5e0 * (pf[i, j+1, k] - pf[i, j-1, k])
        local psig = 0.5e0 * (pf[i, j, k+1] - pf[i, j, k-1])
        local px = ux[i, j] * pxi + vx[i, j] * peta + wx[i, j, k] * psig
        local py = uy[i, j] * pxi + vy[i, j] * peta + wy[i, j, k] * psig
        local pz = wz[i, j, k] * psig

        u[i, j, k, n] = cx[i, j, k] - dte * (@fortGet("qpr", rc_kind) * px + si[i, j, k])
        v[i, j, k, n] = cy[i, j, k] - dte * (@fortGet("qpr", rc_kind) * py + sj[i, j, k])

      end
    end
  end

  for j in 1:NJ
    for i in 1:NI
      for k in 1:NK
        local wfk = 0.5e0 * (wf[i, j, k] + wf[i, j, k-1])
        w[i, j, k, n] = (wfk / Jac[i, j, k] - u[i, j, k, n] * wx[i, j, k] - v[i, j, k, n] * wy[i, j, k]) / (EPS * wz[i, j, k])
      end
    end
  end
end