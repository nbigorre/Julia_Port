

function coriolis(n)
  local fac = EPS * delta
  local fac2 = EPS * lambda
  local ainv = 1e0 / apr
  local betainv = 1e0 / beta

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        si[i, j, k] = -ffc[i, j] * v[i, j, k, n] + fac * bbc[i, j] * w[i, j, k, n] + drpx[i, j, k] - EPS * uvis[i, j, k]
        sj[i, j, k] = ffc[i, j] * u[i, j, k, n] + drpy[i, j, k] - EPS * vvis[i, j, k]
        sk[i, j, k] = fnhhy * (-bbc[i, j] * u[i, j, k, n] - fac2 * (u[i, j, k, n] * u[i, j, k, n] + v[i, j, k, n] * v[i, j, k, n]) * ainv) - betainv * wvis[i, j, k]
      end
    end
  end
end