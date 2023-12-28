

function coriolis(n)
  local fac = EPS * @fortGet("delta", rc_kind)
  local fac2 = EPS * @fortGet("lambda", rc_kind)
  local ainv = 1e0 / apr
  local betainv = 1e0 / @fortGet("beta", rc_kind)

  for k in 1:NK
    for j in 1:NJ
      for i in 1:NI
        si[i, j, k] = -ffc[i, j] * v[i+1, j+1, k+1, n+1] + fac * bbc[i, j] * w[i+1, j+1, k+1, n+1] + drpx[i,j,k] - EPS * uvis[i,j,k]
        sj[i, j, k] =  ffc[i, j] * u[i+1, j+1, k+1, n+1]                                               + drpy[i,j,k] - EPS * vvis[i,j,k]
        sk[i, j, k] = @fortGet("fnhhy", rc_kind) * (-bbc[i, j] * u[i+1, j+1, k+1, n+1] - fac2 * (u[i+1, j+1, k+1, n+1] * u[i+1, j+1, k+1, n+1] + v[i+1, j+1, k+1, n+1] * v[i+1, j+1, k+1, n+1]) * ainv) - betainv * wvis[i,j,k]
      end
    end
  end
end