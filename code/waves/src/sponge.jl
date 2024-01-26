function sponge(n, dum_time)
  #
  #  USE header
  #
  #  implicit none
  #
  #  REAL(kind=rc_kind) ::  dum_time, gamma_T, dummy
  #  REAL(kind=rc_kind) ::  rinn, rout, rrr, rdum, s_relax, u_relax, v_relax
  #  INTEGER, PARAMETER ::  spl=15, splz=60
  #  
  #  INTEGER :: i,j,k,n
  #
  local gamma_T = 1e0 / (dum_time)

  for j in 0:spl+1

    dummy = gamma_T * cos((j - 1) / spl * PI / 2)
    if (j == 0)
      dummy = gamma_T
    end

    for k in splz+1:NK+1
      for i in 0:NI+1
        s[i, j, k, n] = s[i, j, k, n] - dum_time * dummy * (s[i, j, k, n] - rho_refS[k])
        u[i, j, k, n] = u[i, j, k, n] - dum_time * dummy * u[i, j, k, n]
        v[i, j, k, n] = v[i, j, k, n] - dum_time * dummy * v[i, j, k, n]
        w[i, j, k, n] = w[i, j, k, n] - dum_time * dummy * w[i, j, k, n]
      end
    end

  end

  for j in NJ-spl:NJ+1

    local dummy = gamma_T * sin((j - (NJ - spl)) / spl * PI / 2)
    if (j == NJ + 1)
      dummy = gamma_T
    end
    for k in splz+1:NK+1
      for i in 0:NI+1
        s[i, j, k, n] = s[i, j, k, n] - dum_time * dummy * (s[i, j, k, n] - rho_refN[k])
        u[i, j, k, n] = u[i, j, k, n] - dum_time * dummy * u[i, j, k, n]
        v[i, j, k, n] = v[i, j, k, n] - dum_time * dummy * v[i, j, k, n]
        w[i, j, k, n] = w[i, j, k, n] - dum_time * dummy * w[i, j, k, n]
      end
    end

  end

  local rinn = 2e0 * PI / (ffc[NI/2, NJ/2] * FPAR)
  local rout = rinn / 1e3

  for k in 0:splz

    rdum = (rc_kind(splz - k) * rout + rc_kind(k) * rinn) / rc_kind(splz)
    rrr = 0e0
    if (rdum != 0)
      rrr = 1e0 / rdum

    end

    for j in 0:NJ+1
      for i in 0:NI+1

        s_relax = (rc_kind(splz - k) * rho_refB(j) + rc_kind(k) * s[i, j, k, n]) / rc_kind(splz)
        u_relax = (rc_kind(k) * u[i, j, k, n]) / rc_kind(splz)
        v_relax = (rc_kind(k) * v[i, j, k, n]) / rc_kind(splz)

        s[i, j, k, n] = s[i, j, k, n] - dum_time * rrr * (s[i, j, k, n] - s_relax)
        u[i, j, k, n] = u[i, j, k, n] - dum_time * rrr * (u[i, j, k, n] - u_relax)
        v[i, j, k, n] = v[i, j, k, n] - dum_time * rrr * (v[i, j, k, n] - v_relax)

      end
    end

  end

end
