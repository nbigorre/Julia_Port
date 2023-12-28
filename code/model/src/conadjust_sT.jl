function conadjust(stepl, n)
  for j in 1:NJ
    for i in 1:NI
      for k in NK-1:-1:1
        conv[i+1, j+1, k+1] = 0e0
        local rhtop = rho[i, j, k+1]
        local dzupper = zf[i+1, j+1, k+3] - zf[i+1, j+1, k+2]

        local rh = rho[i, j, k]
        local dz = zf[i + 1, j + 1, k + 2] - zf[i + 1, j + 1, k + 1]

        if (rh < rhtop)
          conv[i+1, j+1, k+1] = 1e0
          local zinv = 1e0 / (dzupper + dz)
          s[i+1, j+1, k+1, n+1] = (dzupper * s[i+1, j+1, k+2, n+1] + dz * s[i+1, j+1, k+1, n+1]) * zinv
          s[i+1, j+1, k+2, n+1] = s[i+1, j+1, k+1, n+1]
          T[i+1, j+1, k+1, n+1] = (dzupper * T[i+1, j+1, k+2, n+1] + dz * T[i+1, j+1, k+1, n+1]) * zinv
          T[i+1, j+1, k+2, n+1] = T[i+1, j+1, k+1, n+1]
          rho[i, j, k] = potdens(s[i+1, j+1, k+1, n+1], T[i+1, j+1, k+1, n+1])
          rho[i, j, k+1] = rho[i, j, k]
          for it in 1:ntr
            Tr[it, i+1, j+1, k+1, n+1] = (dzupper * Tr[it, i+1, j+1, k+1, n+1] + dz * Tr[it, i+1, j+1, k+1, n+1]) * zinv
            Tr[it, i+1, j+1, k+2, n+1] = Tr[it, i+1, j+1, k+1, n+1]
          end
        end
      end
    end
  end

  if (mod(stepl-1, 100) == 0)
    con100 .= 0
  else
    con100 .+= conv
  end
end