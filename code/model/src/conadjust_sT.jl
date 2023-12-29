function conadjust(stepl, n)
  for j in 1:NJ
    for i in 1:NI
      for k in NK-1:-1:1
        conv[i, j, k] = 0e0
        local rhtop = rho[i, j, k+1]
        local dzupper = zf[i, j, k+1] - zf[i, j, k]

        local rh = rho[i, j, k]
        local dz = zf[i, j, k] - zf[i, j, k-1]

        if (rh < rhtop)
          conv[i, j, k] = 1e0
          local zinv = 1e0 / (dzupper + dz)
          s[i, j, k, n] = (dzupper * s[i, j, k+1, n] + dz * s[i, j, k, n]) * zinv
          s[i, j, k+1, n] = s[i, j, k, n]
          T[i, j, k, n] = (dzupper * T[i, j, k+1, n] + dz * T[i, j, k, n]) * zinv
          T[i, j, k+1, n] = T[i, j, k, n]
          rho[i, j, k] = potdens(s[i, j, k, n], T[i, j, k, n])
          rho[i, j, k+1] = rho[i, j, k]
          for it in 1:ntr
            Tr[it, i, j, k, n] = (dzupper * Tr[it, i, j, k, n] + dz * Tr[it, i, j, k, n]) * zinv
            Tr[it, i, j+1, k, n] = Tr[it, i, j, k, n]
          end
        end
      end
    end
  end

  if (mod(stepl - 1, 100) == 0)
    con100 .= 0
  else
    con100 .+= conv
  end
end