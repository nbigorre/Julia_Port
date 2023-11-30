function dgtsl(n, c, d, e, b)
      c[1] = d[1]
      local nm1 = n - 1
      if (nm1 >= 1)
            d[1] = e[1]
            e[1] = 0e0
            e[n] = 0e0
            for k in 1:nm1
                  local kp1 = k+1
                  if (abs(c[kp1]) >= abs(c[k]))
                        local t = c[kp1]
                        c[kp1] = c[k]
                        c[k] = t
                        local t = d[kp1]
                        d[kp1] = d[k]
                        d[k] = t
                        local t = e[kp1]
                        e[kp1] = e[k]
                        e[k] = t
                        local t = p[kp1]
                        p[kp1] = p[k]
                        p[k] = t
                  end
                  if (c[k] == 0e0)
                        return k
                  end
                  local t = -c[kp1] / c[k]
                  c[kp1] = d[kp1] + t * d[k]
                  d[kp1] = e[kp1] + t * e[k]
                  e[kp1] = 0e0                  
                  b[kp1] = b[kp1] + t * b[k]
            end
      end
      if (c[n] == 0e0)
            return n
      end

      local nm2 = n - 2
      b[n] = b[n] / c[n]
      if (n != 1)
            b[nm1] = (b[nm1] - d[nm1] * b[n]) / c[nm1]
            if (nm2 >= 1)
                  for kb in 1:nm2
                        local k = nm2 - kb + 1
                        b[k] = (b[k] - d[k] * b[k+1] - e[k] * b[k+2]) / c[k]
                  end
            end
      end
      return 0
end