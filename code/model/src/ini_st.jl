function ini_st() 
  local npl = 33
  local depoff = zeros(rc_kind, npl)
  local dep = [-5500e0, -5000e0, -4500e0, -4000e0, -3500e0,
    -3000e0, -2500e0, -2000e0, -1750e0, -1500e0, -1400e0,
    -1300e0, -1200e0, -1100e0, -1000e0, -900e0, -800e0,
    -700e0, -600e0, -500e0, -400e0, -300e0, -250e0,
    -200e0, -150e0, -125e0, -100e0, -75e0, -50e0,
    -30e0, -20e0, -10e0, -0e0]

  local svert = [34.76, 34.87, 34.91, 34.94, 35.00, 35.05, 35.08,
    35.11, 35.13, 35.13, 35.11, 35.05, 34.99, 34.94,
    34.89, 34.91, 34.89, 34.89, 34.88, 34.89, 34.90,
    34.88, 34.88, 34.88, 34.88, 34.88, 34.87, 34.81,
    34.76, 34.75, 34.75, 34.74, 34.71]

  local tvert = [1.61, 1.72, 1.77, 1.78, 1.87, 2.18, 2.61, 3.05,
    3.28, 3.62, 3.80, 4.06, 4.32, 4.57, 4.84, 5.29, 5.71,
    6.37, 6.96, 7.83, 8.75, 9.82, 10.42, 11.00, 11.71, 12.12,
    12.63, 13.28, 14.46, 16.01, 16.76, 17.04, 17.25]

  global mldepth = 50e0
  local dsdz = 0e0

  for k in 1:npl
    depoff[k] = dep[k] - mldepth
  end

  for j in 0:NJ+1
    for i in 0:NI+1
      for k in 0:NK+1
        local z = DL * zc[i, j, k]

        T[i, j, k, 0] = 29.95e0 + z * 0.0454e0
        s[i, j, k, 0] = 35e0

      end
    end
  end


  for k in 0:NK+1
    for i in 1:NI
      s[i, 0, k, 0] = s[i, 1, k, 0]
      s[i, NJ+1, k, 0] = s[i, NJ, k, 0]
      T[i, 0, k, 0] = T[i, 1, k, 0]
      T[i, NJ+1, k, 0] = T[i, NJ, k, 0]
    end
    for j in 0:NJ+1
      s[0, j, k, 0] = s[NI, j, k, 0]
      s[NI+1, j, k, 0] = s[1, j, k, 0]
      T[0, j, k, 0] = T[NI, j, k, 0]
      T[NI+1, j, k, 0] = T[1, j, k, 0]
    end
  end
end
