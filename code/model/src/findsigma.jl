function findsigma(i, j, z)

  local zs = z
  local sigma = 0
  local sigup = 0e0

  if (zc[i, j, 1] * DL <= z)

    for k in 2:NK
      if (zc[i, j, k] * DL >= z)
        local dsig = (zc[i, j, k] - zc[i, j, k-1]) * DL
        sigup = -(zc[i, j, k-1] * DL - z) / dsig
        local sigdn = (zc[i, j, k] * DL - z) / dsig
        sigma = k
        # zs[i,j]      = sigup*zc[i,j,k]     + sigdn*zc[i,j,k-1] 
        return (sigma, sigup)
      end
    end

  else
    sigma = 1
    sigup = 0e0

  end
  return (sigma, sigup)


end
