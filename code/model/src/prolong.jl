function findn(x, y, z)
      local pn1 = (1e0 - x) * (1e0 - y) * (1e0 - z)
      local pn2 = x * (1e0 - y) * (1e0 - z)
      local pn3 = x * y * (1e0 - z)
      local pn4 = (1e0 - x) * y * (1e0 - z)
      local pn5 = (1e0 - x) * (1e0 - y) * z
      local pn6 = x * (1e0 - y) * z
      local pn7 = x * y * z
      local pn8 = (1e0 - x) * y * z
      local sum = pn1 + pn2 + pn3 + pn4 + pn5 + pn6 + pn7 + pn8
      if (abs(sum - 1e0) > 0.0001)
            println["coefficients wrong"]
      end
      return (pn1, pn2, pn3, pn4, pn5, pn6, pn7, pn8)
end


function prolong(nxm, nym, nzm, cor, fin)
      local x = 0.25e0
      local y = 0.25e0
      local z = 0.25e0
      local an = findn(x, y, z)
      x = x + 5e0
      local bn = findn(x, y, z)
      y = y + 0.5e0
      local cn = findn(x, y, z)
      x = x - 0.5e0
      local dn = findn(x, y, z)
      y = y - 0.5e0
      z = z + 0.5e0
      local en = findn(x, y, z)
      x = x + 0.5e0
      local fn = findn(x, y, z)
      y = y + 0.5e0
      local gn = findn(x, y, z)
      x = x - 0.5e0
      local hn = findn(x, y, z)

      for k in 0:nzm
            local kb = k * 2
            for j in 0:nym
                  local js = j * 2
                  for i in 0:nxm
                        local iw = i * 2
                        fin[iw, js, kb] = fin[iw, js, kb] +
                                                an[1] * cor[i, j, k] + an[2] * cor[i+1, j, k] +
                                                an[3] * cor[i+1, j+1, k] + an[4] * cor[i, j+1, k] +
                                                an[5] * cor[i, j, k+1] + an[6] * cor[i+1, j, k+1] +
                                                an[7] * cor[i+1, j+1, k+1] + an[8] * cor[i, j+1, k+1]
                        fin[iw+1, js, kb] = fin[iw+1, js, kb] +
                                                bn[1] * cor[i, j, k] + bn[2] * cor[i+1, j, k] +
                                                bn[3] * cor[i+1, j+1, k] + bn[4] * cor[i, j+1, k] +
                                                bn[5] * cor[i, j, k+1] + bn[6] * cor[i+1, j, k+1] +
                                                bn[7] * cor[i+1, j+1, k+1] + bn[8] * cor[i, j+1, k+1]
                        fin[iw+1, js+1, kb] = fin[iw+1, js+1, kb] +
                                                cn[1] * cor[i, j, k] + cn[2] * cor[i+1, j, k] +
                                                cn[3] * cor[i+1, j+1, k] + cn[4] * cor[i, j+1, k] +
                                                cn[5] * cor[i, j, k+1] + cn[6] * cor[i+1, j, k+1] +
                                                cn[7] * cor[i+1, j+1, k+1] + cn[8] * cor[i, j+1, k+1]
                        fin[iw, js+1, kb] = fin[iw, js+1, kb] +
                                                dn[1] * cor[i, j, k] + dn[2] * cor[i+1, j, k] +
                                                dn[3] * cor[i+1, j+1, k] + dn[4] * cor[i, j+1, k] +
                                                dn[5] * cor[i, j, k+1] + dn[6] * cor[i+1, j, k+1] +
                                                dn[7] * cor[i+1, j+1, k+1] + dn[8] * cor[i, j+1, k+1]
                        fin[iw, js, kb+1] = fin[iw, js, kb+1] +
                                                en[1] * cor[i, j, k] + en[2] * cor[i+1, j, k] +
                                                en[3] * cor[i+1, j+1, k] + en[4] * cor[i, j+1, k] +
                                                en[5] * cor[i, j, k+1] + en[6] * cor[i+1, j, k+1] +
                                                en[7] * cor[i+1, j+1, k+1] + en[8] * cor[i, j+1, k+1]
                        fin[iw+1, js, kb+1] = fin[iw+1, js, kb+1] +
                                                fn[1] * cor[i, j, k] + fn[2] * cor[i+1, j, k] +
                                                fn[3] * cor[i+1, j+1, k] + fn[4] * cor[i, j+1, k] +
                                                fn[5] * cor[i, j, k+1] + fn[6] * cor[i+1, j, k+1] +
                                                fn[7] * cor[i+1, j+1, k+1] + fn[8] * cor[i, j+1, k+1]
                        fin[iw+1, js+1, kb+1] = fin[iw+1, js+1, kb+1] +
                                                gn[1] * cor[i, j, k] + gn[2] * cor[i+1, j, k] +
                                                gn[3] * cor[i+1, j+1, k] + gn[4] * cor[i, j+1, k] +
                                                gn[5] * cor[i, j, k+1] + gn[6] * cor[i+1, j, k+1] +
                                                gn[7] * cor[i+1, j+1, k+1] + gn[8] * cor[i, j+1, k+1]
                        fin[iw, js+1, kb+1] = fin[iw, js+1, kb+1] +
                                                hn[1] * cor[i, j, k] + hn[2] * cor[i+1, j, k] +
                                                hn[3] * cor[i+1, j+1, k] + hn[4] * cor[i, j+1, k] +
                                                hn[5] * cor[i, j, k+1] + hn[6] * cor[i+1, j, k+1] +
                                                hn[7] * cor[i+1, j+1, k+1] + hn[8] * cor[i, j+1, k+1]
                  end
            end
      end
end