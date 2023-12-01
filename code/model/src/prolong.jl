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
                        local is = i * 2
                        fin[iw+1, js+1, kb+1] = fin[iw+1, js+1, kb+1] +
                                                an[1] * cor[i+1, j+1, k+1] + an[2] * cor[i+2, j+1, k+1] +
                                                an[3] * cor[i+2, j+2, k+1] + an[4] * cor[i+1, j+2, k+1] +
                                                an[5] * cor[i+1, j+1, k+2] + an[6] * cor[i+2, j+1, k+2] +
                                                an[7] * cor[i+2, j+2, k+2] + an[8] * cor[i+1, j+2, k+2]
                        fin[iw+2, js+1, kb+1] = fin[iw+2, js+1, kb+1] +
                                                bn[1] * cor[i+1, j+1, k+1] + bn[2] * cor[i+2, j+1, k+1] +
                                                bn[3] * cor[i+2, j+2, k+1] + bn[4] * cor[i+1, j+2, k+1] +
                                                bn[5] * cor[i+1, j+1, k+2] + bn[6] * cor[i+2, j+1, k+2] +
                                                bn[7] * cor[i+2, j+2, k+2] + bn[8] * cor[i+1, j+2, k+2]
                        fin[iw+2, js+2, kb+1] = fin[iw+2, js+2, kb+1] +
                                                cn[1] * cor[i+1, j+1, k+1] + cn[2] * cor[i+2, j+1, k+1] +
                                                cn[3] * cor[i+2, j+2, k+1] + cn[4] * cor[i+1, j+2, k+1] +
                                                cn[5] * cor[i+1, j+1, k+2] + cn[6] * cor[i+2, j+1, k+2] +
                                                cn[7] * cor[i+2, j+2, k+2] + cn[8] * cor[i+1, j+2, k+2]
                        fin[iw+1, js+2, kb+1] = fin[iw+1, js+2, kb+1] +
                                                dn[1] * cor[i+1, j+1, k+1] + dn[2] * cor[i+2, j+1, k+1] +
                                                dn[3] * cor[i+2, j+2, k+1] + dn[4] * cor[i+1, j+2, k+1] +
                                                dn[5] * cor[i+1, j+1, k+2] + dn[6] * cor[i+2, j+1, k+2] +
                                                dn[7] * cor[i+2, j+2, k+2] + dn[8] * cor[i+1, j+2, k+2]
                        fin[iw+1, js+1, kb+2] = fin[iw+1, js+1, kb+2] +
                                                en[1] * cor[i+1, j+1, k+1] + en[2] * cor[i+2, j+1, k+1] +
                                                en[3] * cor[i+2, j+2, k+1] + en[4] * cor[i+1, j+2, k+1] +
                                                en[5] * cor[i+1, j+1, k+2] + en[6] * cor[i+2, j+1, k+2] +
                                                en[7] * cor[i+2, j+2, k+2] + en[8] * cor[i+1, j+2, k+2]
                        fin[iw+2, js+1, kb+2] = fin[iw+2, js+1, kb+2] +
                                                fn[1] * cor[i+1, j+1, k+1] + fn[2] * cor[i+2, j+1, k+1] +
                                                fn[3] * cor[i+2, j+2, k+1] + fn[4] * cor[i+1, j+2, k+1] +
                                                fn[5] * cor[i+1, j+1, k+2] + fn[6] * cor[i+2, j+1, k+2] +
                                                fn[7] * cor[i+2, j+2, k+2] + fn[8] * cor[i+1, j+2, k+2]
                        fin[iw+2, js+2, kb+2] = fin[iw+2, js+2, kb+2] +
                                                gn[1] * cor[i+1, j+1, k+1] + gn[2] * cor[i+2, j+1, k+1] +
                                                gn[3] * cor[i+2, j+2, k+1] + gn[4] * cor[i+1, j+2, k+1] +
                                                gn[5] * cor[i+1, j+1, k+2] + gn[6] * cor[i+2, j+1, k+2] +
                                                gn[7] * cor[i+2, j+2, k+2] + gn[8] * cor[i+1, j+2, k+2]
                        fin[iw+1, js+2, kb+2] = fin[iw+1, js+2, kb+2] +
                                                hn[1] * cor[i+1, j+1, k+1] + hn[2] * cor[i+2, j+1, k+1] +
                                                hn[3] * cor[i+2, j+2, k+1] + hn[4] * cor[i+1, j+2, k+1] +
                                                hn[5] * cor[i+1, j+1, k+2] + hn[6] * cor[i+2, j+1, k+2] +
                                                hn[7] * cor[i+2, j+2, k+2] + hn[8] * cor[i+1, j+2, k+2]
                  end
            end
      end
end