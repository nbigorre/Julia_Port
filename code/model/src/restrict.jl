function restrict(nxm, nym, nzm, fin, cor)
      local coeff = 0.125 * 4e0
      for k in 1:div(nzm, 2)
            local kk = 2 * k - 1
            for j in 1:div(nym, 2)
                  local jj = 2 * j - 1
                  for i in 1:div(nxm, 2)
                        local ii = 2 * i - 1
                        cor[i,j,k] = coeff*(fin[ii,jj,kk]+fin[ii,jj+1,kk]
                                      +fin[ii+1,jj,kk]+fin[ii+1,jj+1,kk] +
                                      fin[ii,jj,kk+1] + fin[ii,jj+1,kk+1]
                                      + fin[ii+1,jj,kk+1] + fin[ii+1,jj+1,kk+1] )         
                  end
            end
      end
end