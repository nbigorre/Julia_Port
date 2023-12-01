function resid(m, nxm, nym, nzm, cp, p, fn, res)
      local maxres = 0e0

      for k in 1:nzm
            for j in 1:nym
                  for i in 1:nxm
                        res[i, j, k] = fn[i, j, k] -
                                       (cp[1, i, j, k] * p[i+1, j+1, k+1]
                                        + cp[2, i, j, k] * p[i+2, j+1, k+1]
                                        + cp[3, i, j, k] * p[i+1, j+2, k+1]
                                        + cp[4, i, j, k] * p[i, j+1, k+1]
                                        + cp[5, i, j, k] * p[i+1, j, k+1]
                                        + cp[6, i, j, k] * p[i+1, j+1, k+2]
                                        + cp[7, i, j, k] * p[i+1, j+1, k]
                                        + cp[8, i, j, k] * p[i, j+2, k+1]
                                        + cp[9, i, j, k] * p[i, j, k+1]
                                        + cp[10, i, j, k] * p[i+2, j, k+1]
                                        + cp[11, i, j, k] * p[i+2, j+2, k+1]
                                        + cp[12, i, j, k] * p[i, j+1, k]
                                        + cp[13, i, j, k] * p[i+2, j+1, k]
                                        + cp[14, i, j, k] * p[i+2, j+1, k+2]
                                        + cp[15, i, j, k] * p[i, j+1, k+2]
                                        + cp[16, i, j, k] * p[i+1, j, k]
                                        + cp[17, i, j, k] * p[i+1, j+2, k]
                                        + cp[18, i, j, k] * p[i+1, j+2, k+2]
                                        + cp[19, i, j, k] * p[i+1, j, k+2])
                        maxres = max(maxres, abs(res[i, j, k]))
                  end
            end
      end

      if (maxres > 3000e0)
            println("STOP. res too large, i,j,k,maxres=",maxres)
            exit(1)
      end

      return maxres
end