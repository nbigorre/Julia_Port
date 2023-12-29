function cdiv(dtimel, n)
    local maxdiv = 0e0
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                local Udxi = 0.5e0 * ((u[i+1, j, k, n] * ux[i+1, j] + v[i+1, j, k, n] * uy[i+1, j]) * Jac[i+1, j, k] 
                                    - (u[i-1, j, k, n] * ux[i-1, j] + v[i-1, j, k, n] * uy[i-1, j]) * Jac[i-1, j, k])
                local Vdelta = 0.5e0 * ((u[i, j+1, k, n] * vx[i, j+1] + v[i, j+1, k, n] * vy[i, j+1]) * Jac[i, j+1, k]
                                    - (u[i, j-1, k, n] * vx[i, j-1] + v[i, j-1, k, n] * vy[i, j-1]) * Jac[i, j-1, k])
                local Wdsig = 0.5e0 * ((u[i, j, k+1, n] * wx[i, j, k+1] + v[i, j, k+1, n] * wy[i, j, k+1] + EPS * w[i, j, k+1, n] * wz[i, j, k]) * Jac[i, j, k+1]
                                    - (u[i, j, k-1, n] * wx[i, j, k-1] + v[i, j, k-1, n] * wy[i, j, k-1] + EPS * w[i, j, k-1, n] * wz[i, j, k-1]) * Jac[i, j, k-1])
                local div = abs(Udxi + Vdelta + Wdsig)
                if (div > maxdiv)
                    maxdiv = div
                end
            end
        end
    end
    return maxdiv
end