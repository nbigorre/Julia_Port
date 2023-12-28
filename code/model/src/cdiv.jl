function cdiv(dtimel, n)
    local maxdiv = 0e0
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                local Udxi = 0.5e0 * ((u[i+2, j+1, k+1, n+1] * ux[i+1, j] + v[i+2, j+1, k+1, n+1] * uy[i+1, j]) * Jac[i+1, j, k] 
                                    - (u[i, j+1, k+1, n+1] * ux[i-1, j] + v[i, j+1, k+1, n+1] * uy[i-1, j]) * Jac[i-1, j, k])
                local Vdelta = 0.5e0 * ((u[i+1, j+2, k+1, n+1] * vx[i, j+1] + v[i+1, j+2, k+1, n+1] * vy[i, j+1]) * Jac[i, j+1, k]
                                    - (u[i+1, j, k+1, n+1] * vx[i, j-1] + v[i+1, j, k+1, n+1] * vy[i, j-1]) * Jac[i, j-1, k])
                local Wdsig = 0.5e0 * ((u[i+1, j+1, k+2, n+1] * wx[i, j, k+1] + v[i+1, j+1, k+2, n+1] * wy[i, j, k+1] + EPS * w[i+1, j+1, k+2, n+1] * wz[i, j, k]) * Jac[i, j, k+1]
                                    - (u[i+1, j+1, k, n+1] * wx[i, j, k-1] + v[i+1, j+1, k, n+1] * wy[i, j, k-1] + EPS * w[i+1, j+1, k, n+1] * wz[i, j, k-1]) * Jac[i, j, k-1])
                local div = abs(Udxi + Vdelta + Wdsig)
                if (div > maxdiv)
                    maxdiv = div
                end
            end
        end
    end
    return maxdiv
end