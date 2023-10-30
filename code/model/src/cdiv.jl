function cdiv(dtimel, n)
    local maxdiv = 0e0
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                local Udxi = 0.5e0 * ((u[i+2, j+1, k+1, n+1] * ux[i+2, j+1] + v[i+2, j+1, k+1, n+1] * uy[i+2, j+1]) * Jac[i+2, j+1, k+1] 
                                    - (u[i, j+1, k+1, n+1] * ux[i, j+1] + v[i, j+1, k+1, n+1] * uy[i, j+1]) * Jac[i, j+1, k+1])
                local Vdelta = 0.5e0 * ((u[i+1, j+2, k+1, n+1] * vx[i+1, j+2] + v[i+1, j+2, k+1, n+1] * vy[i+1, j+2]) * Jac[i+1, j+2, k+1]
                                    - (u[i+1, j, k+1, n+1] * vx[i+1, j] + v[i+1, j, k+1, n+1] * vy[i+1, j]) * Jac[i+1, j, k+1])
                local Wdsig = 0.5e0 * ((u[i+1, j+1, k+2, n+1] * wx[i+1, j+1, k+2] + v[i+1, j+1, k+2, n+1] * wy[i+1, j+1, k+2] + EPS * w[i+1, j+1, k+2, n+1] * wz[i+1, j+1, k+1]) * Jac[i+1, j+1, k+2]
                                    - (u[i+1, j+1, k, n+1] * wx[i+1, j+1, k] + v[i+1, j+1, k, n+1] * wy[i+1, j+1, k] + EPS * w[i+1, j+1, k, n+1] * wz[i+1, j+1, k]) * Jac[i+1, j+1, k])
                local div = abs(Udxi + Vdelta + Wdsig)
                if (div > maxdiv)
                    maxdiv = div
                end
            end
        end
    end
    return maxdiv
end