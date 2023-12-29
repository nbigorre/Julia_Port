
using .fortVar

function calcfkfc()

    local be2 = @fortGet("beta",rc_kind) * EPS * EPS
    local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

    for i in 1:NI
        for j in 1:NJ
            local hxi = 0.5e0 * (h[i+1, j] - h[i-1, j])
            local heta = 0.5e0 * (h[i, j+1] - h[i, j-1])
            
            local hx = ux[i+1, j+1] * hxi + vx[i, j] * heta
            local hy = uy[i, j] * hxi + vy[i, j] * heta
            
            hx = gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * gradhn[i, j, 1])
            hy = gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * gradhn[i, j, 2])
            

            for k in 1:NK-1
                local wxk = 0.5e0 * (wx[i, j, k+1] + wx[i, j, k])
                local wxsk = wxk * (0.5e0 * (si[i,j,k+1] + si[i,j,k]) + hx)

                local wyk = 0.5e0 * (wy[i, j, k+1] + wy[i, j, k])
                local wysk = wyk * (0.5e0 * (sj[i,j,k+1] + sj[i,j,k]) + hy)
                
                local wzsk = wzk[i, j, k] * (sk[i, j, k+1] + sk[i, j, k]) * 0.5e0
                
                local Jack = 0.5e0 * (Jac[i, j, k+1] + Jac[i, j, k])
                skfc[i, j, k] = (be2 * wzsk + wxsk + wysk) * Jack 
            end

            local k = 0

            local wzsk = wzk[i, j, k] * 0.5e0 * (3e0 * sk[i,j,k+1] - sk[i,j, k+2])
            local wxsk = wx[i, j, k+1] * (si[i, j, k+1] + hx)
            local wysk = wy[i, j, k+1] * (sj[i, j, k+1] + hy)

            local Jack = Jac[i, j, k+1]
            skfc[i, j, k] = (be2 * wzsk + wxsk + wysk) * Jack

            local k = NK

            local wzsk = wzk[i, j, k] * 0.5e0 * (3e0 * sk[i,j,k] - sk[i,j, k-1])
            local wxsk = 0.5e0 * (wx[i, j, k+1] + wx[i, j, k]) * (si[i,j,k] + hx)
            local wysk = 0.5e0 * (wy[i, j, k+1] + wy[i, j, k]) * (sj[i,j,k] + hy)

            local Jack = Jac[i, j, k]
            skfc[i, j, k] = (be2 * wzsk + wxsk + wysk) * Jack
            end
    end
end