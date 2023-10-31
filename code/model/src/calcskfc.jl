
using .fortVar

function calcfkfc()

    local be2 = @fortGet("beta",rc_kind) * EPS * EPS
    local kaph1 = 1e0 - @fortGet("kappah", rc_kind)

    for i in 1:NI
        for j in 1:NJ
            local hxi = 0.5e0 * (h[i+2, j+1] - h[i, j+1])
            local heta = 0.5e0 * (h[i+1, j+2] - h[i+1, j])
            
            local hx = ux[i+1, j+1] * hxi + vx[i+1, j+1] * heta
            local hy = uy[i+1, j+1] * hxi + vy[i+1, j+1] * heta
            
            hx = gpr * (@fortGet("kappah", rc_kind) * hx + kaph1 * gradhn[i+1, j+1, 1])
            hy = gpr * (@fortGet("kappah", rc_kind) * hy + kaph1 * gradhn[i+1, j+1, 2])
            

            for k in 1:NK-1
                local wxk = 0.5e0 * (wx[i+1, j+1, k+2] + wx[i+1, j+1, k+1])
                local wxsk = wxk * (0.5e0 * (si[i+1,j+1,k+2] + si[i+1,j+1,k+1]) + hx)

                local wyk = 0.5e0 * (wy[i+1, j+1, k+2] + wy[i+1, j+1, k+1])
                local wysk = wyk * (0.5e0 * (sj[i+1,j+1,k+2] + sj[i+1,j+1,k+1]) + hy)
                
                local wzsk = wzk[i+1, j+1, k+1] * (sk[i+1, j+1, k+2] + sk[i+1, j+1, k+1]) * 0.5e0
                
                local Jack = 0.5e0 * (Jac[i+1, j+1, k+2] + Jac[i+1, j+1, k+1])
                skfc[i+1, j+1, k+1] = (be2 * wzsk + wxsk + wysk) * Jack 
            end

            local k = 0

            local wzsk = wzk[i+1, j+1, k+1] * 0.5e0 * (3e0 * sk[i+1,j+1,k+2] - sk[i+1,j+1, k+3])
            local wxsk = wx[i+1, j+1, k+2] * (si[i+1, j+1, k+2] + hx)
            local wysk = wy[i+1, j+1, k+2] * (sj[i+1, j+1, k+2] + hy)

            local Jack = Jac[i+1, j+1, k+2]
            skfc[i+1, j+1, k+1] = (be2 * wzsk + wxsk + wysk) * Jack

            local k = NK

            local wzsk = wzk[i+1, j+1, k+1] * 0.5e0 * (3e0 * sk[i+1,j+1,k+1] - sk[i+1,j+1, k])
            local wxsk = 0.5e0 * (wx[i+1, j+1, k+2] + wx[i+1, j+1, k+1]) * (si[i+1,j+1,k+1] + hx)
            local wysk = 0.5e0 * (wy[i+1, j+1, k+2] + wy[i+1, j+1, k+1]) * (sj[i+1,j+1,k+1] + hy)

            local Jack = Jac[i+1, j+1, k+1]
            skfc[i+1, j+1, k+1] = (be2 * wzsk + wxsk + wysk) * Jack
            end
    end
end