
function wind_stress(udif,vdif,step)
    
    @inbounds @views @. udif = 0e0
    @inbounds @views @. vdif = 0e0

    @inbounds @views @. stress_top_x = 0e0
    @inbounds @views @. stress_top_y = 0e0

    local stressmax = 0.2e0

    local iyw = 120
    local iy0 = 40
    
    local yw = rc_kind(iyw) * @fortGet("dy", rc_kind)
    local y0 = rc_kind(iy0) * @fortGet("dy", rc_kind) - @fortGet("dy", rc_kind) / 2e0
    local yset = y0 + yw / 2e0
    
    for j in 1:div(NJ,2)
        stress_top_y[1,j] = stressmax * 0.5 * (1e0 + tanh(((yc[j] * 1e3 - yset) / yw) * PI)) 
    end
    local yset = yc[NJ] * 1e3 - y0 - yw / 2e0
    for j in div(NJ,2)+1:NJ
        stress_top_y[1,j] = stressmax * 0.5 * (1e0 + tanh((-(yc[j] * 1e3 - yset) / yw) * PI)) 
    end
    
    for i in 1:NI
        for j in 1:NJ
            stress_top_y[i,j] = stress_top_y[1,j]
        end
    end

    
    local fac = 1e0/(@fortGet("ul",rc_kind) * DL * @fortGet("delta",rc_kind))
    local fact = DL/@fortGet("ul",rc_kind)
    for j in 1:NJ
        for i in 1:NI
            stress_top[i,j] = sqrt(stress_top_x[i,j]*stress_top_x[i,j] + stress_top_y[i,j]*stress_top_y[i,j])
            local rhoinv = 1e0 / rho[i,j,NK]
            
            local Kdudzt = stress_top_x[i,j] * rhoinv * fact
            local Kdvdzt = stress_top_y[i,j] * rhoinv * fact
            
            udif[i,j,NK] = fac * Jac[i,j,NK] * wz[i,j,NK] * Kdudzt
            vdif[i,j,NK] = fac * Jac[i,j,NK] * wz[i,j,NK] * Kdvdzt
        end
    end


end