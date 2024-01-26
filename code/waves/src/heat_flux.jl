
const global heat_flux_swrd = zeros(rc_kind, (NK))
const global heat_flux_Kdfluxdzz = zeros(rc_kind, (NK))

"""
function heat_flux(Tdif, step)
    # arguments
    - `Tdif` 
"""
function heat_flux(Tdif, step)
    local swrd = heat_flux_swrd
    local Kdfluxdzz = heat_flux_Kdfluxdzz

    # This part deals with the water type coefficients
    local J_lambda1 = 0.6e0
    @fortSet("j_lambda1", J_lambda1, rc_kind)
    local J_lambda2 = 20e0
    @fortSet("j_lambda2", J_lambda2, rc_kind)
    local J_A = 0.62e0
    @fortSet("j_a", J_A, rc_kind)

    # This part deals with short wave radiation
    @views @. swr = 0e0
    @views @. qloss = 0e0


    local fac = 1e0 / (UL * DL * delta)
    for j in 1:NJ
        for i in 1:NI
            local Kdfluxdzt = DL * (swr[j] - qloss[j]) / (R0 * 4187e0)
            for k in 1:NK
                swrd[k] = swr[j] * (J_A * exp(zf[div(NI, 2), div(NJ, 2), k] * DL / J_lambda1) + (1 - J_A) * exp(zf[div(NI, 2), div(NJ, 2), k] * DL / J_lambda2))
                Kdfluxdzz[k] = DL * swrd[k] / (R0 * 4187e0)
            end
            local Kdfluxdzb = 0e0
            Tdif[i, j, NK] = fac * Jac[i, j, NK] * wz[i, j, NK] * Kdfluxdzt
            for k in 2:NK-1
                Tdif[i, j, k] = fac * Jac[i, j, k] * wz[i, j, k] * (Kdfluxdzz[k] - Kdfluxdzz[k-1])
            end
            Tdif[i, j, 1] = fac * Jac[i, j, 1] * wz[i, j, 1] * Kdfluxdzb
        end
    end

    #@ccall "./PSOM_LIB.so".heat_flux_(pointer(Tdif)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid
end