using .header
using .fortVar
using .cppdefs


function mixing_vertical_worker(area, fac, var, vardif)
    local dvardzfc = zeros(rc_kind, NK + 1)
    for counter in area
        local j = div(counter - 1, NI) + 1
        local i = mod(counter - 1, NI) + 1
        for k in 1:NK-1
            dvardzfc[k+1] = wz[i+1, j+1, k+1] * (var[i+1, j+1, k+2] - var[i+1, j+1, k+1])
        end

        dvardzfc[1] = 0e0
        dvardzfc[NK+1] = 0e0

        local dfactor = fac * Jac[i, j, 1] * wz[i+1, j+1, 2]
        @static if (cppdefs.implicit)
            #########################################################################################################
        else
            vardif[i, j, 1] = dfactor * (dvardzfc[2] * Kz[i, j, 1] - dvardzfc[1] * Kz[i, j, 0])
        end
        for k in 2:NK-1
            local dfactor = fac * Jac[i, j, k] * wz[i+1, j+1, k+1]
            @static if (cppdefs.implicit)
                #########################################################################################################
            else
                vardif[i, j, k] = dfactor * (dvardzfc[k+1] * Kz[i, j, k] - dvardzfc[k] * Kz[i, j, k-1])
            end
        end
        local dfactor = fac * Jac[i, j, NK] * wz[i+1, j+1, NK+1]
        @static if (cppdefs.implicit)
            #########################################################################################################
        else
            vardif[i, j, NK] = dfactor * (dvardzfc[NK+1] * Kz[i, j, NK] - dvardzfc[NK] * Kz[i, j, NK-1])
        end

    end
end

function mixing_vertical(var, vardif, m, step, iv_compute_kzl)
    local Kzmax = 1e-3
    local KzmaxTr = 1e-3
    local facb = @fortGet("rr", rc_kind) * DL
    local fac = 1e0 / (@fortGet("ul", rc_kind) * DL * @fortGet("delta", rc_kind))
    local fact = DL / @fortGet("ul", rc_kind)
    @static if (cppdefs.gotm_call)
        ###################################################################################################################################
    else
        if (iv_compute_kzl == 1)
            for j in 1:NJ
                for i in 1:NI
                    local ustar = sqrt(stress_top[i, j] / R0)
                    local ff = ffc[i, j] * FPAR
                    local Ekdepth = 0.4e0 * ustar / ff
                    local zextent = 0.5e0 * Ekdepth
                    local ztransit = -Ekdepth
                    for k in NK:-1:0
                        Kz[i, j, k] = 1e0
                        local zdep = zf[i+1, j+1, k+2] * DL
                        local thy = (1e0 + tanh(((zdep - ztransit) / zextent) * PI)) * 0.5e0
                        Kz[i, j, k] = max(0.01e0, thy) * KzmaxTr
                    end
                end
            end
        end
    end
    
    local chunks = Iterators.partition(1:NJ*NI, div(NJ*NI, Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn mixing_vertical_worker(area, fac, var, vardif)
    end
    wait.(tasks)
    
end

function viscosity(dudz, dvdz, drdz, i, j)
    local n1 = 3
    local grho = 9.81 / R0
    local DLinv = @fortGet("dlinv", rc_kind)
    local fac = (UL^2) / (@fortGet("dl", rc_kind)^2)
    local RiCr = 0.7e0
    Kz[i, j, 0] = 0e0
    Kz[i, j, NK] = 0e0

    for k in 1:NK-1
        local bvfreq = -grho * drdz[k+1] * DLinv
        local vshear = ((dudz[k+1]^2) + (dvdz[k+1]^2)) * fac
        local Ri = vshear == 0 ? 100e0 : bvfreq / vshear
        if (Ri == 0)
            Kz[i, j, k] = 1e0
        elseif (Ri == RiCr)
            Kz[i, j, k] = (1e0 - (Ri^2) / (RiCr^2))^n1
        else
            Kz[i, j, k] = 1e0
        end
    end

end