using .header
using .fortVar

function mixing_horizontal_l1(area, vy, var, Ky_l, vardif)
    local dvardyfc = zeros(rc_kind, NJ + 1)
    for cindex in area
        local i,k = Tuple(cindex)
        for j in 1:NJ-1
            dvardyfc[j+1] = 0.5 * (vy[i+1, j+1] + vy[i+1, j+2]) * (var[i+1, j+2, k+1] - var[i+1, j+1, k+1])
        end
        dvardyfc[1] = 0e0
        dvardyfc[NJ+1] = 0e0
        for j in 1:NJ
            vardif[i, j, k] = Ky_l[j] * vy[i+1, j+1] * (dvardyfc[j+1] - dvardyfc[j])
        end

    end
end

function mixing_horizontal_l2(area, ux, var, Kx_l, vardif)
    local dvardxfc = zeros(rc_kind, NI + 1)
    for cindex in area
        local j,k = Tuple(cindex)

        for i in 1:NI-1
            dvardxfc[i+1] = 0.5 * (ux[i+1, j+1] + ux[i+2, j+1]) * (var[i+2, j+1, k+1] - var[i+1, j+1, k+1])
        end
        dvardxfc[1] = 0.5 * (ux[NI+1, j+1] + ux[2, j+1]) * (var[2, j+1, k+1] - var[NI+1, j+1, k+1])
        dvardxfc[NI+1] = dvardxfc[1]
        for i in 1:NI
            vardif[i, j, k] += Kx_l * ux[i+1, j+1] * (dvardxfc[i+1] - dvardxfc[i])
        end

    end
end

function mixing_horizontal(var, vardif)
    local Kx_l = @fortGet("kx", rc_kind)
    local Ky_l = fill(@fortGet("ky", rc_kind), NJ)
    
    local chunks = Iterators.partition(CartesianIndices((1:NI, 1:NK)), div(NI*NK,Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn mixing_horizontal_l1(area, vy, var, Ky_l, vardif)
    end
    wait.(tasks)

    local chunks = Iterators.partition(CartesianIndices((1:NJ, 1:NK)), div(NJ*NK,Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn mixing_horizontal_l2(area, ux, var, Kx_l, vardif)
    end
    wait.(tasks)

    local fac = 1e0 / (@fortGet("ul", rc_kind) * LEN)
    @inbounds @views @. vardif[1:NI, 1:NJ, 1:NK] .*= fac * Jac[2:NI+1, 2:NJ+1, 2:NK+1]

end