using .header
using .fortVar

function mixing_horizontal_l1(area, vy, var, Ky_l, vardif)
    local dvardyfc = OffsetArray(zeros(rc_kind, NJ + 1), 0:NJ)
    for cindex in area
        local i,k = Tuple(cindex)
        for j in 1:NJ-1
            dvardyfc[j] = 0.5 * (vy[i, j] + vy[i, j+1]) * (var[i, j+1, k] - var[i, j, k])
        end
        dvardyfc[0] = 0e0
        dvardyfc[NJ] = 0e0
        for j in 1:NJ
            vardif[i, j, k] = Ky_l[j] * vy[i, j] * (dvardyfc[j] - dvardyfc[j-1])
        end

    end
end

function mixing_horizontal_l2(area, ux, var, Kx_l, vardif)
    local dvardxfc = OffsetArray(zeros(rc_kind, NI + 1), 0:NI)
    for cindex in area
        local j,k = Tuple(cindex)

        for i in 1:NI-1
            dvardxfc[i] = 0.5 * (ux[i, j] + ux[i+1, j]) * (var[i+1, j, k] - var[i, j, k])
        end
        dvardxfc[0] = 0.5 * (ux[NI, j] + ux[1, j]) * (var[1, j, k] - var[NI, j, k])
        dvardxfc[NI] = dvardxfc[0]
        for i in 1:NI
            vardif[i, j, k] += Kx_l * ux[i, j] * (dvardxfc[i] - dvardxfc[i-1])
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
    @views @. vardif[1:NI, 1:NJ, 1:NK] .*= fac * Jac[1:NI, 1:NJ, 1:NK]

end