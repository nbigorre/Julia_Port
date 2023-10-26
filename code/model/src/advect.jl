using .header
using .fortVar

function ultim(up::rc_kind, ct::rc_kind, dn::rc_kind, gdn::rc_kind, gup::rc_kind, courinv::rc_kind, fc::Ref{rc_kind})
    local del = dn - up
    local adel = abs(del)
    local acurv = abs(gdn - gup)
    if (acurv >= adel)
        fc[] = ct
        return
    else
        local ref = up + (ct - up) * courinv
        if (del > 0)
            fc[] = max(ct, fc[])
            ref = min(ref, dn)
            fc[] = min(fc[], ref)
        else
            fc[] = min(ct, fc[])
            ref = max(ref, dn)
            fc[] = max(fc[], ref)
        end
    end
    return fc[]
end

function advect_x(area, var, uf, varflx)
    local Eighth::rc_kind = 0.125e0
    local Half::rc_kind = 0.5e0
    local courinv::rc_kind = 2e0
    local dvarx::Array{rc_kind} = zeros(rc_kind, NI + 2)
    for cindex in area
        local k,j = Tuple(cindex)
        dvarx[1] = var[2, j+1, k+1] - var[NI+1, j+1, k+1]

        @inbounds @. @views dvarx[2:NI] = var[3:NI+1,j+1,k+1] - var[2:NI, j+1, k+1]

        dvarx[NI+1] = dvarx[1]
        dvarx[NI+2] = dvarx[2]

        var[NI+2, j+1, k+1] = var[2, j+1, k+1]
        var[1, j+1, k+1] = var[NI+1, j+1, k+1]

        for i in 1:NI
            if (uf[i+1, j, k] > 0)
                local left = Eighth * (dvarx[i+1] - dvarx[i])
                local ctr = Half * (var[i+2, j+1, k+1] + var[i+1, j+1, k+1])
                local fc = ctr - left
                let a7 = Ref(fc)
                    ultim(var[i, j+1, k+1], var[i+1, j+1, k+1], var[i+2, j+1, k+1], dvarx[i+1], dvarx[i], courinv, a7)
                    fc = a7[]
                end
                varflx[i+1, j+1, k+1] = uf[i+1, j, k] * fc
            else
                local i2 = i == NI ? 0 : i
                local right = Eighth * (dvarx[i2+2] - dvarx[i2+1])
                local ctr = Half * (var[i2+2, j+1, k+1] + var[i+1, j+1, k+1])
                local fc = ctr - right
                let a7 = Ref(fc)
                    ultim(var[i2+3, j+1, k+1], var[i2+2, j+1, k+1], var[i+1, j+1, k+1], dvarx[i+1], dvarx[i2+2], courinv, a7)
                    fc = a7[]
                end
                varflx[i+1, j+1, k+1] = uf[i+1, j, k] * fc
            end
        end
        varflx[1, j+1, k+1] = varflx[NI+1, j+1, k+1]
    end
end

function advect_y(area, var, vf, varflx)
    local Eighth = 0.125e0
    local Half = 0.5e0
    local courinv = 2e0
    local dvary = zeros(rc_kind, NJ + 1)
    for cindex in area
        local k, i = Tuple(cindex)
        dvary[1] = 0e0
        dvary[NJ+1] = 0e0

        @inbounds @. @views dvary[2:NJ] = var[i+1, 3:NJ+1, k+1] - var[i+1, 2:NJ, k+1]
        
        var[i+1, 1, k+1] = var[i+1, 2, k+1]
        var[i+1, NJ+2, k+1] = var[i+1, NJ+1, k+1]

        for j in 1:NJ-1
            if (vf[i, j+1, k] >= 0)
                local left = Eighth * (dvary[j+1] - dvary[j])
                local ctr = Half * (var[i+1, j+2, k+1] + var[i+1, j+1, k+1])
                local fc = ctr - left
                let a7 = Ref(fc)
                    ultim(var[i+1, j, k+1], var[i+1, j+1, k+1], var[i+1, j+2, k+1], dvary[j+1], dvary[j], courinv, a7)
                    fc = a7[]
                end
                varflx[i+1, j+1, k+1] = vf[i, j+1, k] * fc
            else
                local right = Eighth * (dvary[j+2] - dvary[j+1])
                local ctr = Half * (var[i+1, j+2, k+1] + var[i+1, j+1, k+1])
                local fc = ctr - right
                let a7 = Ref(fc)
                    ultim(var[i+1, j+3, k+1], var[i+1, j+2, k+1], var[i+1, j+1, k+1], dvary[j+1], dvary[j+2], courinv, a7)
                    fc = a7[]
                end
                varflx[i+1, j+1, k+1] = vf[i, j+1, k] * fc
            end
        end
        varflx[i+1, 1, k+1] = 0e0
        varflx[i+1, NJ+1, k+1] = 0e0
    end
end

function advect_z(area, var, wf, varflx)
    local Eighth = 0.125e0
    local Half = 0.5e0
    local courinv = 2e0
    local dvarz = zeros(rc_kind, NK + 2)
    for cindex in area
        local j, i = Tuple(cindex)

        var[i+1, j+1, 1] = var[i+1, j+1, 2]

        dvarz[1] = 0e0

        @inbounds @. @views dvarz[2:NK] = var[i+1, j+1, 3:NK+1] - var[i+1, j+1, 2:NK]
        
        dvarz[NK+1:NK+2] .= 0e0

        var[i+1, j+1, NK+2] = var[i+1, j+1, NK+1]

        for k in 1:NK
            if (wf[i, j, k+1] >= 0)
                local left = Eighth * (dvarz[k+1] - dvarz[k])
                local ctr = Half * (var[i+1, j+1, k+2] + var[i+1, j+1, k+1])
                local fc = ctr - left
                let a7 = Ref(fc)
                    ultim(var[i+1, j+1, k], var[i+1, j+1, k+1], var[i+1, j+1, k+2], dvarz[k+1], dvarz[k], courinv, a7)
                    fc = a7[]
                end
                varflx[i+1, j+1, k+1] = wf[i, j, k+1] * fc
            else
                if (k == NK)
                    varflx[i+1, j+1, k+1] = wf[i, j, k+1] * Half * (var[i+1, j+1, k+2] + var[i+1, j+1, k+1])
                else
                    local right = Eighth * (dvarz[k+2] - dvarz[k+1])
                    local ctr = Half * (var[i+1, j+1, k+2] + var[i+1, j+1, k+1])
                    local fc = ctr - right
                    let a7 = Ref(fc)
                        ultim(var[i+1, j+1, k+3], var[i+1, j+1, k+2], var[i+1, j+1, k+1], dvarz[k+1], dvarz[k+2], courinv, a7)
                        fc = a7[]
                    end
                    varflx[i+1, j+1, k+1] = wf[i, j, k+1] * fc
                end
            end
        end
        varflx[i+1, j+1, 1] = 0e0
    end
end



function advect(var::Array{rc_kind}, uvarx::Array{rc_kind})
    local varflx = zeros(rc_kind, (NI + 1, NJ + 1, NK + 1))

    local chunks = Iterators.partition(CartesianIndices((1:NK, 1:NJ)), div(NJ * NK, Threads.nthreads()))

    local tasks = map(chunks) do area
        Threads.@spawn advect_x(area, var, uf, varflx)
    end
    wait.(tasks)

    
    @inbounds @. @views uvarx[2:NI+1, 2:NJ+1, 2:NK+1] = varflx[2:NI+1, 2:NJ+1, 2:NK+1] - varflx[1:NI, 2:NJ+1, 2:NK+1]
    

    local chunks = Iterators.partition(CartesianIndices((1:NK, 1:NI)), div(NI * NK, Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn advect_y(area, var, vf, varflx)
    end
    wait.(tasks)

    @inbounds @. @views uvarx[2:NI+1, 2:NJ+1, 2:NK+1] += (varflx[2:NI+1, 2:NJ+1, 2:NK+1] - varflx[2:NI+1, 1:NJ, 2:NK+1])




    local chunks = Iterators.partition(CartesianIndices((1:NJ, 1:NI)), div(NI * NJ, Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn advect_z(area, var, wf, varflx)
    end
    wait.(tasks)

    @inbounds @. @views uvarx[2:NI+1, 2:NJ+1, 2:NK+1] += (varflx[2:NI+1, 2:NJ+1, 2:NK+1] - varflx[2:NI+1, 2:NJ+1, 1:NK])



end