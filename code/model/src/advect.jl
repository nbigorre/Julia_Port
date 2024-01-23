using .fortVar

function ultim(up::rc_kind, ct::rc_kind, dn::rc_kind, gdn::rc_kind, gup::rc_kind, courinv::rc_kind, fc::rc_kind)
    local del = dn - up
    local adel = abs(del)
    local acurv = abs(gdn - gup)
    if (acurv >= adel)
        return ct
    else
        local ref = up + (ct - up) * courinv
        if (del > 0)
            fc = max(ct, fc)
            ref = min(ref, dn)
            fc = min(fc, ref)
        else
            fc = min(ct, fc)
            ref = max(ref, dn)
            fc = max(fc, ref)
        end
    end
    return fc
end

function advect_x(area, var, uf, varflx)
    local Eighth::rc_kind = 0.125e0
    local Half::rc_kind = 0.5e0
    local courinv::rc_kind = 2e0
    local dvarx = OffsetArray(zeros(rc_kind, NI+2), 0:NI+1)
    for cindex in area
        local k,j = Tuple(cindex)
        dvarx[0] = var[1, j, k] - var[NI, j, k]

        @. @views dvarx[1:NI-1] = var[2:NI,j,k] - var[1:NI-1, j, k]

        dvarx[NI] = dvarx[0]
        dvarx[NI+1] = dvarx[1]

        var[NI+1, j, k] = var[1, j, k]
        var[0, j, k] = var[NI, j, k]

        for i in 1:NI
            if (uf[i, j, k] > 0)
                local left = Eighth * (dvarx[i] - dvarx[i-1])
                local ctr = Half * (var[i+1, j, k] + var[i, j, k])
                local fc = ctr - left
                fc = ultim(var[i-1, j, k], var[i, j, k], var[i+1, j, k], dvarx[i], dvarx[i-1], courinv, fc)
                varflx[i, j, k] = uf[i, j, k] * fc
            else
                local i2 = i == NI ? 0 : i
                local right = Eighth * (dvarx[i2+1] - dvarx[i2])
                local ctr = Half * (var[i2+1, j, k] + var[i, j, k])
                local fc = ctr - right
                fc = ultim(var[i2+2, j, k], var[i2+1, j, k], var[i, j, k], dvarx[i], dvarx[i2+1], courinv, fc)
                varflx[i, j, k] = uf[i, j, k] * fc
            end
        end
        varflx[0, j, k] = varflx[NI, j, k]
    end
end

function advect_y(area, var, vf, varflx)
    local Eighth = 0.125e0
    local Half = 0.5e0
    local courinv = 2e0
    local dvary = OffsetArray(zeros(rc_kind, NJ+1), 0:NJ)
    for cindex in area
        local k, i = Tuple(cindex)
        dvary[0] = 0e0
        dvary[NJ] = 0e0

        @. @views dvary[1:NJ-1] = var[i, 2:NJ, k] - var[i, 1:NJ-1, k]
        
        var[i, 0, k] = var[i, 1, k]
        var[i, NJ+1, k] = var[i, NJ, k]

        for j in 1:NJ-1
            if (vf[i, j, k] >= 0)
                local left = Eighth * (dvary[j] - dvary[j-1])
                local ctr = Half * (var[i, j+1, k] + var[i, j, k])
                local fc = ctr - left
                fc = ultim(var[i, j-1, k], var[i, j, k], var[i, j+1, k], dvary[j], dvary[j-1], courinv, fc)
                varflx[i, j, k] = vf[i, j, k] * fc
            else
                local right = Eighth * (dvary[j+1] - dvary[j])
                local ctr = Half * (var[i, j+1, k] + var[i, j, k])
                local fc = ctr - right
                fc = ultim(var[i, j+2, k], var[i, j+1, k], var[i, j, k], dvary[j], dvary[j+1], courinv, fc)
                varflx[i, j, k] = vf[i, j, k] * fc
            end
        end
        varflx[i, 0, k] = 0e0
        varflx[i, NJ, k] = 0e0
    end
end

function advect_z(area, var, wf, varflx)
    local Eighth = 0.125e0
    local Half = 0.5e0
    local courinv = 2e0
    local dvarz = OffsetArray(zeros(rc_kind, NK+2), 0:NK+1)
    for cindex in area
        local j, i = Tuple(cindex)

        var[i, j, 0] = var[i, j, 1]

        dvarz[0] = 0e0

        @. @views dvarz[1:NK-1] = var[i, j, 2:NK] - var[i, j, 1:NK-1]
        
        dvarz[NK:NK+1] .= 0e0

        var[i, j, NK+1] = var[i, j, NK]

        for k in 1:NK
            if (wf[i, j, k] >= 0)
                local left = Eighth * (dvarz[k] - dvarz[k-1])
                local ctr = Half * (var[i, j, k+1] + var[i, j, k])
                local fc = ctr - left
                fc = ultim(var[i, j, k-1], var[i, j, k], var[i, j, k+1], dvarz[k], dvarz[k-1], courinv, fc)
                varflx[i, j, k] = wf[i, j, k] * fc
            else
                if (k == NK)
                    varflx[i, j, k] = wf[i, j, k] * Half * (var[i, j, k+1] + var[i, j, k])
                else
                    local right = Eighth * (dvarz[k+1] - dvarz[k])
                    local ctr = Half * (var[i, j, k+1] + var[i, j, k])
                    local fc = ctr - right
                    fc = ultim(var[i, j, k+2], var[i, j, k+1], var[i, j, k], dvarz[k], dvarz[k+1], courinv, fc)
                    varflx[i, j, k] = wf[i, j, k] * fc
                end
            end
        end
        varflx[i, j, 0] = 0e0
    end
end



function advect(var, uvarx)
    local varflx = OffsetArray(zeros(rc_kind, (NI + 1, NJ + 1, NK + 1)), 0:NI, 0:NJ, 0:NK)

    local chunks = Iterators.partition(CartesianIndices((1:NK, 1:NJ)), div(NJ * NK, Threads.nthreads()))

    local tasks = map(chunks) do area
        Threads.@spawn advect_x(area, var, uf, varflx)
    end
    wait.(tasks)

    
    @. @views uvarx[1:NI, 1:NJ, 1:NK] = varflx[1:NI, 1:NJ, 1:NK] - varflx[0:NI-1, 1:NJ, 1:NK]
    

    local chunks = Iterators.partition(CartesianIndices((1:NK, 1:NI)), div(NI * NK, Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn advect_y(area, var, vf, varflx)
    end
    wait.(tasks)

    @. @views uvarx[1:NI, 1:NJ, 1:NK] += (varflx[1:NI, 1:NJ, 1:NK] - varflx[1:NI, 0:NJ-1, 1:NK])




    local chunks = Iterators.partition(CartesianIndices((1:NJ, 1:NI)), div(NI * NJ, Threads.nthreads()))
    local tasks = map(chunks) do area
        Threads.@spawn advect_z(area, var, wf, varflx)
    end
    wait.(tasks)

    @. @views uvarx[1:NI, 1:NJ, 1:NK] += (varflx[1:NI, 1:NJ, 1:NK] - varflx[1:NI, 1:NJ, 0:NK-1])



end