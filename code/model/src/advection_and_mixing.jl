using .header
using .fortVar
using .cppdefs

function advection_and_mixing(m, n, dtimel, step)

    local iv_compute_kz = 1
    local av_comp = ones(rc_kind, 5 + ntr)

    local var = zeros(rc_kind, (NI + 2, NJ + 2, NK + 2))
    local var0 = zeros(rc_kind, (NI + 2, NJ + 2, NK + 2))
    
    if (cppdefs.rhoonly)
        av_comp[1] = 0
    end
    local u_w = zeros(rc_kind, (NI, NJ, NK))
    local v_w = zeros(rc_kind, (NI, NJ, NK))
    local T_h = zeros(rc_kind, (NI, NJ, NK))
    
    heat_flux(T_h, step)
    #@ccall "./PSOM_LIB.so".heat_flux_(pointer(T_h)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid
    
    @ccall "./PSOM_LIB.so".wind_stress_(pointer(u_w)::Ptr{rc_kind}, pointer(v_w)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid
    
    local uvarx = zeros(rc_kind, (NI + 2, NJ + 2, NK + 2))
    local uvarx_advec = zeros(rc_kind, (NI + 2, NJ + 2, NK + 2))
    local vardif = zeros(rc_kind, (NI, NJ, NK))


    for selvar in 1:5+ntr
        @fortSet("selvar", selvar, rc_kind)
        if (cppdefs.implicit)
            @inbounds @views @. mat_A .= 0e0
            @inbounds @views @. mat_B .= 0e0
            @inbounds @views @. mat_C .= 0e0
            @inbounds @views @. mat_D .= 0e0
            @inbounds @views @. mat_E .= 0e0
        end

        @inbounds @views @. uvarx .= 0.0
        @inbounds @views @. uvarx_advec .= 0.0
        @inbounds @views @. vardif .= 0.0
        

        if (av_comp[selvar] == 1)
            if (selvar == 1)
                @inbounds @views @. var[:, :, :] = T[:, :, :, m+1]
                @inbounds @views @. var0[:, :, :] = T[:, :, :, 1]
            elseif (selvar == 2)
                @inbounds @views @. var[:, :, :] = s[:, :, :, m+1]
                @inbounds @views @. var0[:, :, :] = s[:, :, :, 1]
            elseif (selvar == 3)
                @inbounds @views @. var[:, :, :] = u[:, :, :, m+1]
                @inbounds @views @. var0[:, :, :] = u[:, :, :, 1]
            elseif (selvar == 4)
                @inbounds @views @. var[:, :, :] = v[:, :, :, m+1]
                @inbounds @views @. var0[:, :, :] = v[:, :, :, 1]
            elseif (selvar == 5)
                @inbounds @views @. var[:, :, :] = w[:, :, :, m+1]
                @inbounds @views @. var0[:, :, :] = w[:, :, :, 1]
            end
            for it in 1:ntr
                if (selvar == 5 + it)
                    @inbounds @views @. var[:, :, :] = Tr[it, :, :, :, m+1]
                end
            end

            advect(var, uvarx)
            #@ccall "./PSOM_LIB.so".advect_(pointer(var)::Ptr{rc_kind}, pointer(uvarx)::Ptr{rc_kind})::Cvoid

            if (selvar < 5.5)
                vardif .= 0.0
                mixing_horizontal(var, vardif)
                #@ccall "./PSOM_LIB.so".mixing_horizontal_(pointer(var)::Ptr{rc_kind}, pointer(vardif)::Ptr{rc_kind})::Cvoid
                @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2:NK+1] -= vardif[1:NI, 1:NJ, 1:NK]
            end

            @inbounds @views @. uvarx_advec[2:NI+1, 2:NJ+1, 2:NK+1] = uvarx[2:NI+1, 2:NJ+1, 2:NK+1]

            if (selvar < 4.5)
                vardif .= 0.0
                mixing_vertical(var, vardif, m, step, iv_compute_kz)
                #@ccall "./PSOM_LIB.so".mixing_vertical_(pointer(var)::Ptr{rc_kind}, pointer(vardif)::Ptr{rc_kind}, Ref(m)::Ptr{Int}, Ref(step)::Ptr{Int}, Ref(iv_compute_kz)::Ptr{Int})::Cvoid
                @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2:NK+1] -= vardif[1:NI, 1:NJ, 1:NK]
                iv_compute_kz = 0
            end

            if (selvar == 1)
                if (cppdefs.implicit)
                    mat_D[:,:,:,selvar] = T_h[1:NI, 1:NJ, 1:NK]
                else
                    @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2:NK+1] -= T_h[1:NI, 1:NJ, 1:NK]
                end
            end
            if (selvar == 3)
                if (cppdefs.implicit)
                    mat_D[:,:,:,selvar] = u_w[1:NI, 1:NJ, 1:NK]
                else
                    @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2:NK+1] -= u_w[1:NI, 1:NJ, 1:NK]
                end
            end
            if (selvar == 4)
                if (cppdefs.implicit)
                    mat_D[:,:,:,selvar] = v_w[1:NI, 1:NJ, 1:NK]
                else
                    @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2:NK+1] -= v_w[1:NI, 1:NJ, 1:NK]
                end
            end

            if (!cppdefs.implicit)
                if (selvar == 3 || selvar == 4)
                    @inbounds @views @. uvarx[2:NI+1, 2:NJ+1, 2] .+= @fortGet("rr", rc_kind) * (1e0 / (@fortGet("ul", rc_kind) * @fortGet("delta", rc_kind))) * (Jac[2:NI+1, 2:NJ+1, 2] * wz[2:NI+1, 2:NJ+1, 2]) * var[2:NI+1, 2:NJ+1, 2]
                end
            end
            if (cppdefs.implicit)
                if (selvar < 4.5)
                    for i in 1:NI
                        for j in 1:NJ
                            mat_A[i,j,1:NK, selvar] .*= dtimel * JacInv[i+1, j+1, 2:NK+1]
                            mat_B[i,j,1:NK, selvar] = (dtimel * JacInv[i+1, j+1, 2:NK+1] * mat_B[i,j,1:NK, selvar]) .- 1
                            mat_C[i,j,1:NK-1, selvar] .*= dtimel * JacInv[i+1, j+1, 2:NK]
                            mat_D[i,j,1:NK, selvar] .= 0e0 .- var0[i+1, j+1, 2:NK+1] - (dtimel * JacInv[i+1, j+1, 2:NK+1] * mat_D[i,j,1:NK, selvar]) + (dtimel * JacInv[i+1, j+1, 2:NK+1] * uvarx_advec[i+1,j+1,2:NK+1]) 
                            @ccall "./PSOM_LIB.so".solve_tridiag_(pointer(mat_A[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_B[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_C[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_D[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_test[i,j,1:NK, selvar])::Ptr{rc_kind}, Ref(NK)::Ref{rc_kind})::Cvoid
                        end
                    end
                    if (selvar == 1)
                        T[2:NI+1, 2:NJ+1, 2:NK+1, n+1] = mat_test[:,:,:,selvar]
                    end
                    if (selvar == 2)
                        s[2:NI+1, 2:NJ+1, 2:NK+1, n+1] = mat_test[:,:,:,selvar]
                    end
                    if (selvar == 3)
                        cx[2:NI+1, 2:NJ+1, 2:NK+1] = mat_test[:,:,:,selvar]
                    end
                    if (selvar == 4)
                        cy[2:NI+1, 2:NJ+1, 2:NK+1] = mat_test[:,:,:,selvar]
                    end

                end
            else
                local lc_diff = @inbounds @views @. dtimel * JacInv[2:NI+1, 2:NJ+1, 2:NK+1] .* uvarx[2:NI+1, 2:NJ+1, 2:NK+1]
                if (selvar == 1)
                    @inbounds @views @. T[2:NI+1, 2:NJ+1, 2:NK+1, n+1] = T[2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                end
                if (selvar == 2)
                    @inbounds @views @. s[2:NI+1, 2:NJ+1, 2:NK+1, n+1] = s[2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                end
                if (selvar == 3)
                    @inbounds @views @. cx[2:NI+1, 2:NJ+1, 2:NK+1] = u[2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                end
                if (selvar == 4)
                    @inbounds @views @. cy[2:NI+1, 2:NJ+1, 2:NK+1] = v[2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                end
            end
            if (selvar == 5)
                local lc_diff = @inbounds @views @. dtimel * JacInv[2:NI+1, 2:NJ+1, 2:NK+1] .* uvarx[2:NI+1, 2:NJ+1, 2:NK+1]
                
                @inbounds @views @. cz[2:NI+1, 2:NJ+1, 2:NK+1] = w[2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                for it in 1:ntr
                    if (selvar == 5 + it)
                             @views @. Tr[it, 2:NI+1, 2:NJ+1, 2:NK+1, n+1] = Tr[it, 2:NI+1, 2:NJ+1, 2:NK+1, 1] - lc_diff
                    end
                end
            end
        end
    end

end