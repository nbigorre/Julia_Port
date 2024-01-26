using .fortVar
using .cppdefs

function advection_and_mixing(m, n, dtimel, step)

    local iv_compute_kz = 1
    local av_comp = ones(rc_kind, 5 + ntr)

    local var = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
    local var0 = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)

    if (cppdefs.rhoonly)
        av_comp[1] = 0
    end
    @views @. mat_A .= 0e0
    @views @. mat_B .= 0e0
    @views @. mat_C .= 0e0
    @views @. mat_D .= 0e0
    @views @. mat_E .= 0e0
    local u_w = zeros(rc_kind, (NI, NJ, NK))
    local v_w = zeros(rc_kind, (NI, NJ, NK))

    wind_stress(u_w, v_w, step)
    #@ccall "./PSOM_LIB.so".wind_stress_(pointer(u_w)::Ptr{rc_kind}, pointer(v_w)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid

    local uvarx = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
    local uvarx_advec = OffsetArray(zeros(rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
    local vardif = zeros(rc_kind, (NI, NJ, NK))


    for selvar in 1:5+ntr
        @fortSet("selvar", selvar, rc_kind)



        if (av_comp[selvar] == 1)
            if (selvar == 1)
                @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = T[0:NI+1, 0:NJ+1, 0:NK+1, m]
            elseif (selvar == 2)
                @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = s[0:NI+1, 0:NJ+1, 0:NK+1, m]
            elseif (selvar == 3)
                @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = u[0:NI+1, 0:NJ+1, 0:NK+1, m]
            elseif (selvar == 4)
                @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = v[0:NI+1, 0:NJ+1, 0:NK+1, m]
            elseif (selvar == 5)
                @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = w[0:NI+1, 0:NJ+1, 0:NK+1, m]
            end
            for it in 1:ntr
                if (selvar == 5 + it)
                    @views @. var[0:NI+1, 0:NJ+1, 0:NK+1] = Tr[it, 0:NI+1, 0:NJ+1, 0:NK+1, m]
                end
            end

            advect(var, uvarx)
            #@ccall "./PSOM_LIB.so".advect_(pointer(var)::Ptr{rc_kind}, pointer(uvarx)::Ptr{rc_kind})::Cvoid

            if (selvar < 5.5)
                vardif .= 0.0
                mixing_horizontal(var, vardif)
                #@ccall "./PSOM_LIB.so".mixing_horizontal_(pointer(var)::Ptr{rc_kind}, pointer(vardif)::Ptr{rc_kind})::Cvoid
                @views @. uvarx[1:NI, 1:NJ, 1:NK] -= vardif[1:NI, 1:NJ, 1:NK]
            end

            @views @. uvarx_advec[1:NI, 1:NJ, 1:NK] = uvarx[1:NI, 1:NJ, 1:NK]

            if (selvar < 4.5)
                vardif .= 0.0
                mixing_vertical(var, vardif, m, step, iv_compute_kz, selvar)
                #@ccall "./PSOM_LIB.so".mixing_vertical_(pointer(var)::Ptr{rc_kind}, pointer(vardif)::Ptr{rc_kind}, Ref(m)::Ptr{Int}, Ref(step)::Ptr{Int}, Ref(iv_compute_kz)::Ptr{Int})::Cvoid
                @views @. uvarx[1:NI, 1:NJ, 1:NK] -= vardif[1:NI, 1:NJ, 1:NK]
                iv_compute_kz = 0
            end

            if (selvar == 3)
                @views @. uvarx[1:NI, 1:NJ, 1:NK] -= u_w[1:NI, 1:NJ, 1:NK]
                mat_D[1:NI, 1:NJ, 1:NK, selvar] = u_w[1:NI, 1:NJ, 1:NK]

            end
            if (selvar == 4)
                @views @. uvarx[1:NI, 1:NJ, 1:NK] -= v_w[1:NI, 1:NJ, 1:NK]
                mat_D[:, :, :, selvar] = v_w[1:NI, 1:NJ, 1:NK]
            end

            if (selvar == 3 || selvar == 4)
                @views @. uvarx[1:NI, 1:NJ, 1] .+= RR * (1e0 / (UL * delta)) * (Jac[1:NI, 1:NJ, 1] * wz[1:NI, 1:NJ, 1]) * var[1:NI, 1:NJ, 1]
            end

            @static if (cppdefs.gotm_im)
                if (selvar < 6)
                    for i in 1:NI
                        for j in 1:NJ
                            mat_A[i, j, 1:NK, selvar] = (-1e0) * dtimel * Jacinv[i, j, 1:NK] * mat_A[i, j, 1:NK, selvar]
                            mat_B[i, j, 1:NK, selvar] = (1e0 - (dtimel * (Jacinv[i, j, 1:NK] * mat_B[i, j, 1:NK, selvar])))
                            mat_C[i, j, 1:NK-1, selvar] = (-1e0) * (dtimel * (Jacinv[i, j, 1:NK-1] * mat_C[i, j, 1:NK-1, selvar]))
                            mat_D[i, j, 1:NK, selvar] = dtimel * (Jacinv[i, j, 1:NK] * (mat_D[i, j, 1:NK, selvar] - uvarx_advec[i, j, 1:NK])) + var[i, j, 1:NK]
                            solve_tridiag(view(mat_A, i, j, 1:NK, selvar), view(mat_B, i, j, 1:NK, selvar), view(mat_C, i, j, 1:NK, selvar), view(mat_D, i, j, 1:NK, selvar), view(mat_test, i, j, 1:NK, selvar), NK)
                            #@ccall "./PSOM_LIB.so".solve_tridiag_(pointer(mat_A[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_B[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_C[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_D[i,j,1:NK, selvar])::Ptr{rc_kind},pointer(mat_test[i,j,1:NK, selvar])::Ptr{rc_kind}, Ref(NK)::Ref{rc_kind})::Cvoid
                        end
                    end
                end
                if (selvar == 1)
                    T[1:NI, 1:NJ, 1:NK, n] = mat_test[:, :, :, selvar]
                end
                if (selvar == 2)
                    s[1:NI, 1:NJ, 1:NK, n] = mat_test[:, :, :, selvar]
                end
                if (selvar == 3)
                    cx[1:NI, 1:NJ, 1:NK] = mat_test[:, :, :, selvar]
                end
                if (selvar == 4)
                    cy[1:NI, 1:NJ, 1:NK] = mat_test[:, :, :, selvar]
                end
                if (selvar == 5)
                    cz[1:NI, 1:NJ, 1:NK] = mat_test[:, :, :, selvar]
                end

            else
                local lc_diff = @views @. dtimel * JacInv[1:NI, 1:NJ, 1:NK] .* uvarx[1:NI, 1:NJ, 1:NK]
                if (selvar == 1)
                    @views @. T[1:NI, 1:NJ, 1:NK, n] = T[1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
                if (selvar == 2)
                    @views @. s[1:NI, 1:NJ, 1:NK, n] = s[1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
                if (selvar == 3)
                    @views @. cx[1:NI, 1:NJ, 1:NK] = u[1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
                if (selvar == 4)
                    @views @. cy[1:NI, 1:NJ, 1:NK] = v[1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
                if (selvar == 5)
                    @views @. cz[1:NI, 1:NJ, 1:NK] = w[1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
            end

            local lc_diff = @views @. dtimel * JacInv[1:NI, 1:NJ, 1:NK] .* uvarx[1:NI, 1:NJ, 1:NK]

            @views @. cz[1:NI, 1:NJ, 1:NK] = w[1:NI, 1:NJ, 1:NK, 0] - lc_diff
            for it in 1:ntr
                if (selvar == 5 + it)
                    @views @. Tr[it, 1:NI, 1:NJ, 1:NK, n] = Tr[it, 1:NI, 1:NJ, 1:NK, 0] - lc_diff
                end
            end
        end
    end

end