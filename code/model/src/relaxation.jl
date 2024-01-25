module relaxation

    using ..fortVar


    local r_T = OffsetArrays.zeros(rc_kind, 0:NJ+1, 0:NK+1)
    local varbar = OffsetArrays.zeros(rc_kind, 0:NJ+1, 0:NK+1)

    function set_coef()
        local gamma_T = 1e0 / (5e0 * 86400e0 * UL / LEN)
        for k in 0:NK+1
            @views r_T[:, k] .= gamma_T .* r_sponge
        end
    end

    function sponge(n, dtime)
        for j in 0:NJ+1
            for k in 0:NK+1
                varbar[j,k] = sum(T[1:NI, j, k, n])
            end
        end

        for i in 0:NI+1
            for j in 0:NJ+1
                for k in 0:NK+1
                    T[i, j, k, n] -= dtime * r_T[j, k] * (T[i, j, k, n] - T_ref[i, j, k])
                end
            end
        end
    end

    export set_coef, sponge
end