module relaxation

    using ..header : NI, NJ, NK,  LEN, T, T_ref, r_sponge, rc_kind
    using ..fortVar


    local r_T = zeros(rc_kind, (NJ+2, NK+2))
    local varbar = zeros(rc_kind, (NJ+2, NK+2))

    function set_coef()
        local gamma_T = 1e0 / (5e0 * 86400e0 * @fortGet("ul", rc_kind) / LEN)
        r_T .= gamma_T * r_sponge
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
                    T[i, j, k, n] -= dtime * r_T[j+1, k+1] * (T[i, j, k, n] - T_ref[i, j, k])
                end
            end
        end
    end

    export set_coef, sponge
end