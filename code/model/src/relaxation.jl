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
                varbar[j,k] = sum(T[2:NI+1, j+1, k+1, n+1])
            end
        end

        for i in 0:NI+1
            for j in 0:NJ+1
                for k in 0:NK+1
                    T[i+1, j+1, k+1, n+1] -= dtime * r_T[j+1, k+1] * (T[i+1, j+1, k+1, n+1] - T_ref[i+1, j+1, k+1])
                end
            end
        end
    end

    export set_coef, sponge
end