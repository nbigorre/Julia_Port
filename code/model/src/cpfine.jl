function cpfine(dtimel,cpf,fn)
    cpf = reshape(view(cpf, 1:(19*NI*NJ*NK)), 19, NI, NJ, NK)
    fn = reshape(view(fn, 1:(NI*NJ*NK)), NI, NJ, NK)
    local edt = EPS / dtimel
    
    pbc(cpf, fn, dtimel)
    #@ccall "./PSOM_LIB.so".pbc_(pointer(cpf)::Ptr{rc_kind}, pointer(fn)::Ptr{rc_kind}, Ref(dtimel)::Ref{rc_kind})::Cvoid

    for i in 1:NI
        for j in 2:NJ-1
            for k in 2:NK-1
                fn[i, j, k] = edt * (cxf[i, j, k] - cxf[i-1, j, k] + cyf[i, j, k] - cyf[i, j-1, k] + czf[i, j, k] - czf[i, j, k-1]
                                    - (skfc[i, j, k] - skfc[i, j, k-1])
                                    + (sifc[i, j, k] - sifc[i-1, j, k])
                                    + (sjfc[i, j, k] - sjfc[i, j-1, k]))
                cpf[1, i, j, k] = -(gqi[i, j, k, 1] + gqi[i+1, j, k, 1] + gqj[i, j, k, 2] + gqj[i, j+1, k, 2] + gqk[i+1, j+1, k, 3] + gqk[i+1, j+1, k+1, 3])
                cpf[2, i, j, k] = gqi[i+1, j, k, 1] + 0.25e0 * (gqj[i, j+1, k, 1] - gqj[i, j, k, 1] + gqk[i+1, j+1, k+1, 1] - gqk[i+1, j+1, k, 1])
                cpf[3, i, j, k] = gqj[i, j+1, k, 2] + 0.25e0 * (gqi[i+1, j, k, 2] - gqi[i, j, k, 2] + gqk[i+1, j+1, k+1, 2] - gqk[i+1, j+1, k, 2])
                cpf[4, i, j, k] = gqi[i, j, k, 1] + 0.25e0 * (-gqj[i, j+1, k, 1] + gqj[i, j, k, 1] - gqk[i+1, j+1, k+1, 1] + gqk[i+1, j+1, k, 1])
                cpf[5, i, j, k] = gqj[i, j, k, 2] + 0.25e0 * (-gqi[i+1, j, k, 2] + gqi[i, j, k, 2] - gqk[i+1, j+1, k+1, 2] + gqk[i+1, j+1, k, 2])
                cpf[6, i, j, k] = gqk[i+1, j+1, k+1, 3] + 0.25e0 * (gqi3[i, j, k] - gqi3[i-1, j, k] + gqj3[i, j, k] - gqj3[i, j-1, k])
                cpf[7, i, j, k] = gqk[i+1, j+1, k, 3] + 0.25e0 * (- gqi3[i, j, k] + gqi3[i-1, j, k] - gqj3[i, j, k] - gqj3[i, j-1, k])
                cpf[8, i, j, k] = 0.25e0 * (-gqi[i, j, k, 2] - gqj[i, j+1, k, 1])
                cpf[9, i, j, k] = 0.25e0 * (gqi[i, j, k, 2] + gqj[i, j, k, 1])
                cpf[10, i, j, k] = 0.25e0 * (-gqi[i+1, j, k, 2] - gqj[i, j, k, 1])
                cpf[11, i, j, k] = 0.25e0 * (gqi[i+1, j, k, 2] + gqj[i, j+1, k, 1])
                cpf[12, i, j, k] = 0.25e0 * (gqi3[i-1, j, k] + gqk[i+1, j+1, k, 1])
                cpf[13, i, j, k] = 0.25e0 * (-gqi3[i, j, k] - gqk[i+1, j+1, k, 1])
                cpf[14, i, j, k] = 0.25e0 * (gqi3[i, j, k] + gqk[i+1, j+1, k+1, 1])
                cpf[15, i, j, k] = 0.25e0 * (-gqi3[i-1, j, k] - gqk[i+1, j+1, k+1, 1])
                cpf[16, i, j, k] = 0.25e0 * (gqj3[i, j-1, k] + gqk[i+1, j+1, k, 2])
                cpf[17, i, j, k] = 0.25e0 * (-gqj3[i, j, k] - gqk[i+1, j+1, k, 2])
                cpf[18, i, j, k] = 0.25e0 * (gqj3[i, j, k] + gqk[i+1, j+1, k+1, 2])
                cpf[19, i, j, k] = 0.25e0 * (-gqj3[i, j-1, k] - gqk[i+1, j+1, k+1, 2])
            end
        end
    end

end