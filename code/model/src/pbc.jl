function pbc(cpf, fn, dtimel)
    local edt = EPS / dtimel
    for k in 1:NK
        for j in 1:NJ
            for i in 1:NI
                if (j == 1 || j == NJ || k == 1 || k == NK)
                    fn[i, j, k] = 0e0
                    cpf[1:19, i, j, k] .= 0e0

                    fn[i, j, k] += edt * cxf[i, j, k] - sifc[i, j, k]
                    cpf[1, i, j, k] -= gqi[i, j, k, 1]
                    cpf[2, i, j, k] += gqi[i, j, k, 1]
                    cpf[3, i, j, k] += 0.25e0 * gqi[i, j, k, 2]
                    cpf[5, i, j, k] -= 0.25e0 * gqi[i, j, k, 2]
                    cpf[6, i, j, k] += 0.25e0 * gqi3[i, j, k]
                    cpf[7, i, j, k] -= 0.25e0 * gqi3[i, j, k]
                    cpf[10, i, j, k] -= 0.25e0 * gqi[i, j, k, 2]
                    cpf[11, i, j, k] += 0.25e0 * gqi[i, j, k, 2]
                    cpf[13, i, j, k] -= 0.25e0 * gqi3[i, j, k]
                    cpf[14, i, j, k] += 0.25e0 * gqi3[i, j, k]

                    fn[i, j, k] += -edt * cxf[i-1, j, k] + sifc[i-1, j, k]
                    cpf[1, i, j, k] = cpf[1, i, j, k] - gqi[i-1, j, k, 1]
                    cpf[3, i, j, k] = cpf[3, i, j, k] - 0.25e0 * gqi[i-1, j, k, 2]
                    cpf[4, i, j, k] = cpf[4, i, j, k] + gqi[i-1, j, k, 1]
                    cpf[5, i, j, k] = cpf[5, i, j, k] + 0.25e0 * gqi[i-1, j, k, 2]
                    cpf[6, i, j, k] = cpf[6, i, j, k] - 0.25e0 * gqi3[i-1, j, k]
                    cpf[7, i, j, k] = cpf[7, i, j, k] + 0.25e0 * gqi3[i-1, j, k]
                    cpf[8, i, j, k] = cpf[8, i, j, k] - 0.25e0 * gqi[i-1, j, k, 2]
                    cpf[9, i, j, k] = cpf[9, i, j, k] + 0.25e0 * gqi[i-1, j, k, 2]
                    cpf[12, i, j, k] = cpf[12, i, j, k] + 0.25e0 * gqi3[i-1, j, k]
                    cpf[15, i, j, k] = cpf[15, i, j, k] - 0.25e0 * gqi3[i-1, j, k]

                    if (j == NJ)
                        fn[i, j, k] += edt * vfbcn[i, k]
                    else
                        fn[i, j, k] += edt * cyf[i, j, k] - sjfc[i, j, k]
                        cpf[1, i, j, k] = cpf[1, i, j, k] -            gqj[i, j, k, 2]
                        cpf[2, i, j, k] = cpf[2, i, j, k] + 0.25e0 *   gqj[i, j, k, 1]
                        cpf[3, i, j, k] = cpf[3, i, j, k] +            gqj[i, j, k, 2]
                        cpf[4, i, j, k] = cpf[4, i, j, k] - 0.25e0 *   gqj[i, j, k, 1]
                        cpf[6, i, j, k] = cpf[6, i, j, k] + 0.25e0 *  gqj3[i, j, k]
                        cpf[7, i, j, k] = cpf[7, i, j, k] - 0.25e0 *  gqj3[i, j, k]
                        cpf[8, i, j, k] = cpf[8, i, j, k] - 0.25e0 *   gqj[i, j, k, 1]
                        cpf[11, i, j, k] = cpf[11, i, j, k] + 0.25e0 * gqj[i, j, k, 1]
                        cpf[17, i, j, k] = cpf[17, i, j, k] - 0.25e0 *gqj3[i, j, k]
                        cpf[18, i, j, k] = cpf[18, i, j, k] + 0.25e0 *gqj3[i, j, k]
                    end

                    if (j == 1)
                        fn[i, j, k] -= edt * vfbcs[i, k]
                    else
                        fn[i, j, k] += -edt * cyf[i, j-1, k] + sjfc[i, j-1, k]
                        cpf[1, i, j, k] = cpf[1, i, j, k] - gqj[i, j-1, k, 2]
                        cpf[2, i, j, k] = cpf[2, i, j, k] - 0.25e0 * gqj[i, j-1, k, 1]
                        cpf[4, i, j, k] = cpf[4, i, j, k] + 0.25e0 * gqj[i, j-1, k, 1]
                        cpf[5, i, j, k] = cpf[5, i, j, k] + gqj[i, j-1, k, 2]
                        cpf[6, i, j, k] = cpf[6, i, j, k] - 0.25e0 * gqj3[i, j-1, k]
                        cpf[7, i, j, k] = cpf[7, i, j, k] + 0.25e0 * gqj3[i, j-1, k]
                        cpf[9, i, j, k] = cpf[9, i, j, k] + 0.25e0 * gqj[i, j-1, k, 1]
                        cpf[10, i, j, k] = cpf[10, i, j, k] - 0.25e0 * gqj[i, j-1, k, 1]
                        cpf[16, i, j, k] = cpf[16, i, j, k] + 0.25e0 * gqj3[i, j-1, k]
                        cpf[19, i, j, k] = cpf[19, i, j, k] - 0.25e0 * gqj3[i, j-1, k]
                    end

                    if (k == NK)
                        fn[i, j, k] += edt * czf[i, j, k] - skfc[i, j, k]
                        cpf[1, i, j, k] -= 2e0 * gqk[i, j, k, 3]
                    else
                        fn[i, j, k] += edt * czf[i, j, k] - skfc[i, j, k]
                        cpf[1, i, j, k] = cpf[1, i, j, k] - gqk[i, j, k, 3]
                        cpf[2, i, j, k] = cpf[2, i, j, k] + 0.25e0 * gqk[i, j, k, 1]
                        cpf[3, i, j, k] = cpf[3, i, j, k] + 0.25e0 * gqk[i, j, k, 2]
                        cpf[4, i, j, k] = cpf[4, i, j, k] - 0.25e0 * gqk[i, j, k, 1]
                        cpf[5, i, j, k] = cpf[5, i, j, k] - 0.25e0 * gqk[i, j, k, 2]
                        cpf[6, i, j, k] = cpf[6, i, j, k] + gqk[i, j, k, 3]
                        cpf[14, i, j, k] = cpf[14, i, j, k] + 0.25e0 * gqk[i, j, k, 1]
                        cpf[15, i, j, k] = cpf[15, i, j, k] - 0.25e0 * gqk[i, j, k, 1]
                        cpf[18, i, j, k] = cpf[18, i, j, k] + 0.25e0 * gqk[i, j, k, 2]
                        cpf[19, i, j, k] = cpf[19, i, j, k] - 0.25e0 * gqk[i, j, k, 2]
                    end

                    if (k == 1)
                        fn[i, j, k] -= wfbcb[i, j]
                    else
                        fn[i, j, k] += -edt * czf[i, j, k-1] + skfc[i, j, k-1]
                        cpf[1, i, j, k] =  cpf[1, i, j, k]- gqk[i, j, k-1, 3]
                        cpf[2, i, j, k] =  cpf[2, i, j, k]- 0.25e0 * gqk[i, j, k-1, 1]
                        cpf[3, i, j, k] =  cpf[3, i, j, k]- 0.25e0 * gqk[i, j, k-1, 2]
                        cpf[4, i, j, k] =  cpf[4, i, j, k]+ 0.25e0 * gqk[i, j, k-1, 1]
                        cpf[5, i, j, k] =  cpf[5, i, j, k]+ 0.25e0 * gqk[i, j, k-1, 2]
                        cpf[7, i, j, k] =  cpf[7, i, j, k]+ gqk[i, j, k-1, 3]
                        cpf[12, i, j, k] = cpf[12, i, j, k] + 0.25e0 * gqk[i, j, k-1, 1]
                        cpf[13, i, j, k] = cpf[13, i, j, k] - 0.25e0 * gqk[i, j, k-1, 1]
                        cpf[16, i, j, k] = cpf[16, i, j, k] + 0.25e0 * gqk[i, j, k-1, 2]
                        cpf[17, i, j, k] = cpf[17, i, j, k] - 0.25e0 * gqk[i, j, k-1, 2]
                    end

                end
            end
        end
    end

end