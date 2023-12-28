function mgpfill(dtloc, pf)
      local edt = 0e0
      if (dtloc > 1e-16)
            edt = EPS / dtloc
      end
      if (@fortGet("qpr", rc_kind) != 0e0)


            local k = 0
            for i in 1:NI
                  for j in 1:NJ
                        local gin = 1e0 / gqk[i+1, j+1, k+1, 3]
                        local dpdxi = 0.25e0 * (pf[i+2, j+1, k+2] + pf[i+2, j+1, k+1] - pf[i, j+1, k+2] - pf[i, j+1, k+1])
                        local dpdelta = 0.25e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+2, k+1] - pf[i+1, j, k+2] - pf[i+1, j, k+1])
                        local psig = gin * (-skfc[i, j, k] - gqk[i+1, j+1, k+1, 1] * dpdxi - gqk[i+1, j+1, k+1, 2] * dpdelta + edt * (czf[i, j, k] - wfbcb[i, j]))
                        pf[i+1, j+1, 1] = pf[i+1, j+1, 2] - psig
                  end
            end

            local j = 0
            for i in 1:NJ
                  local gin = 1e0 / gqk[i+1, j+1, k+1, 3]
                  local dpdxi = 0.25e0 * (pf[i+2, j+1, k+2] + pf[i+2, j+1, k+1] - pf[i, j+1, k+2] - pf[i, j+1, k+1])
                  local dpdelta = 0.5e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+2, k+1] - pf[i+1, j+1, k+2] - pf[i+1, j+1, k+1])
                  local psig = gin * (-skfc[i, j, k] - gqk[i+1, j+1, k+1, 1] * dpdxi - gqk[i+1, j+1, k+1, 2] * dpdelta)
                  pf[i+1, j+1, 1] = pf[i+1, j+1, 2] - psig
            end
            local j = NJ + 1
            for i in 1:NJ
                  local gin = 1e0 / gqk[i+1, j+1, k+1, 3]
                  local dpdxi = 0.25e0 * (pf[i+2, j+1, k+2] + pf[i+2, j+1, k+1] - pf[i, j+1, k+2] - pf[i, j+1, k+1])
                  local dpdelta = 0.5e0 * (pf[i+1, j+1, k+2] + pf[i+1, j+1, k+1] - pf[i+1, j, k+2] - pf[i+1, j, k+1])
                  local psig = gin * (-skfc[i, j, k] - gqk[i+1, j+1, k+1, 1] * dpdxi - gqk[i+1, j+1, k+1, 2] * dpdelta)
                  pf[i+1, j+1, 1] = pf[i+1, j+1, 2] - psig
            end

            local k = NK
            for i in 1:NI
                  for j in 0:NJ+1
                        pf[i+1, j+1, NK+2] = -pf[i+1, j+1, NK+1]
                  end
            end



            local j = 0
            for i in 1:NI
                  for k in 1:NK
                        local gin = 1e0 / gqj[i, j+1, k, 2]
                        local dpdxi = 0.25e0 * (pf[i+2, j+2, k+1] + pf[i+2, j+1, k+1] - pf[i, j+2, k+1] - pf[i, j+1, k+1])
                        local dpdsig = 0.25e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+1, k+2] - pf[i+1, j+2, k] - pf[i+1, j+1, k])
                        local peta = gin * (-gqj[i, j+1, k, 1] * dpdxi - gqj3[i, j, k] * dpdsig + (edt * (cyf[i, j, k] - vfbcs[i, k]) - sjfc[i, j, k]))
                        pf[i+1, 1, k+1] = pf[i+1, 2, k+1] - peta
                  end
            end


            local j = NJ
            for i in 1:NI
                  for k in 1:NK
                        local gin = 1e0 / gqj[i, j+1, k, 2]
                        local dpdxi = 0.25e0 * (pf[i+2, j+2, k+1] + pf[i+2, j+1, k+1] - pf[i, j+2, k+1] - pf[i, j+1, k+1])
                        local dpdsig = 0.25e0 * (pf[i+1, j+2, k+2] + pf[i+1, j+1, k+2] - pf[i+1, j+2, k] - pf[i+1, j+1, k])
                        local peta = gin * (-gqj[i, j+1, k, 1] * dpdxi - gqj3[i, j, k] * dpdsig + (edt * (cyf[i, j, k] - vfbcn[i, k]) - sjfc[i, j, k]))
                        pf[i+1, NJ+2, k+1] = pf[i+1, NJ+1, k+1] + peta
                  end
            end
      else
            local k = NK
            for i in 1:NI
                  for j in 0:NJ+1
                        pf[i+1, j+1, NK+2] = -pf[i+1, j+1, NK+1]
                  end
            end
      end

      local i = 0
      for j in 0:NJ+1
            for k in 0:NK+1
                  pf[1, j+1, k+1] = pf[NI+1, j+1, k+1]
                  pf[NI+2, j+1, k+1] = pf[2, j+1, k+1]
            end
      end

end