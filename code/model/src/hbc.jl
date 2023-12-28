function hbc(chf,fn,dtimel) 
  local kaph1 = 1e0 - @fortGet("kappah", rc_kind)
  local edt = EPS / dtimel
  local eg = edt / (gpr * @fortGet("kappah", rc_kind))
  local gprinv = 1e0 / (gpr * @fortGet("kappah", rc_kind))
  local constv = @fortGet("kaphinv", rc_kind) - 1e0
  local dtinv = @fortGet("hdl", rc_kind) / (dtimel * @fortGet("kappah", rc_kind))

  for i in 1:NI
    for j in 1:NJ
      if (j == 1 || j == NJ)
        fn[i, j] = eg * (@fortGet("kaphinv", rc_kind) * wfbcb[i, j] - J2d[i, j] * oldh[i+1, j+1] * dtinv)
        chf[1, i, j] = -eg * J2d[i, j] * dtinv
        chf[2:9, i, j] .= 0e0
        local sumsif = sum(sifc[i, j, 1:NK])
        local sumcxf = sum(cxf[i, j, 1:NK])
        local sumuf = sum(uf[i, j, 1:NK])
        local sumhxn = sum(hxn[i, j, 1:NK])
        local sumgi = (sum(gi[i, j, 1:NK, 1]), sum(gi[i, j, 1:NK, 2]))
        
        fn[i, j] = fn[i, j] + eg * (sumcxf + constv * sumuf) - gprinv * sumsif - constv * sumhxn
        chf[1, i, j] -= sumgi[1]
        chf[2, i, j] += sumgi[1]
        chf[4, i, j] += 0.25e0 * sumgi[2]
        chf[5, i, j] -= 0.25e0 * sumgi[2]
        chf[6, i, j] += 0.25e0 * sumgi[2]
        chf[8, i, j] -= 0.25e0 * sumgi[2]
        
        local im1 = (i == 1 ? NI : i - 1)
        
        local sumsif = sum(sifc[im1, j, 1:NK])
        local sumcxf = sum(cxf[im1, j, 1:NK])
        local sumuf = sum(uf[im1, j, 1:NK])
        local sumhxn = sum(hxn[im1, j, 1:NK])
        local sumgi = (sum(gi[im1, j, 1:NK, 1]), sum(gi[im1, j, 1:NK, 2]))
        
        fn[i, j] = fn[i, j] - eg * (sumcxf + constv * sumuf) + gprinv * sumsif + constv * sumhxn
        chf[1, i, j] -= sumgi[1]
        chf[3, i, j] += sumgi[1]
        chf[4, i, j] -= 0.25e0 * sumgi[2]
        chf[5, i, j] += 0.25e0 * sumgi[2]
        chf[7, i, j] -= 0.25e0 * sumgi[2]
        chf[9, i, j] += 0.25e0 * sumgi[2]
        
        if (j == NJ)
          local sumvf = sum(vfbcn[i, 1:NK])
          fn[i, j] = fn[i, j] + eg * sumvf * @fortGet("kaphinv", rc_kind)
        else
          local sumsjf = sum(sjfc[i, j, 1:NK])
          local sumcyf = sum(cyf[i, j, 1:NK])
          local sumvf = sum(vf[i, j, 1:NK])
          local sumhyn = sum(hyn[i, j, 1:NK])
          local sumgj = (sum(gj[i, j, 1:NK, 1]), sum(gj[i, j, 1:NK, 2]))
          fn[i, j] = fn[i, j] + eg * (sumcyf + constv *sumvf) - gprinv * sumsjf - constv * sumhyn
          chf[1, i, j] -=          sumgj[2]
          chf[2, i, j] += 0.25e0 * sumgj[1]
          chf[3, i, j] -= 0.25e0 * sumgj[1]
          chf[4, i, j] +=          sumgj[2]
          chf[6, i, j] += 0.25e0 * sumgj[1]
          chf[7, i, j] -= 0.25e0 * sumgj[1]
        end
        
        if (j == 1)
          local sumvf = sum(vfbcs[i, 1:NK])
          fn[i, j] = fn[i, j] - eg * sumvf * @fortGet("kaphinv", rc_kind)
        else
          local sumsjf = sum(sjfc[i, j-1, 1:NK])
          local sumcyf = sum(cyf[i, j-1, 1:NK])
          local sumvf = sum(vf[i, j-1, 1:NK])
          local sumhyn = sum(hyn[i, j-1, 1:NK])
          local sumgj = (sum(gj[i, j-1, 1:NK, 1]), sum(gj[i, j-1, 1:NK, 2]))
          fn[i, j] = fn[i, j] - eg * (sumcyf + constv *sumvf) + gprinv * sumsjf + constv * sumhyn
          chf[1, i, j] -=          sumgj[2]
          chf[2, i, j] -= 0.25e0 * sumgj[1]
          chf[3, i, j] += 0.25e0 * sumgj[1]
          chf[5, i, j] +=          sumgj[2]
          chf[8, i, j] -= 0.25e0 * sumgj[1]
          chf[9, i, j] += 0.25e0 * sumgj[1]
        end
      end
    end
  end

end
