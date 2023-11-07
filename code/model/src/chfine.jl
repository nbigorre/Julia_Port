function chfine(dtimel,chf,fn)
   local sumsif = zeros(rc_kind, (NI,NJ))
   local sumsjf = zeros(rc_kind, (NI,NJ))
   local sumcxf = zeros(rc_kind, (NI,NJ))
   local sumcyf = zeros(rc_kind, (NI,NJ))
   local sumuf  = zeros(rc_kind, (NI,NJ))
   local sumvf  = zeros(rc_kind, (NI,NJ))
   local sumhxn = zeros(rc_kind, (NI,NJ))
   local sumhyn = zeros(rc_kind, (NI,NJ))
   local sumgi = zeros(rc_kind, (NI,NJ,2))
   local sumgj = zeros(rc_kind, (NI,NJ,2))
   
   local edt = EPS / dtimel
   local edtg = edt / (gpr * @fortGet("kappah", rc_kind))
   local gpkinv = 1e0 / (gpr * @fortGet("kappah", rc_kind))
   local dtinv = @fortGet("hdl", rc_kind) / (dtimel * @fortGet("kappah", rc_kind))
   local constv = @fortGet("kaphinv", rc_kind) - 1e0
   
   hbc(chf, fn, dtimel)
   #@ccall "./PSOM_LIB.so".hbc_(pointer(chf)::Ptr{rc_kind}, pointer(fn)::Ptr{rc_kind}, Ref(dtimel)::Ref{rc_kind})::Cvoid
   
   
   for i in 1:NI
      for j in 1:NJ-1
         sumsif[i, j] = 0e0
         sumsjf[i, j] = 0e0
         sumcxf[i, j] = 0e0
         sumcyf[i, j] = 0e0
         sumuf[i, j]  = 0e0
         sumvf[i, j]  = 0e0
         sumhxn[i, j] = 0e0
         sumhyn[i, j] = 0e0
         sumgi[i, j, :] .= 0e0
         sumgj[i, j, :] .= 0e0
         for k in 1:NK
            sumsif[i, j] += sifc[i+1, j, k]
            sumsjf[i, j] += sjfc[i, j+1, k]
            sumcxf[i, j] += cxf[i+1, j, k]
            sumcyf[i, j] += cyf[i, j+1, k]
            sumuf[i, j]  += uf[i+1, j, k]
            sumvf[i, j]  += vf[i, j+1, k]
            sumhxn[i, j] += hxn[i+1, j, k]
            sumhyn[i, j] += hyn[i, j+1, k]
            for l in 1:2
               sumgi[i, j, l] += gi[i+1, j, k, l]
               sumgj[i, j, l] += gj[i, j+1, k, l]
            end
         end
      end
   end
   
   
   
   local i = 1
   local im1 = NI
   for j in 2:NJ-1

      fn[i, j] = (edtg * (@fortGet("kaphinv", rc_kind) * wfbcb[i, j] - J2d[i+1, j+1] * oldh[i+1, j+1] * dtinv
      + sumcxf[i, j] - sumcxf[im1, j] + sumcyf[i, j] - sumcyf[i, j - 1]
      + constv * (sumuf[i, j] - sumuf[im1, j] + sumvf[i, j] - sumvf[i, j-1]))
      - gpkinv * (sumsif[i, j] - sumsif[im1, j] + sumsjf[i, j] - sumsjf[i, j-1])
      - constv * (sumhxn[i, j] - sumhxn[im1, j] + sumhyn[i, j] - sumhyn[i, j-1]))
      
      chf[1, i, j] = -edtg * J2d[i+1, j+1] * dtinv - (sumgi[i, j, 1] + sumgi[im1, j, 1] + sumgj[i,j,2] + sumgj[i, j-1, 2])
      chf[2, i, j] = sumgi[i,j,1] + 0.25e0 * ( sumgj[i,j,1] - sumgj[i,j-1,1])
      chf[3, i, j] = sumgi[im1,j,1] + 0.25e0 * (-sumgj[i,j,1] + sumgj[i,j-1,1])
      chf[4, i, j] = sumgj[i,j,2] + 0.25e0 * (-sumgi[im1,j,2] + sumgi[i,j,2])
      chf[5, i, j] = sumgj[i,j-1,2] + 0.25e0 * (-sumgi[i,j,2] + sumgi[im1,j,2])
      chf[6, i, j] = 0.25e0 * (sumgi[i,j,2] + sumgj[i,j,1])
      chf[7, i, j] = -0.25e0 * (sumgi[im1,j,2] + sumgj[i,j,1])
      chf[8, i, j] = -0.25e0 * (sumgi[i,j,2] + sumgj[i,j-1,1]) 
      chf[9, i, j] = 0.25e0 * (sumgi[im1,j,2] + sumgj[i,j-1,1])
   end
   
   for i in 2:NI
      for j in 2:NJ-1
         fn[i, j] = (edtg * (@fortGet("kaphinv", rc_kind) * wfbcb[i, j] - J2d[i+1, j+1] * oldh[i+1, j+1] * dtinv
         + sumcxf[i, j] - sumcxf[i-1, j] + sumcyf[i, j] - sumcyf[i, j - 1]
         + constv * (sumuf[i, j] - sumuf[i-1, j] + sumvf[i, j] - sumvf[i, j-1]))
         - gpkinv * (sumsif[i, j] - sumsif[i-1, j] + sumsjf[i, j] - sumsjf[i, j-1])
         - constv * (sumhxn[i, j] - sumhxn[i-1, j] + sumhyn[i, j] - sumhyn[i, j-1]))
         
         chf[1, i, j] = -edtg * J2d[i+1, j+1] * dtinv - (sumgi[i, j, 1] + sumgi[i-1, j, 1] + sumgj[i,j,2] + sumgj[i, j-1, 2])
         chf[2, i, j] = sumgi[i,j,1] + 0.25e0 * ( sumgj[i,j,1] - sumgj[i,j-1,1])
         chf[3, i, j] = sumgi[i-1,j,1] + 0.25e0 * (-sumgj[i,j,1] + sumgj[i,j-1,1])
         chf[4, i, j] = sumgj[i,j,2] + 0.25e0 * (-sumgi[i-1,j,2] + sumgi[i,j,2])
         chf[5, i, j] = sumgj[i,j-1,2] + 0.25e0 * (-sumgi[i,j,2] + sumgi[i-1,j,2])
         chf[6, i, j] = 0.25e0 * (sumgi[i,j,2] + sumgj[i,j,1])
         chf[7, i, j] = -0.25e0 * (sumgi[i-1,j,2] + sumgj[i,j,1])
         chf[8, i, j] = -0.25e0 * (sumgi[i,j,2] + sumgj[i,j-1,1]) 
         chf[9, i, j] = 0.25e0 * (sumgi[i-1,j,2] + sumgj[i,j-1,1])
      end
   end
   #=
   =#
end