global const chfine_sumsif = zeros(rc_kind, (NI,NJ))
global const chfine_sumsjf = zeros(rc_kind, (NI,NJ))
global const chfine_sumcxf = zeros(rc_kind, (NI,NJ))
global const chfine_sumcyf = zeros(rc_kind, (NI,NJ))
global const chfine_sumuf  = zeros(rc_kind, (NI,NJ))
global const chfine_sumvf  = zeros(rc_kind, (NI,NJ))
global const chfine_sumhxn = zeros(rc_kind, (NI,NJ))
global const chfine_sumhyn = zeros(rc_kind, (NI,NJ))
global const chfine_sumgi = zeros(rc_kind, (NI,NJ,2))
global const chfine_sumgj = zeros(rc_kind, (NI,NJ,2))

function chfine(dtimel,chf,fn)
   chfine_sumsjf .= 0e0
   chfine_sumcyf .= 0e0
   chfine_sumvf .= 0e0
   chfine_sumhyn .= 0e0
   chfine_sumgj .= 0e0
   chfine_sumsif .= 0e0
   chfine_sumcxf .= 0e0
   chfine_sumuf .= 0e0
   chfine_sumhxn .= 0e0
   chfine_sumgi .= 0e0
   local sumsif = chfine_sumsif
   local sumsjf = chfine_sumsjf
   local sumcxf = chfine_sumcxf
   local sumcyf = chfine_sumcyf
   local sumuf  = chfine_sumuf
   local sumvf  = chfine_sumvf
   local sumhxn = chfine_sumhxn
   local sumhyn = chfine_sumhyn
   local sumgi = chfine_sumgi
   local sumgj = chfine_sumgj

   
   local edt = EPS / dtimel
   local edtg = edt / (gpr * kappah)
   local gpkinv = 1e0 / (gpr * kappah)
   local dtinv = HDL / (dtimel * kappah)
   local constv = kaphinv - 1e0
   
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
            sumsif[i, j] += sifc[i, j, k]
            sumsjf[i, j] += sjfc[i, j, k]
            sumcxf[i, j] += cxf[i, j, k]
            sumcyf[i, j] += cyf[i, j, k]
            sumuf[i, j]  += uf[i, j, k]
            sumvf[i, j]  += vf[i, j, k]
            sumhxn[i, j] += hxn[i, j, k]
            sumhyn[i, j] += hyn[i, j, k]
            for l in 1:2
               sumgi[i, j, l] += gi[i, j, k, l]
               sumgj[i, j, l] += gj[i, j, k, l]
            end
         end
      end
   end
   
   
   
   local i = 1
   local im1 = NI
   for j in 2:NJ-1

      fn[i, j] = (edtg * (kaphinv * wfbcb[i, j] - J2d[i, j] * oldh[i, j] * dtinv
      + sumcxf[i, j] - sumcxf[im1, j] + sumcyf[i, j] - sumcyf[i, j - 1]
      + constv * (sumuf[i, j] - sumuf[im1, j] + sumvf[i, j] - sumvf[i, j-1]))
      - gpkinv * (sumsif[i, j] - sumsif[im1, j] + sumsjf[i, j] - sumsjf[i, j-1])
      - constv * (sumhxn[i, j] - sumhxn[im1, j] + sumhyn[i, j] - sumhyn[i, j-1]))
      
      chf[1, i, j] = -edtg * J2d[i, j] * dtinv - (sumgi[i, j, 1] + sumgi[im1, j, 1] + sumgj[i,j,2] + sumgj[i, j-1, 2])
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
         fn[i, j] = (edtg * (kaphinv * wfbcb[i, j] - J2d[i, j] * oldh[i, j] * dtinv
         + sumcxf[i, j] - sumcxf[i-1, j] + sumcyf[i, j] - sumcyf[i, j - 1]
         + constv * (sumuf[i, j] - sumuf[i-1, j] + sumvf[i, j] - sumvf[i, j-1]))
         - gpkinv * (sumsif[i, j] - sumsif[i-1, j] + sumsjf[i, j] - sumsjf[i, j-1])
         - constv * (sumhxn[i, j] - sumhxn[i-1, j] + sumhyn[i, j] - sumhyn[i, j-1]))
         
         chf[1, i, j] = -edtg * J2d[i, j] * dtinv - (sumgi[i, j, 1] + sumgi[i-1, j, 1] + sumgj[i,j,2] + sumgj[i, j-1, 2])
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
end