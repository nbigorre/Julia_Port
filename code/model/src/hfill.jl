function hfill(dtimel,hf, sums) 
      local edt = EPS / dtimel
      local dnk = rc_kind(NK)
      local cons = (1e0 - @fortGet("kappah", rc_kind)) * gpr
      local gkn = gpr * @fortGet("kappah", rc_kind)

      for i in 1:NI
            local hxi = 0.25e0 * (hf[i+2, 1] + hf[i+2, 2] - hf[i, 1] - hf[i, 2])
            @views local sumcyf = sums.sumcyfj[i]
            @views local sumvf = sums.sumvf[i]
            @views local sumsjf = sums.sumsjfj[i]
            @views local sumhyn = sums.sumhynj[i]
            @views local sumgj = sums.sumgjj[i]
            @views local gradh = (edt * (sumcyf - sumvf) - sumsjf - cons * sumhyn) / gkn
            hf[i+1, 1] = hf[i+1, 2] + (sumgj[1] * hxi - gradh) / sumgj[2]
      end


      for i in 1:NI
            local hxi = 0.25e0 * (hf[i+2, NJ+1] + hf[i+2, NJ+2] - hf[i, NJ+1] - hf[i, NJ+2])
            @views local sumcyf = sums.sumcyfnj[i]
            @views local sumvf =  sums.sumvf[i]
            @views local sumsjf = sums.sumsjfnj[i]
            @views local sumhyn = sums.sumhynnj[i]
            @views local sumgj =  sums.sumgjnj[i]
            @views local gradh = (edt * (sumcyf - sumvf) - sumsjf - cons * sumhyn) / gkn
            hf[i+1, NJ+2] = hf[i+1, NJ+1] + (-sumgj[1] * hxi + gradh) / sumgj[2]
      end
      
      hf[1, 1:NJ+2] = hf[NI+1, 1:NJ+2]
      hf[NI+2, 1:NJ+2] = hf[2, 1:NJ+2]

      
end