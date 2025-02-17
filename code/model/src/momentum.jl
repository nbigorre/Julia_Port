using .fortVar


function momentum(pcorr, step)
    oldh .= h

    rho_old .= rho
    
    local dtim::Ref{rc_kind} = Ref(0e0)
    local cfcdiv::Ref{rc_kind} = Ref(0e0)
    local edt::Ref{rc_kind} = Ref(0e0)
    local fdiv::Ref{rc_kind} = Ref(0e0)
    local ctrdiv::Ref{rc_kind} = Ref(0e0)

    for ivb in 1:3
        local ivs = ivb == 1 ? 0 : 1
        local ivf = ivb == 3 ? 0 : 1
        dtim[] = dtf / (ivb == 1 ? 3e0 : (ivb == 2 ? 2 : 1))

        
        findzall()
        #@ccall "./PSOM_LIB.so".findzall_()::Cvoid

        sigma()
        #@ccall "./PSOM_LIB.so".sigma_()::Cvoid

        advection_and_mixing(ivs[], ivf[], dtim[], step)
        #@ccall "./PSOM_LIB.so".advection_and_mixing_(ivs::Ref{Int}, ivf::Ref{Int}, dtim::Ref{rc_kind}, Ref(step)::Ref{Int})::Cvoid

        intpol()
        #@ccall "./PSOM_LIB.so".intpol_()::Cvoid

        if (use_shchepetkin == 0)
            rpevalgrad_Song(ivf[])
            #@ccall "./PSOM_LIB.so".rpevalgrad_song_(ivf::Ref{Int})::Cvoid
        else
            rpevalgrad_Song(ivf[])
            #@ccall "./PSOM_LIB.so".rpevalgrad_song_(ivf::Ref{Int})::Cvoid
            
            rpevalgrad_Sche(ivf[])
            #@ccall "./PSOM_LIB.so".rpevalgrad_sche_(ivf::Ref{Int})::Cvoid
            @. @views drpx[1:NI, 1:NJ, 1:NK] = ru4_Sche[1:NI, 1:NJ, 1:NK]

            @. @views drpy[1:NI, 1:NJ, 1:NK] = rv4_Sche[1:NI, 1:NJ, 1:NK]

            @. @views grpifc[0:NI, 1:NJ, 1:NK] = ru2_Sche[0:NI, 1:NJ, 1:NK]

            @. @views grpjfc[1:NI, 0:NJ, 1:NK] = rv2_Sche[1:NI, 0:NJ, 1:NK]
        end
        
        coriolis(ivs[])
        #@ccall "./PSOM_LIB.so".coriolis_(ivs::Ref{Int})::Cvoid

        srcface(ivs[], step)
        #@ccall "./PSOM_LIB.so".srcface_(ivs::Ref{Int}, Ref(step)::Ref{Int})::Cvoid

        hsolve(h, oldh, hdt, dtim[])
        #@ccall "./PSOM_LIB.so".hsolve_(@lkGet("h", rc_kind)::Ref{rc_kind}, @lkGet("oldh", rc_kind)::Ref{rc_kind}, @lkGet("hdt", rc_kind)::Ref{rc_kind}, dtim::Ref{rc_kind})::Cvoid
        
        calcfkfc()
        #@ccall "./PSOM_LIB.so".calcskfc_()::Cvoid
        
        vhydro(dtim[])
        #@ccall "./PSOM_LIB.so".vhydro_(dtim::Ref{rc_kind})::Cvoid
        
        cfdiv(cfcdiv)

        newcor(dtim[], ivs[])
        #@ccall "./PSOM_LIB.so".newcor_(dtim::Ref{rc_kind}, ivs::Ref{Int})::Cvoid

        newsrc()
        #@ccall "./PSOM_LIB.so".newsrc_()::Cvoid

        edt[] = EPS / dtim[]
        mgrid(pcorr, dtim[], edt[], cfcdiv[])
        #@ccall "./PSOM_LIB.so".mgrid_(pointer(pcorr)::Ptr{rc_kind}, dtim::Ref{rc_kind}, edt::Ref{rc_kind}, cfcdiv::Ref{rc_kind})::Cvoid
        
        vface(pcorr, dtim[])
        #@ccall "./PSOM_LIB.so".vface_(pointer(pcorr)::Ptr{rc_kind}, dtim::Ref{rc_kind})::Cvoid
        
        vcenter(OffsetArrays.reshape(view(pcorr, 1:(NI+2)*(NJ+2)*(NK+2)), (0:NI+1), (0:NJ+1), (0:NK+1)), dtim[], ivf[])
        #@ccall "./PSOM_LIB.so".vcenter_(pointer(pcorr)::Ptr{rc_kind}, dtim::Ref{rc_kind}, ivf::Ref{Int})::Cvoid
        
        if (fnhhy != 0e0)
            pcorrect(OffsetArrays.reshape(view(pcorr, 1:(NI+2)*(NJ+2)*(NK+2)), (0:NI+1), (0:NJ+1), (0:NK+1)))
            #@ccall "./PSOM_LIB.so".pcorrect_(pointer(pcorr)::Ptr{rc_kind})::Cvoid
        end
        fdiv[] = facediv(dtim[])
        #@ccall "./PSOM_LIB.so".facediv_(dtim::Ref{rc_kind}, fdiv::Ref{rc_kind})::Cvoid

        ctrdiv[] = cdiv(dtim[], ivf[])
        #@ccall "./PSOM_LIB.so".cdiv_(dtim::Ref{rc_kind}, ctrdiv::Ref{rc_kind}, ivf::Ref{Int})::Cvoid

    end
    evalrho(rho, 0)
    #@ccall "./PSOM_LIB.so".evalrho_(@lkGet("rho", rc_kind)::Ref{rc_kind}, Ref(0)::Ref{Int})::Cvoid
    
    conadjust(step, 0)
    #@ccall "./PSOM_LIB.so".conadjust_(Ref(step)::Ref{Int}, Ref(0)::Ref{Int})::Cvoid

    diag_n2()
    #@ccall "./PSOM_LIB.so".diag_n2_()::Cvoid
end

