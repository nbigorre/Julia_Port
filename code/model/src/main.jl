include("../inc/includeFile.jl")

using .fortVar
using .cppdefs


#preproc relaxation
@static if (cppdefs.relaxation)
    include("relaxation.jl")
    using .relaxation
end

#preproc allow_particle
@static if (cppdefs.allow_particle)
    include("particles.jl")
    using .particles
end

# Main decl                                      
function main()
    step = 0::Int

    fdiv = 0e0::rc_kind
    ctrdiv = 0e0::rc_kind


    pcorr = zeros(rc_kind, maxout)

    include("../inc/ini_param.jl")

    @static if (cppdefs.relaxation)
        set_coef()
    end

    @static if (cppdefs.allow_particle)
        # ALLOCATE ARRAY
    end

    ini_setup(pcorr)
    #@ccall "./PSOM_LIB.so".ini_setup_(pointer(pcorr)::Ptr{rc_kind})::Cvoid

    global hmean = meanh(NI, NJ, h)
    @fortSet("hmean", hmean, rc_kind)
    #@ccall "./PSOM_LIB.so".meanh_(Ref(NI)::Ptr{Int}, Ref(NJ)::Ptr{Int}, @lkGet("h", rc_kind)::Ptr{rc_kind}, @lkGet("hmean", rc_kind)::Ptr{rc_kind})::Cvoid

    sigma()
    #@ccall "./PSOM_LIB.so".sigma_()::Cvoid

    staticsigma()
    #@ccall "./PSOM_LIB.so".staticsigma_()::Cvoid

    tracerinit(0)
    #@ccall "./PSOM_LIB.so".tracerinit_(Ref(0)::Ptr{Int})::Cvoid

    hsave()
    #@ccall "./PSOM_LIB.so".hsave_()::Cvoid

    ini_uv(0)
    #@ccall "./PSOM_LIB.so".ini_uv_(Ref(0)::Ptr{Int})::Cvoid

    fdiv = facediv(EPS)
    #@ccall "./PSOM_LIB.so".facediv_(Ref(EPS)::Ptr{rc_kind}, Ref(fdiv)::Ptr{rc_kind})::Cvoid

    ctrdiv = cdiv(EPS, 0)
    #@ccall "./PSOM_LIB.so".cdiv_(Ref(EPS)::Ptr{rc_kind}, Ref(ctrdiv)::Ptr{rc_kind}, Ref(0)::Ptr{Int})::Cvoid

    vort(0)
    #@ccall "./PSOM_LIB.so".vort_(Ref(0)::Ptr{Int})::Cvoid

    step = 0
    global time_nondim = dtf * step
    @fortSet("time_nondim", time_nondim, rc_kind)
    global time_seconds = time_nondim * TL
    @fortSet("time_seconds", time_seconds, rc_kind)

    @static if (cppdefs.gotm_call)
        @ccall "./PSOM_LIB.so".initial_tke_()::Cvoid
    end
    @static if (cppdefs.implicit)
        println("mixing implicit")
    else
        println("mixing explicit")
    end

    diag_n2()
    #@ccall "./PSOM_LIB.so".diag_n2_()::Cvoid
    
    diag_n2budget(step)
    #@ccall "./PSOM_LIB.so".diag_n2budget_(Ref(step)::Ptr{Int})::Cvoid

    #tim = dtime

    setbc(step)
    #@ccall "./PSOM_LIB.so".setbc_(Ref(step)::Ptr{Int})::Cvoid
    
    correctbc()
    #@ccall "./PSOM_LIB.so".correctbc_()::Cvoid

    checks()
    #@ccall "./PSOM_LIB.so".checks_()::Cvoid

    local lv_test_output_bin = cppdefs.file_output && cppdefs.file_output_bin

    local initial_step = 0
    if (pickup_step < -0.5)
        initial_step = 0
        println("No Pickup")
    else
        if (! lv_test_output_bin)
            error("Error: A pickup is required but the binary I/O is not activated")
            exit(1)
        else
            println("The simulation will start from step ", pickup_step)
            initial_step = pickup_step
            step = initial_step
        end

    end

    @static if (cppdefs.file_output)
        @static if (cppdefs.file_output_cdf)
            write_cdf(step, 0)
            #@ccall "./PSOM_LIB.so".write_cdf_(Ref(step)::Ptr{Int}, Ref(0)::Ptr{Int})::Cvoid
        end
        @static if (cppdefs.file_output_bin)
            #@ccall "./PSOM_LIB.so".write_bin_(Ref(step)::Ptr{Int})::Cvoid
        end
    end

    diag_energy(step)
    #@ccall "./PSOM_LIB.so".diag_energy_(Ref(step)::Ptr{Int})::Cvoid


    #########################################################
    ##               END OF INITIALISATION                 ##
    #########################################################

    istep = initial_step
    nstep = nsteps

    for step in (istep+1:istep+nstep)
        global time_nondim = dtf * step
        @fortSet("time_nondim", time_nondim, rc_kind)
        global time_seconds = time_nondim * TL
        @fortSet("time_seconds", time_seconds, rc_kind)

        #preproc allow_particle
        @static if (cppdefs.allow_particle)
            #                       IMPLEM JULIA
        end

        diag_n2()
        #@ccall "./PSOM_LIB.so".diag_n2_()::Cvoid

        #@ccall "./PSOM_LIB.so".momentum_(pointer(pcorr)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid
        momentum(pcorr, step)


        diag_n2budget(step)
        #@ccall "./PSOM_LIB.so".diag_n2budget_(Ref(step)::Ptr{Int})::Cvoid

        diag_energy(step)
        #@ccall "./PSOM_LIB.so".diag_energy_(Ref(step)::Ptr{Int})::Cvoid

        global hmean = meanh(NI, NJ, h)
        @fortSet("hmean", hmean, rc_kind)
        #@ccall "./PSOM_LIB.so".meanh_(Ref(NI)::Ptr{Int}, Ref(NJ)::Ptr{Int}, @lkGet("h", rc_kind)::Ptr{rc_kind}, @lkGet("hmean", rc_kind)::Ptr{rc_kind})::Cvoid

        @static if (cppdefs.allow_particle)
            #                       IMPLEM JULIA
        end

        @static if (cppdefs.file_output)
            @static if (cppdefs.file_output_cdf)
                write_cdf(step, 0)
                #@ccall "./PSOM_LIB.so".write_cdf_(Ref(step)::Ptr{Int}, Ref(0)::Ptr{Int})::Cvoid
            end
            @static if (cppdefs.file_output_bin)
                #@ccall "./PSOM_LIB.so".write_bin_(Ref(step)::Ptr{Int})::Cvoid
            end
        end

    end
end

#main()
@time main()
#@profview main()
