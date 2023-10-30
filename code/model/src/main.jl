include("../inc/includeFile.jl")

using .fortVar
using .header
using .cppdefs


#preproc relaxation
if (cppdefs.relaxation)
    include("relaxation.jl")
end

#preproc allow_particle
if (cppdefs.allow_particle)
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

    if (cppdefs.relaxation)
        set_coef()
    end

    if (cppdefs.allow_particle)
        # ALLOCATE ARRAY
    end

    @ccall "./PSOM_LIB.so".ini_setup_(pointer(pcorr)::Ptr{rc_kind})::Cvoid

    @ccall "./PSOM_LIB.so".meanh_(Ref(NI)::Ptr{Int}, Ref(NJ)::Ptr{Int}, @lkGet("h", rc_kind)::Ptr{rc_kind}, @lkGet("hmean", rc_kind)::Ptr{rc_kind})::Cvoid

    @ccall "./PSOM_LIB.so".sigma_()::Cvoid
    @ccall "./PSOM_LIB.so".staticsigma_()::Cvoid

    @ccall "./PSOM_LIB.so".tracerinit_(Ref(0)::Ptr{Int})::Cvoid

    @ccall "./PSOM_LIB.so".hsave_()::Cvoid

    @ccall "./PSOM_LIB.so".ini_uv_(Ref(0)::Ptr{Int})::Cvoid

    @ccall "./PSOM_LIB.so".facediv_(Ref(EPS)::Ptr{rc_kind}, Ref(fdiv)::Ptr{rc_kind})::Cvoid

    @ccall "./PSOM_LIB.so".cdiv_(Ref(EPS)::Ptr{rc_kind}, Ref(ctrdiv)::Ptr{rc_kind}, Ref(0)::Ptr{Int})::Cvoid

    @ccall "./PSOM_LIB.so".vort_(Ref(0)::Ptr{Int})::Cvoid

    step = 0
    @fortSet("time_nondim", @fortGet("dtf", rc_kind) * step, rc_kind)
    @fortSet("time_seconds", @fortGet("time_nondim", rc_kind) * @fortGet("tl", rc_kind), rc_kind)

    if (cppdefs.gotm_call)
        @ccall "./PSOM_LIB.so".initial_tke_()::Cvoid
    end
    if (cppdefs.implicit)
        println("mixing implicit")
    else
        println("mixing explicit")
    end

    @ccall "./PSOM_LIB.so".diag_n2_()::Cvoid
    @ccall "./PSOM_LIB.so".diag_n2budget_(Ref(step)::Ptr{Int})::Cvoid

    #tim = dtime

    @ccall "./PSOM_LIB.so".setbc_(Ref(step)::Ptr{Int})::Cvoid
    @ccall "./PSOM_LIB.so".correctbc_()::Cvoid

    @ccall "./PSOM_LIB.so".checks_()::Cvoid

    @fortSet("lv_test_output_bin", cppdefs.file_output && cppdefs.file_output_bin, Int)


    if (@fortGet("pickup_step", Int) < -0.5)
        @fortSet("initial_step", 0, Int)
        println("No Pickup")
    else
        if (! @fortGet("lv_test_output_bin", Bool))
            error("Error: A pickup is required but the binary I/O is not activated")
            exit(1)
        else
            println("The simulation will start from step ", @fortGet("pickup_step", Int))
            @fortSet("initial_step", @fortGet("pickup_step", Int), Int)
            step = @fortGet("initial_step", Int)
        end

    end

    if (cppdefs.file_output)
        if (cppdefs.file_output_cdf)
            @ccall "./PSOM_LIB.so".write_cdf_(Ref(step)::Ptr{Int}, Ref(0)::Ptr{Int})::Cvoid
        end
        if (cppdefs.file_output_bin)
            @ccall "./PSOM_LIB.so".write_bin_(Ref(step)::Ptr{Int})::Cvoid
        end
    end

    @ccall "./PSOM_LIB.so".diag_energy_(Ref(step)::Ptr{Int})::Cvoid


    #########################################################
    ##               END OF INITIALISATION                 ##
    #########################################################

    istep = @fortGet("initial_step", Int)
    nstep = @fortGet("nsteps", Int)

    for step in (istep+1:istep+nstep)
        @fortSet("time_nondim", @fortGet("dtf", rc_kind) * step, rc_kind)
        @fortSet("time_seconds", @fortGet("time_nondim", rc_kind) * @fortGet("tl", rc_kind), rc_kind)

        #preproc allow_particle
        if (cppdefs.allow_particle)
            #                       IMPLEM JULIA
        end

        @ccall "./PSOM_LIB.so".diag_n2_()::Cvoid

        #@ccall "./PSOM_LIB.so".momentum_(pointer(pcorr)::Ptr{rc_kind}, Ref(step)::Ptr{Int})::Cvoid
        momentum(pcorr, step)


        @ccall "./PSOM_LIB.so".diag_n2budget_(Ref(step)::Ptr{Int})::Cvoid

        @ccall "./PSOM_LIB.so".diag_energy_(Ref(step)::Ptr{Int})::Cvoid

        @ccall "./PSOM_LIB.so".meanh_(Ref(NI)::Ptr{Int}, Ref(NJ)::Ptr{Int}, @lkGet("h", rc_kind)::Ptr{rc_kind}, @lkGet("hmean", rc_kind)::Ptr{rc_kind})::Cvoid

        if (cppdefs.allow_particle)
            #                       IMPLEM JULIA
        end

        if (cppdefs.file_output)
            if (cppdefs.file_output_cdf)
                @ccall "./PSOM_LIB.so".write_cdf_(Ref(step)::Ptr{Int}, Ref(0)::Ptr{Int})::Cvoid
            end
            if (cppdefs.file_output_bin)
                @ccall "./PSOM_LIB.so".write_bin_(Ref(step)::Ptr{Int})::Cvoid
            end
        end

    end
end

#main()
@time main()
#@profview main()
