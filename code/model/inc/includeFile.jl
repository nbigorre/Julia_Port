if (!@isdefined(INCLUDE_All))
    global INCLUDE_All = true

    import Pkg
    if (Base.find_package("ExportAll") === nothing)
        Pkg.add("ExportAll")
    end

    include("../inc/fortVar.jl")
    
    include("../inc/cppdefs.jl")
    include("../src/header.jl")
    include("../src/momentum.jl")
    include("../src/cfdiv.jl")
    include("../src/advection_and_mixing.jl")
    include("../src/advect.jl")
    include("../src/mixing_horizontal.jl")
    include("../src/mixing_vertical.jl")
end