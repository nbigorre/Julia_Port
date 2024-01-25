if (!@isdefined(INCLUDE_All))
    global INCLUDE_All = true

    import Pkg
    if (Base.find_package("OffsetArrays") === nothing)
        Pkg.add("OffsetArrays")
    end
    using OffsetArrays
    if (Base.find_package("NetCDF") === nothing)
        Pkg.add("NetCDF")
    end
    using NetCDF
    if (Base.find_package("YAML") === nothing)
        Pkg.add("YAML")
    end
    using YAML

    include("../inc/fortVar.jl")
    
    include("../inc/cppdefs.jl")
    include("../src/header.jl")
    include("../src/momentum.jl")
    include("../src/cfdiv.jl")
    include("../src/advection_and_mixing.jl")
    include("../src/advect.jl")
    include("../src/mixing_horizontal.jl")
    include("../src/mixing_vertical.jl")
    include("../src/heat_flux.jl")
    include("../src/wind_stress.jl")
    include("../src/intpol_3rdorder.jl")
    include("../src/sigma_toplayer.jl")
    include("../src/findz_topmoves.jl")
    include("../src/coriolis.jl")
    include("../src/calcskfc.jl")
    include("../src/evalrho.jl")
    include("../src/potdens.jl")
    include("../src/evalrho_sT.jl")
    include("../src/evalrho_rho.jl")
    include("../src/rpevalgrad_Song.jl")
    include("../src/rpevalgrad_Sche.jl")
    include("../src/cdiv.jl")
    include("../src/srcface_3rdorder.jl")
    include("../src/newcor.jl")
    include("../src/newsrc.jl")
    include("../src/hsolve.jl")
    include("../src/vhydro.jl")
    include("../src/mgrid.jl")
    include("../src/chfine.jl")
    include("../src/hbc.jl")
    include("../src/hfill.jl")
    include("../src/mprove.jl")
    include("../src/vface.jl")
    include("../src/cpfine.jl")
    include("../src/cpcors.jl")
    include("../src/pbc.jl")
    include("../src/linerelax_periodicew.jl")
    include("../src/dgtsl.jl")
    include("../src/mgpfill.jl")
    include("../src/resid.jl")
    include("../src/restrict.jl")
    include("../src/sor.jl")
    include("../src/efill.jl")
    include("../src/prolong.jl")
    include("../src/uvchy.jl")
    include("../src/vcenter.jl")
    include("../src/solve_tridiag2.jl")
    include("../src/pcorrect.jl")
    include("../src/facediv.jl")
    include("../src/conadjust_sT.jl")
    include("../src/diag_n2.jl")
    include("../src/ini_setup.jl")
    include("../src/ini_topog.jl")
    include("../src/spline.jl")
    include("../src/seval.jl")
    include("../src/smooth.jl")
    include("../src/ini_st.jl")
    include("../src/ini_h.jl")
    include("../src/staticsigma.jl")
    include("../src/findsigma.jl")
    include("../src/meanh.jl")
    include("../src/ini_tracer.jl")
    include("../src/hsave.jl")
    include("../src/ini_uv.jl")
    include("../src/srcface_nopy.jl")
    include("../src/setbc.jl")
    include("../src/correctbc.jl")
    include("../src/checks.jl")
    include("../src/write_cdf.jl")
    include("../src/write_cdf_1D_mooring.jl")
    include("../src/diag_pv.jl")
    include("../src/write_cdf_2D_sigma.jl")
    include("../src/write_cdf_2D_x.jl")
    include("../src/uv_geostrophic.jl")
    include("../src/diag_streamfunction.jl")
    include("../src/write_cdf_3D.jl")
    include("../src/write_cdf_3D_strain.jl")
    include("../src/write_bin_mod.jl")
    include("../src/write_bin.jl")
    include("../src/diag_energy.jl")
    include("../src/diag_n2budget.jl")
    include("../src/diag_vort.jl")
    include("../src/sigma2z.jl")
    
    
    
end