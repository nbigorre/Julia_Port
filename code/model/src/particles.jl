module particles
#using ..header : NI, NJ, NK, uf, vf, wf, Jifc
# needs struct fortran translation
using ..fortVar
using ..header

struct particle
    i::rc_kind
    j::rc_kind
    k::rc_kind
    x::rc_kind
    y::rc_kind
    z::rc_kind
    u::rc_kind
    v::rc_kind
    w::rc_kind
    s::rc_kind
    t::rc_kind
    u0::rc_kind
    v0::rc_kind
    w0::rc_kind
    id::rc_kind
    vor::rc_kind
    strain::rc_kind
    rho::rc_kind
    time::rc_kind
end

parti::Array{particle} = fill(particle(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), 0)
dz::rc_kind = 0e0
swap1::Ref{rc_kind} = Ref{0e0}
swap2::Ref{rc_kind} = Ref{0e0}
swap3::Ref{rc_kind} = Ref{0e0}
file_id::Array{IOStream} = Int[]
NPR_eachfile::Int = 0
file_id_char::string = ""

function open_parti_files()
    file_id = zeros(IOStream, length(parti_file_num))
    if (mod(NPR, parti_file_num) != 0)
        println("Error: Please make sure NPR/file_num = integer in mod_particles.jl")
        println("Stop Model")
        exit(1)
    end

    NPR_eachfile = div(NPR, parti_file_num)
    local dirout = string(@fortGetArr("dirout", Char, 151))
    dirout = rstrip(dirout, ' ')
    for fi in 1:parti_file_num
        file_id_char = lpad(string(fi), 3, "0")
        file_id[fi] = open(dirout * "/op.parti-" * file_id_char * ".bin", "w+")
    end
end

function save_parti()
    println("SAVE PARTICLES")

    for i_file in 1:parti_file_num
        for i in (i_file-1)*NPR_eachfile+1:i_file*NPR_eachfile
            write(file_id[i_file], parti[i])
        end
    end
end

function ini_particles(time)
    open_parti_files()

    for i in 1:NPR
        parti[i].time = rc_kind(time)
        parti[i].u0 = 0e0
        parti[i].v0 = 0e0
        parti[i].w0 = 0e0
        parti[i].u = 0e0
        parti[i].v = 0e0
        parti[i].w = 0e0
        parti[i].t = 0e0
        parti[i].s = 0e0
        parti[i].rho = 0e0
        parti[i].vor = 0e0
        parti[i].shear = 0e0
        parti[i].strain = 0e0

        parti[i].i = rc_kind(i * rc_kind(NI) / rc_kind(NPR))
        parti[i].j = 90e0
        parti[i].k = 2e0
    end

end

function interp_trilinear(di, dj, dk, var)
    local i1 = (var[2, 1, 1] - var[1, 1, 1]) * di + var[1, 1, 1]
    local i2 = (var[2, 1, 2] - var[1, 1, 2]) * di + var[1, 1, 2]
    local i3 = (var[2, 2, 2] - var[1, 2, 2]) * di + var[1, 2, 2]
    local i4 = (var[2, 2, 1] - var[1, 2, 1]) * di + var[1, 2, 1]

    local j1 = (i3 - i2) * dj + i2
    local j2 = (i4 - i1) * dj + i1

    return (j1 - j2) * dk + j2
end


const global get_parti_vel_wzf = zeros(rc_kind, (NI + 2, NJ + 2, NK + 1))
const global get_parti_vel_wfp = zeros(rc_kind, (NI, NJ, NK + 1))
const global get_parti_vel_ufp = zeros(rc_kind, (NI + 1, NJ, NK))
const global get_parti_vel_ufp = zeros(rc_kind, (NI, NJ + 1, NK))
const global get_parti_vel_uf_ex = zeros(rc_kind, (NI + 1, NJ + 2, NK + 2))
const global get_parti_vel_vf_ex = zeros(rc_kind, (NI + 2, NJ + 1, NK + 2))
const global get_parti_vel_wf_ex = zeros(rc_kind, (NI + 2, NJ + 2, NK + 1))
function get_parti_vel(time)
    @views @. wzf = 0.5d0 * (wz[:, :, 0:NK] + wz[:, :, 1:NK+1])

    k = 0
    wfp[:, :, 1] = wf[:, :, 1] / J2d[1:NI, 1:NJ] * wzf[2:NI+1, 2:NJ+1, 1]
    for k in 1:NK
        @views @. ufp[:, :, k] = uf[:, :, k] / Jifc[:, :, k]
        @views @. vfp[:, :, k] = vf[:, :, k] / Jjfc[:, :, k]
        @views @. wfp[:, :, k+1] = wf[:, :, k+1] / J2d[1:NI, 1:NJ] * wzf[2:NI+1, 2:NJ+1, k+1]
    end
    uf_ex .= 0e0
    @views @. uf_ex[:, 2:NJ+1, 2:NK+1] = ufp
    # === vertical extrapolation
    @views @. uf_ex[:, :, NK+2] = 2 * uf_ex[:, :, NK+1] - uf_ex[:, :, NK]

    vf_ex .= 0e0
    @views @. vf_ex[2:NI+1, :, 2:NK+1] = vfp
    # === zonally periodic
    @views @. vf_ex[1, :, :] = vf_ex[NI+1, :, :]
    @views @. vf_ex[NI+2, :, :] = vf_ex[2, :, :]
    # === vertical extrapolation
    vf_ex[:, :, NK+2] = 2 * vf_ex[:, :, NK+1] - vf_ex[:, :, NK]

    wf_ex .= 0e0
    wf_ex[2:NI+1, 2:NJ+2, :] = wfp
    # === zonally periodic
    wf_ex[1, :, :] = wf_ex[NI+1, :, :]
    wf_ex[NI+2, :, :] = wf_ex[2, :, :]
    # ===
    for ip in 1:NPR
        parti[ip].time = rc_kind(time)
        if (parti[ip].j < NJ && parti[ip].j > 0 && parti[ip].k < NK && parti[ip].k > 0)
            #ic, jc, kc, is the integer index of the particle relative to 
            #the grids center. Use these values for variables with the ghost points.
            #ifc, jfc, and kfc is the index relative to the coordinates of grid faces. 
            #Use these values for variables on faces.

            local ic = Int(parti[ip].i + 0.5e0)
            local jc = Int(parti[ip].j + 0.5e0)
            local kc = Int(parti[ip].k + 0.5e0)

            local ifc = Int(parti[ip].i)
            local jfc = Int(parti[ip].j)
            local kfc = Int(parti[ip].k)

            local dif = parti[ip].i - ifc
            local djf = parti[ip].j - jfc
            local dkf = parti[ip].k - kfc

            local dic = parti[ip].i - ic + 0.5e0
            local djc = parti[ip].j - jc + 0.5e0
            local dkc = parti[ip].k - kc + 0.5e0
            #calculate the zonal velocity
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].k + 0.5e0)::Ref{rc_kind}, swap1::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(kc)::Ref{rc_kind}, swap2::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(kc + 1e0)::Ref{rc_kind}, swap3::Ref{rc_kind})::Cvoid
            local dz = (swap1[] - swap2[]) / (swap3[] - swap2[])

            parti[ip].u = interp_trilinear(dif, djc, dz, uf_ex[ifc+1:ifc+2, jc+1:jc+2, kc+1:kc+2])

            #calculate the meridional velocity
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(parti[ip].k + 0.5e0)::Ref{rc_kind}, swap1::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(kc)::Ref{rc_kind}, swap2::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(kc + 1e0)::Ref{rc_kind}, swap3::Ref{rc_kind})::Cvoid
            dz = (swap1[] - swap2[]) / (swap3[] - swap2[])
            parti[ip].v = interp_trilinear(dic, djf, dz, vf_ex[ic+1:ic+2, jfc+1:jfc+2, kc+1:kc+2])

            #calculate the vertical velocity
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(parti[ip].k)::Ref{rc_kind}, swap1::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(kfc)::Ref{rc_kind}, swap2::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i)::Ref{rc_kind}, Ref(parti[ip].j)::Ref{rc_kind}, Ref(kfc + 1e0)::Ref{rc_kind}, swap3::Ref{rc_kind})::Cvoid
            dz = (swap1[] - swap2[]) / (swap3[] - swap2[])
            parti[ip].z = swap1[] * DL
            parti[ip].w = interp_trilinear(dic, djc, dz, wf_ex[ic+1:ic+2, jc+1:jc+2, kfc+1:kfc+2])

            #diagnose other properties
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].k + 0.5e0)::Ref{rc_kind}, swap1::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(kc)::Ref{rc_kind}, swap2::Ref{rc_kind})::Cvoid
            @ccall "PSOM_LIB.so".sigma2z_(Ref(parti[ip].i + 0.5e0)::Ref{rc_kind}, Ref(parti[ip].j + 0.5e0)::Ref{rc_kind}, Ref(kc + 1e0)::Ref{rc_kind}, swap3::Ref{rc_kind})::Cvoid
            dz = (swap1[] - swap2[]) / (swap3[] - swap2[])
            parti[ip].vor = interp_trilinear(dic, djc, dz, vor[ic:ic+1, jc:jc+1, kc:kc+1])
            parti[ip].rho = interp_trilinear(dic, djc, dz, rho[ic:ic+1, jc:jc+1, kc:kc+1])
            parti[ip].shear = interp_trilinear(dic, djc, dz, shear[ic:ic+1, jc:jc+1, kc:kc+1])
            parti[ip].strain = interp_trilinear(dic, djc, dz, strain[ic:ic+1, jc:jc+1, kc:kc+1])
        else
            parti[ip].u = 0e0
            parti[ip].v = 0e0
            parti[ip].w = 0e0
        end
    end
end

function parti_forward()
    for i in 1:NPR
        parti[i].i = parti[i].i + 0.5e0 * dtf * (3e0 * parti[i].u - parti[i].u0)
        if (parti[i].i > NI)
            parti[i].i = parti[i].i - rc_kind(NI)
        end
        if (parti[i].i < 0e0)
            parti[i].i = parti[i].i + rc_kind(NI)
        end

        if (parti[i].j > NJ - 1 && parti[i].v > 0e0)
            parti[i].j = parti[i].j + parti[i].v * @fortGet("dtf", rc_kind) / (1e0 + (parti[i].v * @fortGet("dtf", rc_kind)) / (rc_kind(NJ) - parti[i].j))
        elseif (parti[i].j < 1 && parti[i].v < 0e0)
            parti[i].j = parti[i].j + parti[i].v * @fortGet("dtf", rc_kind) / (1e0 - @fortGet("dtf", rc_kind) / parti[i].j)
        else
            parti[i].j = parti[i].j + 0.5e0 * @fortGet("dtf", rc_kind) * (3e0 * parti[i].v - parti[i].v0)
        end

        if (parti[i].k > NK - 1 && parti[i].w > 0)
            parti[i].k = parti[i].k + parti[i].w * @fortGet("dtf", rc_kind) / (1e0 + (parti[i].w * @fortGet("dtf", rc_kind)) / (rc_kind(NK) - parti[i].k))
        elseif (parti[i].k < 1 && parti[i].w < 0)
            parti[i].k = parti[i].k + parti[i].w * @fortGet("dtf", rc_kind) / (1e0 - @fortGet("dtf", rc_kind) / parti[i].k)
        else
            parti[i].k = parti[i].k + 0.5e0 * @fortGet("dtf", rc_kind) * (3e0 * parti[i].w - parti[i].w0)
        end

        if (parti[i].j < 0e0 || parti[i].j > NJ || parti[i].k > NK || parti[i].k < 0e0)
            println("particles coordinates are wrong, iPR=", i, "j,k", parti[i].j, parti[i].k)
        end

        parti[i].u0 = parti[i].u
        parti[i].v0 = parti[i].v
        parti[i].w0 = parti[i].w
    end
end

function get_parti_vel_ana()
    for ip in 1:NPR
        parti[ip].u = -1e0 * sin(PI * parti[ip].i / rc_kind(NI)) * cos(PI * parti[ip].j / rc_kind(NJ))
        parti[ip].v = cos(PI * parti[ip].i / rc_kind(NI)) * sin(PI * parti[ip].j / rc_kind(NJ))
    end
end

end