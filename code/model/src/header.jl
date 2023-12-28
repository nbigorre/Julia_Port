
module header

using OffsetArrays
using ExportAll
using ..fortVar
using ..cppdefs

global const rc_kind = Float64

# SIZE
include("../inc/size.jl")


###############     CONSTANTS     ###############

global const EPS::rc_kind = 0.1e0
global const FPAR::rc_kind=1e-4
global const AL::rc_kind=1.e7

global const LEN::rc_kind=1e5
global const DL::rc_kind=1e3
global const R0::rc_kind=1027e0
global const PI::rc_kind=3.14159265358979323846e0
global const periodicew=true
global const gpr::rc_kind = 0.981e0
global const apr::rc_kind=0.6371e0
global const rect=true

###############      ARRAYS       ###############

global const conv = @fortGetArr("conv", Int, (NI + 2, NJ + 2, NK + 2))
global const con100 = @fortGetArr("con100", Int, (NI + 2, NJ + 2, NK + 2))

global const ru_Sche = @fortGetArr("ru_sche", rc_kind, (NI+2,NJ+2,NK))
global const rv_Sche = @fortGetArr("rv_sche", rc_kind, (NI+2,NJ+2,NK))
global const ru2_Sche = @fortGetArr("ru2_sche", rc_kind, (NI+2,NJ+2,NK))
global const rv2_Sche = @fortGetArr("rv2_sche", rc_kind, (NI+2,NJ+2,NK))
global const ru3_Sche = @fortGetArr("ru3_sche", rc_kind, (NI+2,NJ+2,NK))
global const rv3_Sche = @fortGetArr("rv3_sche", rc_kind, (NI+2,NJ+2,NK))
global const ru4_Sche = @fortGetArr("ru4_sche", rc_kind, (NI+2,NJ+2,NK))
global const rv4_Sche = @fortGetArr("rv4_sche", rc_kind, (NI+2,NJ+2,NK))
global const grpifc = @fortGetArr("grpifc", rc_kind, (NI+1, NJ, NK))
global const grpjfc = @fortGetArr("grpjfc", rc_kind, (NI, NJ+1, NK))
global const T_ref = @fortGetArr("t_ref", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const Tr = @fortGetArr("tr", rc_kind, (ntr, NI + 2, NJ + 2, NK + 2, 2))
global const h = @fortGetArr("h", rc_kind, (NI+2, NJ+2))
global const oldh = @fortGetArr("oldh", rc_kind, (NI+2, NJ+2))
global const r_sponge = @fortGetArr("r_sponge", rc_kind, (NJ+2))
global const stress_top = @fortGetArr("stress_top", rc_kind, (NI, NJ))
global const stress_top_x = @fortGetArr("stress_top_x", rc_kind, (NI, NJ))
global const stress_top_y = @fortGetArr("stress_top_y", rc_kind, (NI, NJ))

global const u = @fortGetArr("u", rc_kind, (NI + 2, NJ + 2, NK + 2, 2))
global const v = @fortGetArr("v", rc_kind, (NI + 2, NJ + 2, NK + 2, 2))
global const w = @fortGetArr("w", rc_kind, (NI + 2, NJ + 2, NK + 2, 2))
global const s = @fortGetArr("s", rc_kind, (NI + 2, NJ + 2, NK + 2, 2))
global const T = @fortGetArr("t", rc_kind, (NI + 2, NJ + 2, NK + 2, 2))

global const gqk = @fortGetArr("gqk", rc_kind, (NI+2, NJ+2, NK+1, 3))

global const gi = @fortGetArr("gi", rc_kind, (NI + 1, NJ, NK, 2))
global const gqi = @fortGetArr("gqi", rc_kind, (NI + 1, NJ, NK, 2))

global const gj = @fortGetArr("gj", rc_kind, (NI, NJ + 1, NK, 2))
global const gqj = @fortGetArr("gqj", rc_kind, (NI, NJ + 1, NK, 2))

global const zf = @fortGetArr("zf", rc_kind, (NI + 2, NJ + 2, NK + 3))

global const zc = @fortGetArr("zc", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const wx = @fortGetArr("wx", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const wy = @fortGetArr("wy", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const wz = @fortGetArr("wz", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const p = OffsetArray(@fortGetArr("p", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const shear = OffsetArray(@fortGetArr("shear", rc_kind, (NI+2,NJ+2,NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const strain = OffsetArray(@fortGetArr("strain", rc_kind, (NI+2,NJ+2,NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const Jac = OffsetArray(@fortGetArr("jac", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const JacInv = OffsetArray(@fortGetArr("jacinv", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const vor = OffsetArray(@fortGetArr("vor", rc_kind, (NI+2,NJ+2,NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rho = OffsetArray(@fortGetArr("rho", rc_kind, (NI+2,NJ+2,NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rho_old = OffsetArray(@fortGetArr("rho_old", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)

global const freqby = OffsetArray(@fortGetArr("freqby", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const freqbz = OffsetArray(@fortGetArr("freqbz", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const freqN2 = OffsetArray(@fortGetArr("freqn2", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const si = OffsetArray(@fortGetArr("si", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const sj = OffsetArray(@fortGetArr("sj", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const sk = OffsetArray(@fortGetArr("sk", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const rp = OffsetArray(@fortGetArr("rp", rc_kind, (NI+2, NJ+2, NK+2)), 0:NI+1, 0:NJ+1, 0:NK+1)

global const wt = OffsetArray(@fortGetArr("wt", rc_kind, (NI+2, NJ+2, NK+1)), 0:NI+1, 0:NJ+1, 0:NK)
global const wzk = OffsetArray(@fortGetArr("wzk", rc_kind, (NI+2, NJ+2, NK+1)), 0:NI+1, 0:NJ+1, 0:NK)
global const skfc = OffsetArray(@fortGetArr("skfc", rc_kind, (NI+2, NJ+2, NK+1)), 0:NI+1, 0:NJ+1, 0:NK)

global const czf = OffsetArray(@fortGetArr("czf", rc_kind, (NI, NJ, NK + 1)), 1:NI, 1:NJ, 0:NK)
global const Kz = OffsetArray(@fortGetArr("kz", rc_kind, (NI, NJ, NK + 1)), 1:NI, 1:NJ, 0:NK)
global const wf = OffsetArray(@fortGetArr("wf", rc_kind, (NI, NJ, NK + 1)), 1:NI, 1:NJ, 0:NK)

global const cx = OffsetArray(@fortGetArr("cx", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const cy = OffsetArray(@fortGetArr("cy", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)
global const cz = OffsetArray(@fortGetArr("cz", rc_kind, (NI + 2, NJ + 2, NK + 2)), 0:NI+1, 0:NJ+1, 0:NK+1)

global const uvis = @fortGetArr("uvis", rc_kind, (NI, NJ, NK))
global const vvis = @fortGetArr("vvis", rc_kind, (NI, NJ, NK))
global const wvis = @fortGetArr("wvis", rc_kind, (NI, NJ, NK))

global const ufbcw = @fortGetArr("ufbcw", rc_kind, (NJ, NK))
global const ufbce = @fortGetArr("ufbce", rc_kind, (NJ, NK))

global const vfbcn = @fortGetArr("vfbcn", rc_kind, (NI, NK))
global const vfbcs = @fortGetArr("vfbcs", rc_kind, (NI, NK))

global const hyn = OffsetArray(@fortGetArr("hyn", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const sjfc = OffsetArray(@fortGetArr("sjfc", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const Jjfc = OffsetArray(@fortGetArr("jjfc", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const cyf = OffsetArray(@fortGetArr("cyf", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const gj3 = OffsetArray(@fortGetArr("gj3", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const gqj3 = OffsetArray(@fortGetArr("gqj3", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)
global const vf = OffsetArray(@fortGetArr("vf", rc_kind, (NI, NJ + 1, NK)), 1:NI, 0:NJ, 1:NK)

global const hxn = OffsetArray(@fortGetArr("hxn", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const sifc = OffsetArray(@fortGetArr("sifc", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const Jifc = OffsetArray(@fortGetArr("jifc", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const cxf = OffsetArray(@fortGetArr("cxf", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const gi3 = OffsetArray(@fortGetArr("gi3", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const gqi3 = OffsetArray(@fortGetArr("gqi3", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)
global const uf = OffsetArray(@fortGetArr("uf", rc_kind, (NI + 1, NJ, NK)), 0:NI, 1:NJ, 1:NK)

global const gradhn = OffsetArray(@fortGetArr("gradhn", rc_kind, (NI + 2, NJ + 2, 2)), 0:NI+1, 0:NJ+1, 1:2)

global const ux = OffsetArray(@fortGetArr("ux", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const uy = OffsetArray(@fortGetArr("uy", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const vx = OffsetArray(@fortGetArr("vx", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const vy = OffsetArray(@fortGetArr("vy", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const ffc = OffsetArray(@fortGetArr("ffc", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const bbc = OffsetArray(@fortGetArr("bbc", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const hdt = OffsetArray(@fortGetArr("hdt", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const g11 = OffsetArray(@fortGetArr("g11", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const g12 = OffsetArray(@fortGetArr("g12", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const g22 = OffsetArray(@fortGetArr("g22", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const D = OffsetArray(@fortGetArr("d", rc_kind, (NI + 2, NJ + 2)),0:NI+1, 0:NJ+1)
global const J2d = OffsetArray(@fortGetArr("j2d", rc_kind, (NI+2, NJ+2)),0:NI+1, 0:NJ+1)



global const drpx = @fortGetArr("drpx", rc_kind, (NI, NJ, NK))
global const drpy = @fortGetArr("drpy", rc_kind, (NI, NJ, NK))

global const ffj = OffsetArray(@fortGetArr("ffj", rc_kind, (NI, NJ+1)), 1:NI, 0:NJ)
global const bbj = OffsetArray(@fortGetArr("bbj", rc_kind, (NI, NJ+1)), 1:NI, 0:NJ)

global const ffi = OffsetArray(@fortGetArr("ffi", rc_kind, (NI+1, NJ)), 0:NI, 1:NJ)
global const bbi = OffsetArray(@fortGetArr("bbi", rc_kind, (NI+1, NJ)), 0:NI, 1:NJ)

global const wfbcb = @fortGetArr("ufbce", rc_kind, (NI, NJ))

global const yc = OffsetArray(@fortGetArr("yc", rc_kind, (NJ+2)), 0:NJ+1)

global const xc = OffsetArray(@fortGetArr("xc", rc_kind, (NI+2)), 0:NI+1)


global const swr = OffsetArray(@fortGetArr("swr", rc_kind, (NJ+2)), 0:NJ+1)
global const qloss = OffsetArray(@fortGetArr("qloss", rc_kind, (NJ+2)), 0:NJ+1)


@static if (cppdefs.implicit)
    global const mat_A = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_B = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_C = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_D = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_test = zeros(rc_kind, (NI, NJ, NK, 5))
end

@static if (cppdefs.allow_particle)
    global parti_file_num::Int = 0
    global NPR::Int = 0
end


@exportAll()


end