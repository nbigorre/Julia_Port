
module header

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
global const p = @fortGetArr("p", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const shear = @fortGetArr("shear", rc_kind, (NI+2,NJ+2,NK+2))
global const strain = @fortGetArr("strain", rc_kind, (NI+2,NJ+2,NK+2))
global const Jac = @fortGetArr("jac", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const JacInv = @fortGetArr("jacinv", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const vor = @fortGetArr("vor", rc_kind, (NI+2,NJ+2,NK+2))
global const rho = @fortGetArr("rho", rc_kind, (NI+2,NJ+2,NK+2))
global const rho_old = @fortGetArr("rho_old", rc_kind, (NI+2, NJ+2, NK+2))

global const freqby = @fortGetArr("freqby", rc_kind, (NI+2, NJ+2, NK+2))
global const freqbz = @fortGetArr("freqbz", rc_kind, (NI+2, NJ+2, NK+2))
global const freqN2 = @fortGetArr("freqn2", rc_kind, (NI+2, NJ+2, NK+2))
global const si = @fortGetArr("si", rc_kind, (NI+2, NJ+2, NK+2))
global const sj = @fortGetArr("sj", rc_kind, (NI+2, NJ+2, NK+2))
global const sk = @fortGetArr("sk", rc_kind, (NI+2, NJ+2, NK+2))
global const rp = @fortGetArr("rp", rc_kind, (NI+2, NJ+2, NK+2))

global const wt = @fortGetArr("wt", rc_kind, (NI+2, NJ+2, NK+1))
global const wzk = @fortGetArr("wzk", rc_kind, (NI+2, NJ+2, NK+1))
global const skfc = @fortGetArr("skfc", rc_kind, (NI+2, NJ+2, NK+1))

global const czf = @fortGetArr("czf", rc_kind, (NI, NJ, NK + 1))
global const Kz = @fortGetArr("kz", rc_kind, (NI, NJ, NK + 1))
global const wf = @fortGetArr("wf", rc_kind, (NI, NJ, NK + 1))

global const cx = @fortGetArr("cx", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const cy = @fortGetArr("cy", rc_kind, (NI + 2, NJ + 2, NK + 2))
global const cz = @fortGetArr("cz", rc_kind, (NI + 2, NJ + 2, NK + 2))

global const uvis = @fortGetArr("uvis", rc_kind, (NI, NJ, NK))
global const vvis = @fortGetArr("vvis", rc_kind, (NI, NJ, NK))
global const wvis = @fortGetArr("wvis", rc_kind, (NI, NJ, NK))

global const ufbcw = @fortGetArr("ufbcw", rc_kind, (NJ, NK))
global const ufbce = @fortGetArr("ufbce", rc_kind, (NJ, NK))

global const vfbcn = @fortGetArr("vfbcn", rc_kind, (NI, NK))
global const vfbcs = @fortGetArr("vfbcs", rc_kind, (NI, NK))

global const hyn = @fortGetArr("hyn", rc_kind, (NI, NJ + 1, NK))
global const sjfc = @fortGetArr("sjfc", rc_kind, (NI, NJ + 1, NK))
global const Jjfc = @fortGetArr("jjfc", rc_kind, (NI, NJ + 1, NK))
global const cyf = @fortGetArr("cyf", rc_kind, (NI, NJ + 1, NK))
global const gj3 = @fortGetArr("gj3", rc_kind, (NI, NJ + 1, NK))
global const gqj3 = @fortGetArr("gqj3", rc_kind, (NI, NJ + 1, NK))
global const vf = @fortGetArr("vf", rc_kind, (NI, NJ + 1, NK))

global const hxn = @fortGetArr("hxn", rc_kind, (NI + 1, NJ, NK))
global const sifc = @fortGetArr("sifc", rc_kind, (NI + 1, NJ, NK))
global const Jifc = @fortGetArr("jifc", rc_kind, (NI + 1, NJ, NK))
global const cxf = @fortGetArr("cxf", rc_kind, (NI + 1, NJ, NK))
global const gi3 = @fortGetArr("gi3", rc_kind, (NI + 1, NJ, NK))
global const gqi3 = @fortGetArr("gqi3", rc_kind, (NI + 1, NJ, NK))
global const uf = @fortGetArr("uf", rc_kind, (NI + 1, NJ, NK))

global const gradhn = @fortGetArr("gradhn", rc_kind, (NI + 2, NJ + 2, 2))

global const ux = @fortGetArr("ux", rc_kind, (NI + 2, NJ + 2))
global const uy = @fortGetArr("uy", rc_kind, (NI + 2, NJ + 2))
global const vx = @fortGetArr("vx", rc_kind, (NI + 2, NJ + 2))
global const vy = @fortGetArr("vy", rc_kind, (NI + 2, NJ + 2))
global const ffc = @fortGetArr("ffc", rc_kind, (NI + 2, NJ + 2))
global const bbc = @fortGetArr("bbc", rc_kind, (NI + 2, NJ + 2))
global const hdt = @fortGetArr("hdt", rc_kind, (NI + 2, NJ + 2))
global const g11 = @fortGetArr("g11", rc_kind, (NI + 2, NJ + 2))
global const g12 = @fortGetArr("g12", rc_kind, (NI + 2, NJ + 2))
global const g22 = @fortGetArr("g22", rc_kind, (NI + 2, NJ + 2))
global const D = @fortGetArr("d", rc_kind, (NI + 2, NJ + 2))



global const drpx = @fortGetArr("drpx", rc_kind, (NI, NJ, NK))
global const drpy = @fortGetArr("drpy", rc_kind, (NI, NJ, NK))

global const ffj = @fortGetArr("ffj", rc_kind, (NI, NJ+1))
global const bbj = @fortGetArr("bbj", rc_kind, (NI, NJ+1))

global const ffi = @fortGetArr("ffi", rc_kind, (NI+1, NJ))
global const bbi = @fortGetArr("bbi", rc_kind, (NI+1, NJ))

global const wfbcb = @fortGetArr("ufbce", rc_kind, (NI, NJ))


global const swr = @fortGetArr("swr", rc_kind, (NJ+2))
global const qloss = @fortGetArr("qloss", rc_kind, (NJ+2))

global const yc = @fortGetArr("yc", rc_kind, (NJ+2))

if (cppdefs.implicit)
    global const mat_A = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_B = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_C = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_D = zeros(rc_kind, (NI, NJ, NK, 5))
    global const mat_test = zeros(rc_kind, (NI, NJ, NK, 5))
end

if (cppdefs.allow_particle)
    global parti_file_num::Int = 0
    global NPR::Int = 0
end
global const J2d = @fortGetArr("j2d", rc_kind, (NI+2, NJ+2))

#=

S0=35.7e0::rc_kind
T0=15e0::rc_kind

OMEGA=7.272e-5::rc_kind

nsteps=0::Int
dtf=0e0::rc_kind
dtime_dim::rc_kind=0e0
dtime=0e0::rc_kind
time_nondim=0e0::rc_kind
time_seconds=0e0::rc_kind
kappah=0e0::rc_kind
kaphinv=0e0::rc_kind


dx=0e0::rc_kind
dy=0e0::rc_kind
lv_flat_bottom=true
total_depth=0e0::rc_kind
use_Shchepetkin=true

const z0= 0.2e0::rc_kind
const zm= 0.1e0::rc_kind
dztop=0e0::rc_kind
dztop_dim=0e0::rc_kind

# PREPROC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 83

fnhhy=0e0::rc_kind

phi0deg=0e0::rc_kind
fplane=0::Int

Kx=0e0::rc_kind
Ky=0e0::rc_kind

bottom_linear_drag=true
RR=0e0::rc_kind

lambda=0e0::rc_kind
DLinv=0e0::rc_kind
delta=0e0::rc_kind
delinv=0e0::rc_kind
qpr=0e0::rc_kind
beta=0e0::rc_kind
P1=0e0::rc_kind
UL=0e0::rc_kind
WL=0e0::rc_kind
TL=0e0::rc_kind
HL=0e0::rc_kind
HDL=0e0::rc_kind

user1=0e0::rc_kind
user2=0e0::rc_kind
user3=0e0::rc_kind
user4=0e0::rc_kind
user5=0e0::rc_kind
user6=0e0::rc_kind
user7=0e0::rc_kind
user8=0e0::rc_kind

dirout=fill(UInt8(32), 151);
out1d_int=0::Int
out2d_int=0::Int
out3d_int=0::Int

const verbose=true

# PREPROC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 181

ivb=0::Int
ivs=0::Int
ivf=0::Int

pfac=0e0::rc_kind

stress_top=zeros(rc_kind, NI, NJ)
stress_top_x=zeros(rc_kind, NI, NJ)
stress_top_y=zeros(rc_kind, NI, NJ)

yfront=0e0::rc_kind
dyfront=0e0::rc_kind
sclwidth=0e0::rc_kind
tightness=0e0::rc_kind
mldepth=0e0::rc_kind
sigrelease=zeros(rc_kind, ntr)
drho=0e0::rc_kind
conv=zeros(Int, NI+2, NJ+2, NK+2)
con100=zeros(Int, NI+2, NJ+2, NK+2)

mldn2=zeros(rc_kind,3)
mldn2init=zeros(rc_kind,3)
zbtop=zeros(rc_kind,3)
zbbot=zeros(rc_kind,3)
zbbotinit=zeros(rc_kind,3)
zbtopinit=zeros(rc_kind,3)
fricbot=zeros(rc_kind,3)
frictop=zeros(rc_kind,3)
diabot=zeros(rc_kind,3)
diatop=zeros(rc_kind,3)
advecpv=zeros(rc_kind,3)
friction=zeros(rc_kind,3)
diabatic=zeros(rc_kind,3)

sbackgrnd=0e0::rc_kind
r_sponge=zeros(rc_kind, NJ+2)
r_T=zeros(rc_kind, NJ+2, NK+1)
rho_refS=zeros(rc_kind, NK+2)
rho_refN=zeros(rc_kind, NK+2)
rho_refB=zeros(rc_kind,NJ+2)

hmean=0e0::rc_kind
it=0::Int
dum=0e0::rc_kind

iddatfile=0::Int
idigit=0::Int
idudx=0::Int
idudy=0::Int
idudz=0::Int
idvbysection=0::Int
idvby=0::Int
idvbz=0::Int
idvcon=0::Int
idvcy=0::Int
idvbx=0::Int
idvcz=0::Int
idvc=0::Int
idvdivfeddy=0::Int
idvdivfreyn=0::Int
idvdx=0::Int
idvdz=0::Int
idvd=0::Int
idvfb=0::Int
idvh=0::Int
idvn2bar=0::Int
idvn2=0::Int
idvnsq100m=0::Int
idvnsq30m=0::Int
idvpe=0::Int
idvpsiv=0::Int
idvpsiw=0::Int
idvpv=0::Int
idvp=0::Int
idvrhbar=0::Int
idvrho=0::Int
idvrnk=0::Int
idvstrain=0::Int
idvstress=0::Int
idvstr=0::Int
idvs=0::Int
idvtbar=0::Int
idvtemp=0::Int
idvtim=0::Int
idvtr=0::Int
idvt=0::Int
idvu=0::Int
idvvb=0::Int
idvvc=0::Int
idvvor=0::Int
idvv=0::Int
idvwb=0::Int
idvwc=0::Int
idvwpv=0::Int
idvw=0::Int
idvy=0::Int
idvzsave=0::Int
idvz=0::Int
idwdx=0::Int
idwdy=0::Int
idwdz=0::Int
iimday=0::Int
ipos=0::Int
frame_int=0::Int
ngraph2d=0::Int

pickup_step=0::Int
pickup_int=0::Int

lv_test_output_bin=false

initial_step=0::Int

Tr=zeros(rc_kind, ntr, NI+2, NJ+2, NK+2, 2)
wtr=zeros(rc_kind, NI+2, NJ+2, NK+2)

u=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
v=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
w=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
s=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
T=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)

u_bar=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
v_bar=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
w_bar=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
s_bar=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
T_bar=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)

u_pri=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
v_pri=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
w_pri=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
s_pri=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)
T_pri=zeros(rc_kind, NI+2, NJ+2, NK+2, 2)

gqk=zeros(rc_kind, NI+2, NJ+2, NK+1, 3)

gi = zeros(rc_kind, NI+1, NJ, NK, 2)
gqi = zeros(rc_kind, NI+1, NJ, NK, 2)

gj = zeros(rc_kind, NI, NJ+1, NK, 2)
gqj = zeros(rc_kind, NI, NJ+1, NK, 2)

=#

@exportAll()


end