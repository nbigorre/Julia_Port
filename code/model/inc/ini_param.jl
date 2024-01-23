#import Pkg
#Pkg.add("YAML")

using .fortVar

data = YAML.load_file(ARGS[1])

nsteps = Int(data["PARAM"]["nsteps"])
@fortSet("nsteps", nsteps, Int)

dtime_dim = rc_kind(data["PARAM"]["dtime_dim"])
@fortSet("dtime_dim", dtime_dim, rc_kind)

fnhhy = rc_kind(data["PARAM"]["fnhhy"])
@fortSet("fnhhy", fnhhy, rc_kind)

fplane = Int(data["PARAM"]["fplane"])
@fortSet("fplane", fplane, Int)

phi0deg = rc_kind(data["PARAM"]["phi0deg"])
@fortSet("phi0deg", phi0deg, rc_kind)

dx = rc_kind(data["PARAM"]["dx"])
@fortSet("dx", dx, rc_kind)

dy = rc_kind(data["PARAM"]["dy"])
@fortSet("dy", dy, rc_kind)

dztop_dim = rc_kind(data["PARAM"]["dztop_dim"])
@fortSet("dztop_dim", dztop_dim, rc_kind)

lv_flat_bottom = Int(data["PARAM"]["lv_flat_bottom"])
@fortSet("lv_flat_bottom", lv_flat_bottom, Int)

total_depth = rc_kind(data["PARAM"]["total_depth"])
@fortSet("total_depth", total_depth, rc_kind)

use_shchepetkin = Int(data["PARAM"]["use_Shchepetkin"])
@fortSet("use_shchepetkin", use_shchepetkin, Int)

Kx = rc_kind(data["PARAM"]["Kx"])
@fortSet("kx", Kx, rc_kind)

Ky = rc_kind(data["PARAM"]["Ky"])
@fortSet("ky", Ky, rc_kind)

RR = rc_kind(data["PARAM"]["RR"])
@fortSet("rr", RR, rc_kind)

pickup_int = Int(data["PARAM"]["pickup_int"])
@fortSet("pickup_int", pickup_int, Int)

pickup_step = Int(data["PARAM"]["pickup_step"])
@fortSet("pickup_step", pickup_step, Int)

out1d_int = Int(data["PARAM"]["out1d_int"])
@fortSet("out1d_int", out1d_int, Int)

out2d_int = Int(data["PARAM"]["out2d_int"])
@fortSet("out2d_int", out2d_int, Int)

out3d_int = Int(data["PARAM"]["out3d_int"])
@fortSet("out3d_int", out3d_int, Int)

dirout = string(data["PARAM"]["dirout"])
fortSetStr("dirout", data["PARAM"]["dirout"], Ptr{UInt8}, 151)

user1 = rc_kind(data["USER"]["user1"])
@fortSet("user1", user1, rc_kind)
user2 = rc_kind(data["USER"]["user2"])
@fortSet("user2", user2, rc_kind)
user3 = rc_kind(data["USER"]["user3"])
@fortSet("user3", user3, rc_kind)
user4 = rc_kind(data["USER"]["user4"])
@fortSet("user4", user4, rc_kind)
user5 = rc_kind(data["USER"]["user5"])
@fortSet("user5", user5, rc_kind)
user6 = rc_kind(data["USER"]["user6"])
@fortSet("user6", user6, rc_kind)
user7 = rc_kind(data["USER"]["user7"])
@fortSet("user7", user7, rc_kind)
user8 = rc_kind(data["USER"]["user8"])
@fortSet("user8", user8, rc_kind)




kappah = 0.65e0
@fortSet("kappah", kappah, rc_kind)
@fortSet("kaphinv", 1.e0/fortGet("kappah", rc_kind), rc_kind)

lambda = rc_kind(LEN / AL)
@fortSet("lambda", lambda, rc_kind)
delta = rc_kind(DL/ LEN)
@fortSet("delta", delta, rc_kind)
delinv = rc_kind(LEN/DL)
@fortSet("delinv", delinv, rc_kind)

dztop = rc_kind(dztop_dim/DL)
@fortSet("dztop", dztop, rc_kind)

qpr = rc_kind(fnhhy*delta)
@fortSet("qpr", qpr, rc_kind)

beta = rc_kind(1e0/(EPS*EPS*delta))
@fortSet("beta", beta , rc_kind)

UL = rc_kind(FPAR *LEN*EPS)
@fortSet("ul", UL, rc_kind)

WL = rc_kind(EPS*delta*UL)
@fortSet("wl", WL, rc_kind)

P1 = rc_kind(R0*UL*UL*(EPS^(-1)))
@fortSet("p1", P1, rc_kind)

HL = rc_kind(P1/(R0*10e0))
@fortSet("hl", HL, rc_kind)

HDL = rc_kind(HL/DL)
@fortSet("hdl", HDL, rc_kind)

TL = rc_kind(LEN/UL)
@fortSet("tl", TL, rc_kind)

dtf = rc_kind(dtime_dim / TL)
@fortSet("dtf", dtf, rc_kind)