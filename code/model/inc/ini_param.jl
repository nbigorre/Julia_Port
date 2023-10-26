#import Pkg
#Pkg.add("YAML")

using YAML
using .header
using .fortVar

data = YAML.load_file(ARGS[1])

@fortSet("nsteps", Int(data["PARAM"]["nsteps"]), Int)
@fortSet("dtime_dim", rc_kind(data["PARAM"]["dtime_dim"]), rc_kind)
@fortSet("fnhhy", rc_kind(data["PARAM"]["fnhhy"]), rc_kind)
@fortSet("fplane", Int(data["PARAM"]["fplane"]), Int)
@fortSet("phi0deg", rc_kind(data["PARAM"]["phi0deg"]), rc_kind)
@fortSet("dx", rc_kind(data["PARAM"]["dx"]), rc_kind)
@fortSet("dy", rc_kind(data["PARAM"]["dy"]), rc_kind)
@fortSet("dztop_dim", rc_kind(data["PARAM"]["dztop_dim"]), rc_kind)
@fortSet("lv_flat_bottom", Int(data["PARAM"]["lv_flat_bottom"]), Int)
@fortSet("total_depth", rc_kind(data["PARAM"]["total_depth"]), rc_kind)
@fortSet("use_shchepetkin", Int(data["PARAM"]["use_Shchepetkin"]), Int)
@fortSet("kx", rc_kind(data["PARAM"]["Kx"]), rc_kind)
@fortSet("ky", rc_kind(data["PARAM"]["Ky"]), rc_kind)
@fortSet("rr", rc_kind(data["PARAM"]["RR"]), rc_kind)
@fortSet("pickup_int", Int(data["PARAM"]["pickup_int"]), Int)
@fortSet("pickup_step", Int(data["PARAM"]["pickup_step"]), Int)
@fortSet("out1d_int", Int(data["PARAM"]["out1d_int"]), Int)
@fortSet("out2d_int", Int(data["PARAM"]["out2d_int"]), Int)
@fortSet("out3d_int", Int(data["PARAM"]["out3d_int"]), Int)
fortSetStr("dirout", data["PARAM"]["dirout"], Ptr{UInt8}, 151)

@fortSet("user1", rc_kind(data["USER"]["user1"]), rc_kind)
@fortSet("user2", rc_kind(data["USER"]["user2"]), rc_kind)
@fortSet("user3", rc_kind(data["USER"]["user3"]), rc_kind)
@fortSet("user4", rc_kind(data["USER"]["user4"]), rc_kind)
@fortSet("user5", rc_kind(data["USER"]["user5"]), rc_kind)
@fortSet("user6", rc_kind(data["USER"]["user6"]), rc_kind)
@fortSet("user7", rc_kind(data["USER"]["user7"]), rc_kind)
@fortSet("user8", rc_kind(data["USER"]["user8"]), rc_kind)




@fortSet("kappah", 0.65e0, rc_kind)
@fortSet("kaphinv", 1.e0/fortGet("kappah", rc_kind), rc_kind)

@fortSet("lambda", rc_kind(LEN / AL), rc_kind)
@fortSet("delta", rc_kind(DL/ LEN), rc_kind)
@fortSet("delinv", rc_kind(LEN/DL), rc_kind)

@fortSet("dztop", rc_kind(fortGet("dztop_dim",rc_kind)/DL), rc_kind)

@fortSet("qpr", rc_kind(fortGet("fnhhy",rc_kind)*fortGet("delta",rc_kind)), rc_kind)
@fortSet("beta", rc_kind(1e0/(EPS*EPS*fortGet("delta", rc_kind))) , rc_kind)
@fortSet("ul", rc_kind(FPAR *LEN*EPS), rc_kind)
@fortSet("wl", rc_kind(EPS*fortGet("delta", rc_kind)*fortGet("ul", rc_kind)), rc_kind)
@fortSet("p1", rc_kind(R0*fortGet("ul", rc_kind)*fortGet("ul", rc_kind)*(EPS^(-1))), rc_kind)
@fortSet("hl", rc_kind(fortGet("p1", rc_kind)/(R0*10e0)), rc_kind)
@fortSet("hdl", rc_kind(fortGet("hl", rc_kind)/DL), rc_kind)
@fortSet("tl", rc_kind(LEN/fortGet("ul", rc_kind)), rc_kind)

@fortSet("dtf", rc_kind(fortGet("dtime_dim", rc_kind) / fortGet("tl", rc_kind)), rc_kind)