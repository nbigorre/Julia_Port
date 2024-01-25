function save3d(Nx,Ny,Nz,var,filename)
   #use header, only : rc_kind
   #integer :: Nx,Ny,Nz
   #character(len=*) :: filename
   #REAL(kind=rc_kind) ::  var(Nx,Ny,Nz)
   open(filename, "w") do file
      write(file, var)
   end
end

function save2d(Nx,Ny,var,filename)
   #use header, only : rc_kind
   #!Nx,Ny are not necessary in zonal and meridional directions.
   #integer :: Nx,Ny
   #character(len=*) :: filename
   #REAL(kind=rc_kind) ::  var(Nx,Ny)
   open(filename, "w") do file
      write(file, var)
   end
end 

function w_pickup(filepath)
   #use header, only : h,u,v,w,uf,vf,wf,T,S,ntr,Tr,p,gradhn,hxn,hyn,rc_kind
   #character(len=*) :: filename

      local dntr = NcDim("ntr", ntr)
      local dni2 = NcDim("NI2", NI+2)
      local dnj2 = NcDim("NJ2", NJ+2)
      local dnk2 = NcDim("NK2", NK+2)
      local dn = NcDim("n", 2)
      local dni1 = NcDim("NI1", NI+1)
      local dnj1 = NcDim("NJ1", NJ+1)
      local dnk1 = NcDim("NK1", NK+1)
      local dni = NcDim("NI", NI)
      local dnj = NcDim("NJ", NJ)
      local dnk = NcDim("NK", NK)
   
   local varlist = [
      NcVar("h", [dni2, dnj2], t=Float64),
      NcVar("u", [dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("v", [dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("w", [dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("uf", [dni1, dnj, dnk], t=Float64),
      NcVar("vf", [dni, dnj1, dnk], t=Float64),
      NcVar("wf", [dni, dnj, dnk1], t=Float64),
      NcVar("T", [dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("s", [dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("Tr", [dntr, dni2, dnj2, dnk2, dn], t=Float64),
      NcVar("p", [dni2, dnj2, dnk2], t=Float64),
      NcVar("gradhn", [dni2, dnj2, dn], t=Float64),
      NcVar("hxn", [dni1, dnj, dnk], t=Float64),
      NcVar("hyn", [dni, dnj1, dnk], t=Float64)
   ]


   NetCDF.create(filepath, varlist, mode=NC_CLASSIC_MODEL)

   local ncfile = NetCDF.open(filepath, mode=NC_WRITE)

   NetCDF.putvar(ncfile["h"], OffsetArrays.no_offset_view(h))
   NetCDF.putvar(ncfile["u"], OffsetArrays.no_offset_view(u))
   NetCDF.putvar(ncfile["v"], OffsetArrays.no_offset_view(v))
   NetCDF.putvar(ncfile["w"], OffsetArrays.no_offset_view(w))
   NetCDF.putvar(ncfile["uf"], OffsetArrays.no_offset_view(uf))
   NetCDF.putvar(ncfile["vf"], OffsetArrays.no_offset_view(vf))
   NetCDF.putvar(ncfile["wf"], OffsetArrays.no_offset_view(wf))
   NetCDF.putvar(ncfile["T"], OffsetArrays.no_offset_view(T))
   NetCDF.putvar(ncfile["s"], OffsetArrays.no_offset_view(s))
   NetCDF.putvar(ncfile["Tr"], OffsetArrays.no_offset_view(Tr))
   NetCDF.putvar(ncfile["p"], OffsetArrays.no_offset_view(p))
   NetCDF.putvar(ncfile["gradhn"], OffsetArrays.no_offset_view(gradhn))
   NetCDF.putvar(ncfile["hxn"], OffsetArrays.no_offset_view(hxn))
   NetCDF.putvar(ncfile["hyn"], OffsetArrays.no_offset_view(hyn))

   NetCDF.sync(ncfile)

end

function r_pickup(filepath)
   #use header, only : h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn,rc_kind,dirout
   #integer :: step
   #character(len=10) :: stepchar

   local ncfile = NetCDF.open(filepath, mode=NC_NOWRITE)
   @views OffsetArrays.no_offset_view(h) .= (NetCDF.readvar(ncfile["h"]))
   @views OffsetArrays.no_offset_view(u) .= NetCDF.readvar(ncfile["u"])
   @views OffsetArrays.no_offset_view(v) .= NetCDF.readvar(ncfile["v"])
   @views OffsetArrays.no_offset_view(w) .= NetCDF.readvar(ncfile["w"])
   @views OffsetArrays.no_offset_view(uf) .= NetCDF.readvar(ncfile["uf"])
   @views OffsetArrays.no_offset_view(vf) .= NetCDF.readvar(ncfile["vf"])
   @views OffsetArrays.no_offset_view(wf) .= NetCDF.readvar(ncfile["wf"])
   @views OffsetArrays.no_offset_view(T) .= NetCDF.readvar(ncfile["T"])
   @views OffsetArrays.no_offset_view(s) .= NetCDF.readvar(ncfile["s"])
   @views OffsetArrays.no_offset_view(Tr) .= NetCDF.readvar(ncfile["Tr"])
   @views OffsetArrays.no_offset_view(p) .= NetCDF.readvar(ncfile["p"])
   @views OffsetArrays.no_offset_view(gradhn) .= NetCDF.readvar(ncfile["gradhn"])
   @views OffsetArrays.no_offset_view(hxn) .= NetCDF.readvar(ncfile["hxn"])
   @views OffsetArrays.no_offset_view(hyn) .= NetCDF.readvar(ncfile["hyn"])

   println("# pickup at step ", step)
end
