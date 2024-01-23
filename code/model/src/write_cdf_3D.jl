function write_cdf_3D(stepl,n)

#----------------------------------------------------------             
# 3D output routine.
# It writes both center and face values to be written in .cdf files.
#----------------------------------------------------------             

# USE header, only : NI,NJ,NK,xc,yc,zc,p,h,consump,T,s,rho,Tr,u,v,w,vor,pv,uf,vf,wf,Kz,conv,con100,nconsume,dirout,rc_kind,ntr
#include "netcdf.inc"                                                                        

#  integer i,j,k,n,stepl,nstp,it 
# 
#  character(LEN=150) outname_data, outname_face
#
#  REAL(kind=rc_kind) ::  Trwrite(0:NI+1,0:NJ+1,0:NK+1,ntr)
#  REAL(kind=rc_kind) :: z(0:NI+1,0:NJ+1,0:NK+1) 
#  REAL(kind=rc_kind) :: smax,Tmax,umax,vmax,wmax,smin,Tmin,umin,vmin,wmin 
#  
#  integer start(3), count(3), start2d(2), count2d(2),count2dtr(3),count4(4),start4(4),start2dtr(3),count4consump(4)
#  integer countuf(3), countvf(3), countwf(3)
#  
#  integer dims(3), dims2d(2),dims2dtr(3),dims4(4),dimsconsump(4)                             
#  integer dimuf(3), dimvf(3), dimwf(3)
#
#  integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
#  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
#  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
#  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
#  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idvz3,idwdx,idwdy,idwdz,iimday,ipos
#  integer :: idvx,idvcon100,idFaceFile,idvuf, idvvf, idvwf, idvKf
#  REAL(kind=rc_kind) ::  rcode
#                                                                        
#  DATA start /1, 1, 1/ 
#  DATA start4 /1, 1, 1, 1/ 
#  DATA start2d /1, 1/ 
#  DATA start2dtr /1, 1, 1/ 
local start = [1, 1, 1]
local start4 = [1, 1, 1, 1]
local start2d = [1, 1]
local start2dtr = [1, 1, 1]
local z = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1)
local Trwrite = OffsetArrays.zeros(rc_kind, 0:NI+1, 0:NJ+1, 0:NK+1, ntr)
                                                                        
  for k in 0:NK+1 
    for j in 0:NJ+1 
      for i in 0:NI+1 
        z[i,j,k]= zc[i,j,k] 
        for it in 1:ntr 
          Trwrite[i,j,k,it]= Tr[it,i,j,k,n] 
        end
      end
    end
  end

                                                                        
  # For the sake of better plots 
  for j in 0:NJ+1 
    for i in 0:NI+1 
      z[i,j,NK+1]= 0e0
      z[i,j,0]= 0.5*(z[i,j,0]+z[i,j,1])
      s[i,j,NK+1,n]= s[i,j,NK,n]
       s[i,j,0,n]= s[i,j,1,n] 
      T[i,j,NK+1,n]= T[i,j,NK,n]
       T[i,j,0,n]= T[i,j,1,n] 
      rho[i,j,NK+1]= rho[i,j,NK]
       rho[i,j,0]= rho[i,j,1] 
      u[i,j,NK+1,n]= u[i,j,NK,n]
       u[i,j,0,n]= u[i,j,1,n] 
      v[i,j,NK+1,n]= v[i,j,NK,n]
       v[i,j,0,n]= v[i,j,1,n] 
      w[i,j,NK+1,n]= w[i,j,NK,n]
       w[i,j,0,n]= w[i,j,1,n] 
      vor[i,j,NK+1]= vor[i,j,NK]
       vor[i,j,0]= vor[i,j,1] 
      pv[i,j,NK+1] =  pv[i,j,NK]
       pv[i,j,0] = pv[i,j,1]
      for it in 1:ntr 
        Trwrite[i,j,NK+1,it]= Trwrite[i,j,NK,it] 
        Trwrite[i,j,0,it]= Trwrite[i,j,1,it] 
      end
    end
  end
                              

  # Output file names    
  local outname_data = string("full_", lpad(stepl, 5, "0"), ".cdf") # Cell centers
  local outname_face = string("face_", lpad(stepl, 5, "0"), ".cdf") # Cell faces                                        


  # ---- END OF THE INITIALIZATION PART -----------------------------------------------
  

  # --- START WRITING

  #--------------------------------------- 
  # Write values at the cell centers:


  # 3D dimensions for most variables
  count = [NI+2, NJ+2, NK+2]
  dims = [
    NcDim("x", NI+2),
    NcDim("y", NJ+2),
    NcDim("sigma", NK+2)
  ]

  # 2D dimensions (for h) 
  local count2d = [NI+2, NJ+2]
  local dims2d = [dims[1], dims[2]]

  # 3D dimensions (for tracer)   
  local count4 = [NI+2, NJ+2, NK+2, ntr]
  local dims4 = [dims..., NcDim("ntr", ntr)]                                       
  # 2D dimensions for tracer (currently not used)
  local count2dtr = [NI+2, NJ+2, ntr]
  local dims2dtr = [ dims[1], dims[2], dims4[4]]

  # 3D dimensions (for consump)
  local count4consump = [NI+2, NJ+2, NK+2, nconsume]
  local dimsconsump = [dims..., NcDim("ntrcon", nconsume)]

  local varlist = [
    NcVar("xc", dims[1], t=Float64),
    NcVar("yc", dims[2], t=Float64),
    NcVar("zc", dims, t=Float64),
    NcVar("h", dims2d, t=Float64),
    NcVar("consump", dimsconsump, t=Float64),
    NcVar("tr", dims4, t=Float64),
    NcVar("s", dims, t=Float64),
    NcVar("temp", dims, t=Float64),
    NcVar("rho", dims, t=Float64),
    NcVar("p", dims, t=Float64),
    NcVar("u", dims, t=Float64),
    NcVar("v", dims, t=Float64),
    NcVar("w", dims, t=Float64),
    NcVar("vor", dims, t=Float64),
    NcVar("pv", dims, t=Float64),
    NcVar("conv", dims, t=Int16),
    NcVar("con100", dims, t=Int16)
  ]

  NetCDF.create(string(dirout, "/", outname_data), varlist, mode=NC_CLASSIC_MODEL)
  local ncfile = NetCDF.open(string(dirout, "/", outname_data), mode=NC_WRITE)

  NetCDF.putvar(ncfile["xc"],  OffsetArrays.no_offset_view(xc), start= [start[1]], count= [count[1]])
  NetCDF.putvar(ncfile["yc"],  OffsetArrays.no_offset_view(yc), start= [start[2]], count= [count[2]])
  NetCDF.putvar(ncfile["zc"],  OffsetArrays.no_offset_view(z), start= start, count= count)
  NetCDF.putvar(ncfile["h"],  OffsetArrays.no_offset_view(h), start= start2d, count= count2d)
  NetCDF.putvar(ncfile["consump"], OffsetArrays.no_offset_view(consump), start=start4, count=count4consump)
  NetCDF.putvar(ncfile["tr"],  OffsetArrays.no_offset_view(Trwrite), start= start4, count= count4)
  NetCDF.putvar(ncfile["s"],  OffsetArrays.no_offset_view(s), start= start, count= count)
  NetCDF.putvar(ncfile["temp"],  OffsetArrays.no_offset_view(T), start= start, count= count)
  NetCDF.putvar(ncfile["rho"],  OffsetArrays.no_offset_view(rho), start= start, count= count)
  NetCDF.putvar(ncfile["p"],  OffsetArrays.no_offset_view(p), start= start, count= count)
  NetCDF.putvar(ncfile["u"],  OffsetArrays.no_offset_view(u), start= start, count= count)
  NetCDF.putvar(ncfile["v"],  OffsetArrays.no_offset_view(v), start= start, count= count)
  NetCDF.putvar(ncfile["w"],  OffsetArrays.no_offset_view(w), start= start, count= count)
  NetCDF.putvar(ncfile["vor"],  OffsetArrays.no_offset_view(vor), start= start, count= count)
  NetCDF.putvar(ncfile["pv"],  OffsetArrays.no_offset_view(pv), start= start, count= count)
  NetCDF.putvar(ncfile["conv"],  OffsetArrays.no_offset_view(conv), start= start, count= count)
  NetCDF.putvar(ncfile["con100"],  OffsetArrays.no_offset_view(con100), start= start, count= count)

  NetCDF.sync(ncfile)

  # End of writing of cell centers values.
  #---------------------------------------

  
  #---------------------------------------
  # Write face velocities, uf,vf,wf  

  # 3D dimensions for face values  
  countuf = [NI+1, NJ, NK]
  countvf = [NI, NJ+1, NK]
  countwf = [NI, NJ, NK+1]
              
  dimuf = [NcDim("xi-ew", NI+1), NcDim("eta-ew", NJ), NcDim("sigma-ew", NK)]
  dimvf = [NcDim("xi-ns", NI), NcDim("eta-ns", NJ+1), NcDim("sigma-ns", NK)]
  dimwf = [NcDim("xi-tb", NI), NcDim("eta-tb", NJ), NcDim("sigma-tb", NK+1)]

  local varlist_face = [
    NcVar("uf", dimuf, t=Float64),
    NcVar("vf", dimvf, t=Float64),
    NcVar("wf", dimwf, t=Float64),
    NcVar("Kz", dimwf, t=Float64)
  ]
  NetCDF.create(string(dirout, "/", outname_face), varlist_face, mode=NC_CLASSIC_MODEL)
  local ncfile_face = NetCDF.open(string(dirout, "/", outname_face), mode=NC_WRITE)


  NetCDF.putvar(ncfile_face["uf"],  OffsetArrays.no_offset_view(uf), start= start, count= countuf)
  NetCDF.putvar(ncfile_face["vf"],  OffsetArrays.no_offset_view(vf), start= start, count= countvf)
  NetCDF.putvar(ncfile_face["wf"],  OffsetArrays.no_offset_view(wf), start= start, count= countwf)
  NetCDF.putvar(ncfile_face["Kz"],  OffsetArrays.no_offset_view(Kz), start= start, count= countwf)

  NetCDF.sync(ncfile_face)


  # End of writing of cell faces values.
  #-------------------------------------

  return
end
      
