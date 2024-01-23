function write_cdf_2D_sigma(ksurf, counter_2d, n)

   #  !     ------------------------------------------------------------------
   #USE header,ONLY : NI,NJ,NK,ntr,h,consump,Tr,s,T,rho,u,v,w,p,pvt,vor,strain,freqN2,   &
   # & xc,yc,zc,DL,LEN,Jac,nconsume,time_seconds,dirout,rc_kind
   #  !     periodically writes out the data on the k-th surface              
   ##include "netcdf.inc"                                                                        
   #!      INCLUDE '/sw/lib/netcdf-gfortran/include/netcdf.inc'             
   #                                                                        
   #      integer i,j,k,n,it,keuphotic,kdum,ksurf,k30,k100  
   #      REAL(kind=rc_kind) :: zmax,zeuphotic,dep
   #                                                                        
   #      REAL(kind=4) ::  xsurf(NI),ysurf(NJ),zsurf(NI,NJ),ssurf(NI,NJ),              &
   #     &     tempsurf(NI,NJ),Trsurf(NI,NJ,ntr),rsurf(NI,NJ),rNK(NI,NJ),     &
   #     &     usurf(NI,NJ),vsurf(NI,NJ),wsurf(NI,NJ),hsurf(NI,NJ),          &
   #     &psurf(NI,NJ), pvsurf(NI,NJ),        &
   #     &     vorsurf(NI,NJ),strainsurf(NI,NJ),n2surf(NI,NJ),              &
   #     &     nsq30m(NI,NJ),nsq100m(NI,NJ),                                &
   #     &     integcons(NI,NJ,nconsume)                             
   #!     ,integconsBBL(NI,NJ,ntr)                                          
   #!                                          
   #      REAL :: time_days                             
   #       INTEGER :: counter_2d
   #
   # integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
   #  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
   #  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
   #  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
   #  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,ipos
   # 
   #  REAL(kind=rc_kind) ::  rcode,vol
   #  integer :: idvx
   #
   #                                                                        
   #      character(LEN=150) :: outname 
   #                                                                        
   #      integer start(3),count(3),dims(3),start2d(2),count2d(2),dims2d(2),&
   #     &     dimstr(4),dimsconstr(4),starttr(4),counttr(4),countconstr(4) &
   #     &     ,dimstim,counttim   

   local xsurf = zeros(Float32, NI)
   local ysurf = zeros(Float32, NJ)
   local zsurf = zeros(Float32, NI, NJ)
   local ssurf = zeros(Float32, NI, NJ)
   local tempsurf = zeros(Float32, NI, NJ)
   local Trsurf = zeros(Float32, NI, NJ, ntr)
   local rsurf = zeros(Float32, NI, NJ)
   local rNK = zeros(Float32, NI, NJ)
   local usurf = zeros(Float32, NI, NJ)
   local vsurf = zeros(Float32, NI, NJ)
   local wsurf = zeros(Float32, NI, NJ)
   local hsurf = zeros(Float32, NI, NJ)
   local psurf = zeros(Float32, NI, NJ)
   local pvsurf = zeros(Float32, NI, NJ)
   local vorsurf = zeros(Float32, NI, NJ)
   local strainsurf = zeros(Float32, NI, NJ)
   local n2surf = zeros(Float32, NI, NJ)
   local nsq30m = zeros(Float32, NI, NJ)
   local nsq100m = zeros(Float32, NI, NJ)
   local integcons = zeros(Float32, NI, NJ, nconsume)
   #local integconsBBL = zeros(Float32, NI, NJ, ntr)
   local time_days = Float32(@fortGet("time_seconds", rc_kind) / 86400e0)
   local zeuphotic = -100e0 / DL
   local keuphotic = 1
   for k in 1:NK
      if (zc[div(NI, 2), div(NJ, 2), k] > zeuphotic)
         keuphotic = k
         break
      end
   end

   for j in 1:NJ
      for i in 1:NI
         for it in 1:nconsume
            integcons[i, j, it] = 0.0
            #               integconsBBL[i,j,it]= 0.0                               
         end
         #    sum over euphotic zone (trinit is zero here)                      
         #            for kdum in keuphotic:NK                                      
         for kdum in keuphotic:NK
            for it in 1:nconsume
               integcons[i, j, it] = integcons[i, j, it] + 0.5e0 * (consump[i, j, kdum, it] + abs(consump[i, j, kdum, it])) * Jac[i, j, kdum] * DL * LEN * LEN
               #     integcons is the consump in milimols/timestep if                  
               #     consump is in mili-mols/m^3/timestep                              
            end
         end
         #            for kdum in 2:3                                                
         #               for it in 1:ntr                                             
         #                  integconsBBL[i,j,it]= integconsBBL[i,j,it] +0.5e0*(consump[i,j,kdum,it]+dabs(consump[i,j,kdum,it]))*Jac[i,j,kdum]*DL*LEN*LEN                           
         #     integcons is the consump in milimols/s if consump is in mili-mols/
         #               end
         #     end
      end
   end

   for j in 1:NJ
      for i in 1:NI
         zsurf[i, j] = zc[i, j, ksurf] * DL
         #xsurf[i, j] = xc[i]
         #ysurf[i, j] = yc[j]
         hsurf[i, j] = h[i, j]
         ssurf[i, j] = s[i, j, ksurf, n]
         tempsurf[i, j] = T[i, j, ksurf, n]
         rsurf[i, j] = rho[i, j, ksurf]
         rNK[i, j] = rho[i, j, NK]
         usurf[i, j] = u[i, j, ksurf, n]
         vsurf[i, j] = v[i, j, ksurf, n]
         wsurf[i, j] = w[i, j, ksurf, n]
         psurf[i, j] = p[i, j, ksurf]
         pvsurf[i, j] = pvt[i, j, ksurf]
         vorsurf[i, j] = vor[i, j, ksurf]
         strainsurf[i, j] = strain[i, j, ksurf]
         n2surf[i, j] = freqN2[i, j, ksurf]
         for it in 1:ntr
            Trsurf[i, j, it] = Tr[it, i, j, ksurf, n]
         end
      end
   end
   for i in 1:NI
      xsurf[i] = xc[i]
   end
   for j in 1:NJ
      ysurf[j] = yc[j]
   end

   #     write average N2 in upper 0-30m and in 30-130 m                   
   #     write to a netcdf file                                            

   local k30 = 0
   local k100 = 0
   for k in NK:-1:1
      dep = -DL * zc[div(NI, 2), div(NJ, 2), k]
      if (dep > 30.0)
         k30 = k + 1
         break
      end
   end
   for k in NK:-1:1
      dep = -DL * zc[div(NI, 2), div(NJ, 2), k]
      if (dep > 100.0)
         k100 = k + 1
         break
      end
   end
   for j in 1:NJ
      for i in 1:NI
         vol = 0e0
         nsq30m[i, j] = 0e0
         for k in NK:-1:k30
            nsq30m[i, j] = nsq30m[i, j] + freqN2[i, j, k] * Jac[i, j, k]
            vol = vol + Jac[i, j, k]
         end
         nsq100m[i, j] = nsq30m[i, j]
         nsq30m[i, j] = nsq30m[i, j] / vol
         for k in k30+1:k100
            nsq100m[i, j] = nsq100m[i, j] + freqN2[i, j, k] * Jac[i, j, k]
            vol = vol + Jac[i, j, k]
         end
         nsq100m[i, j] = nsq100m[i, j] / vol
      end
   end

   local outname = string("zslice_", lpad(ksurf, 3, "0"), ".cdf")

   if (counter_2d == 1)
      #local dims2d = [
      #   NcDim("xi", NI),
      #   NcDim("eta", NJ)
      #]

      local dims = [
         NcDim("x", NI),
         NcDim("y", NJ),
         NcDim("time", 0, unlimited=true)
      ]

      #local dimstr = [
      #   NcDim("x", NI),
      #   NcDim("y", NJ),
      #   NcDim("ntr", ntr),
      #   NcDim("time", 0, unlimited=true)
      #]

      local dimstr = [
         dims[1],
         dims[2],
         NcDim("ntr", ntr),
         dims[3]
      ]

      local dimsconstr = [
         dims[1],
         dims[2],
         NcDim("nconsume", nconsume),
         dims[3]
      ]

      local dimstim = dims[3]
      local varlist = [
         NcVar("xc", dims[1], t=Float32),
         NcVar("yc", dims[2], t=Float32),
         NcVar("zc", dims[1:2], t=Float32),
         NcVar("day", dimstim, t=Float32),
         NcVar("h", dims, t=Float32),
         NcVar("Tr", dimstr, t=Float32),
         NcVar("s", dims, t=Float32),
         NcVar("temp", dims, t=Float32),
         NcVar("rho", dims, t=Float32),
         NcVar("rhosurf", dims, t=Float32),
         NcVar("u", dims, t=Float32),
         NcVar("v", dims, t=Float32),
         NcVar("w", dims, t=Float32),
         NcVar("p", dims, t=Float32),
         NcVar("pv", dims, t=Float32),
         NcVar("vor", dims, t=Float32),
         NcVar("strain", dims, t=Float32),
         NcVar("n2", dims, t=Float32),
         NcVar("nsq30m", dims, t=Float32),
         NcVar("nsq100m", dims, t=Float32),
         NcVar("integcons", dimsconstr, t=Float32)
         #NcVar("integconsBBL", dimstr, t=Float32)
      ]
      NetCDF.create(string(dirout, "/", outname), varlist, mode=NC_CLASSIC_MODEL)


   end

   local counttim = 1

   local count2d = [NI, NJ]

   local count = [NI, NJ, 1]

   local counttr = [NI, NJ, ntr, 1]

   local countconstr = [NI, NJ, nconsume, 1]

   local start2d = [1, 1]

   local start = [1, 1, counter_2d]

   local starttr = [1, 1, 1, counter_2d]


   local ncfile = NetCDF.open(string(dirout, "/", outname), mode=NC_WRITE)
   if (counter_2d == 1)

      NetCDF.putvar(ncfile["xc"], xsurf, start=[start[1]], count=[count[1]])
      NetCDF.putvar(ncfile["yc"], ysurf, start=[start[2]], count=[count[2]])
      NetCDF.putvar(ncfile["zc"], zsurf, start=start2d, count=count2d)

   end
   NetCDF.putvar(ncfile["day"], [time_days], start=[start[3]], count=[counttim])
   NetCDF.putvar(ncfile["h"], [time_days], start=start, count=count)
   NetCDF.putvar(ncfile["h"], [time_days], start=start, count=count)
   NetCDF.putvar(ncfile["Tr"], Trsurf, start=starttr, count=counttr)

   NetCDF.putvar(ncfile["s"], ssurf, start=start, count=count)
   NetCDF.putvar(ncfile["temp"], tempsurf, start=start, count=count)
   NetCDF.putvar(ncfile["rho"], rsurf, start=start, count=count)
   NetCDF.putvar(ncfile["rhosurf"], rNK, start=start, count=count)
   NetCDF.putvar(ncfile["u"], usurf, start=start, count=count)
   NetCDF.putvar(ncfile["v"], vsurf, start=start, count=count)
   NetCDF.putvar(ncfile["w"], wsurf, start=start, count=count)
   #      NetCdf.putVar(ncfile["p"], psurf, start=start, count=count)          
   NetCDF.putvar(ncfile["pv"], pvsurf, start=start, count=count)
   NetCDF.putvar(ncfile["vor"], vorsurf, start=start, count=count)
   NetCDF.putvar(ncfile["strain"], strainsurf, start=start, count=count)
   NetCDF.putvar(ncfile["n2"], n2surf, start=start, count=count)
   NetCDF.putvar(ncfile["nsq30m"], nsq30m, start=start, count=count)
   NetCDF.putvar(ncfile["nsq100m"], nsq100m, start=start, count=count)
   NetCDF.putvar(ncfile["integcons"], nsq100m, start=starttr, count=countconstr)
   #NetCdf.putVar(ncfile["integconsBBL"], nsq100m, start=starttr, count=counttr) 

   NetCDF.sync(ncfile)
end