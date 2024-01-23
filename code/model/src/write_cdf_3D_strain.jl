function write_cdf_3D_strain(step, n)
   #!---------------------
   #! writes out the terms for the strainrate tenstor, namely ux,vx,wx, uy,vy,wy,
   #! uz,vz,wz
   #  USE header, ONLY : NI,NJ,NK,xc,yc,zc,u,v,w,uf,vf,wf,Jac,rho,UL,LEN,WL,DL,DLinv,gpr,R0,dirout,rc_kind,&
   #   & ux,uy,vx,vy,wx,wy,wz,time_seconds
   #
   ##include "netcdf.inc"
   #  integer dims(3),start(3),count(3),nstp,step
   #
   #  integer i,idv,idvdy,j,k,n,kp1,km1
   #  REAL(kind=4) ::  dudx(NI,NJ,NK),dudy(NI,NJ,NK),dudz(NI,NJ,NK),dvdx(NI,NJ,NK),     &
   #       dvdy(NI,NJ,NK),dvdz(NI,NJ,NK),dwdx(NI,NJ,NK),dwdy(NI,NJ,NK),       &
   #       dwdz(NI,NJ,NK),rdx(NI,NJ,NK),rdy(NI,NJ,NK),rdz(NI,NJ,NK)
   #  REAL(kind=rc_kind) ::  LENinv,cst,fac
   #
   #  REAL(kind=rc_kind) ::  rcode
   #  integer :: idvx
   #  REAL ::time_days
   #
   #
   #  character(LEN=150) outname
   #
   #  DATA start /1, 1, 1/
   #
   #
   # integer :: iddatfile, idigit, idudx, idudy, idudz, idvbysection, idvby, idvbz, idvcon, idvcy, idvbx
   #  integer :: idvcz, idvc, idvdivfeddy, idvdivfreyn, idvdx, idvdz, idvd, idvfb, idvh, idvn2bar
   #  integer :: idvn2, idvnsq100m, idvnsq30m, idvpe,idvpsiv,idvpsiw,idvpv, idvp,idvrhbar,idvrho,idvrnk
   #  integer :: idvstrain,idvstress,idvstr,idvs,idvtbar,idvtemp,idvtim,idvtr,idvt,idvu,idvvb,idvvc,idvvor,idvv
   #  integer :: idvwb,idvwc,idvwpv,idvw,idvy,idvzsave,idvz,idwdx,idwdy,idwdz,ipos
   #
   #  REAL(kind=rc_kind) ::  vol
   #
   local start = [1, 1, 1]

   local dudx = zeros(rc_kind, NI, NJ, NK)
   local dudy = zeros(rc_kind, NI, NJ, NK)
   local dudz = zeros(rc_kind, NI, NJ, NK)
   local dvdx = zeros(rc_kind, NI, NJ, NK)
   local dvdy = zeros(rc_kind, NI, NJ, NK)
   local dvdz = zeros(rc_kind, NI, NJ, NK)
   local dwdx = zeros(rc_kind, NI, NJ, NK)
   local dwdy = zeros(rc_kind, NI, NJ, NK)
   local dwdz = zeros(rc_kind, NI, NJ, NK)
   local rdx = zeros(rc_kind, NI, NJ, NK)
   local rdy = zeros(rc_kind, NI, NJ, NK)
   local rdz = zeros(rc_kind, NI, NJ, NK)

   local time_days = time_seconds / 86400e0

   #!     name the output file
   #! First caculate ux,vx,wx

   for k in 1:NK
      for j in 1:NJ
         for i in 1:NI
            dudx[i, j, k] = (uf[i, j, k] - uf[i-1, j, k]) / Jac[i, j, k]
            dvdy[i, j, k] = (vf[i, j, k] - vf[i, j-1, k]) / Jac[i, j, k]
            dwdz[i, j, k] = (wf[i, j, k] - wf[i, j, k-1]) / Jac[i, j, k]
         end
      end
   end
   # Now convert to dimensional form
   local dudx = dudx * UL / LEN
   local dvdy = dvdy * UL / LEN
   local dwdz = dwdz * WL / DL


   #     periodic_e-w boundaries
   for k in 0:NK+1
      for j in 0:NJ+1
         u[0, j, k, n] = u[NI, j, k, n]
         v[0, j, k, n] = v[NI, j, k, n]
         w[0, j, k, n] = w[NI, j, k, n]
         rho[0, j, k] = rho[NI, j, k]

         u[NI+1, j, k, n] = u[1, j, k, n]
         v[NI+1, j, k, n] = v[1, j, k, n]
         w[NI+1, j, k, n] = w[1, j, k, n]
         rho[NI+1, j, k] = rho[1, j, k]
      end
   end

   local DLinv = 1.0 / DL
   local LENinv = 1.0 / LEN
   local cst = 2e0 * gpr * 10e0 / R0

   #     for diff in j (solid bndry)
   for k in 1:NK
      for i in 1:NI
         u[i, 0, k, n] = 2e0 * u[i, 1, k, n] - u[i, 2, k, n]
         v[i, 0, k, n] = 2e0 * v[i, 1, k, n] - v[i, 2, k, n]
         w[i, 0, k, n] = 2e0 * w[i, 1, k, n] - w[i, 2, k, n]
         rho[i, 0, k] = 2e0 * rho[i, 1, k] - rho[i, 2, k]

         u[i, NJ+1, k, n] = 2e0 * u[i, NJ, k, n] - u[i, NJ-1, k, n]
         v[i, NJ+1, k, n] = 2e0 * v[i, NJ, k, n] - v[i, NJ-1, k, n]
         w[i, NJ+1, k, n] = 2e0 * w[i, NJ, k, n] - w[i, NJ-1, k, n]
         rho[i, NJ+1, k] = 2e0 * rho[i, NJ, k] - rho[i, NJ-1, k]
      end
   end


   for k in 1:NK
      local kp1 = k + 1
      local km1 = k - 1
      if (k == NK)
         kp1 = NK
      end
      if (k == 1)
         km1 = 1
      end
      local fac = 1e0 / rc_kind(kp1 - km1)
      for j in 1:NJ
         for i in 1:NI
            dvdx[i, j, k] = LENinv * UL * (0.5 * (v[i+1, j, k, n] - v[i-1, j, k, n]) * ux[i, j] + 0.5 * (v[i, j+1, k, n] - v[i, j-1, k, n]) * vx[i, j] + fac * (v[i, j, kp1, n] - v[i, j, km1, n]) * wx[i, j, k])
            dwdx[i, j, k] = LENinv * WL * (0.5 * (w[i+1, j, k, n] - w[i-1, j, k, n]) * ux[i, j] + 0.5 * (w[i, j+1, k, n] - w[i, j-1, k, n]) * vx[i, j] + fac * (w[i, j, kp1, n] - w[i, j, km1, n]) * wx[i, j, k])
            dudy[i, j, k] = LENinv * UL * (0.5 * (u[i+1, j, k, n] - u[i-1, j, k, n]) * uy[i, j] + 0.5 * (u[i, j+1, k, n] - u[i, j-1, k, n]) * vy[i, j] + fac * (u[i, j, kp1, n] - u[i, j, km1, n]) * wy[i, j, k])
            dwdy[i, j, k] = LENinv * WL * (0.5 * (w[i+1, j, k, n] - w[i-1, j, k, n]) * uy[i, j] + 0.5 * (w[i, j+1, k, n] - w[i, j-1, k, n]) * vy[i, j] + fac * (w[i, j, kp1, n] - w[i, j, km1, n]) * wy[i, j, k])

            rdx[i, j, k] = LENinv * (0.5 * (rho[i+1, j, k] - rho[i-1, j, k]) * ux[i, j] + 0.5 * (rho[i, j+1, k] - rho[i, j-1, k]) * vx[i, j] + fac * (rho[i, j, kp1] - rho[i, j, km1]) * wx[i, j, k])

            rdy[i, j, k] = LENinv * (0.5 * (rho[i+1, j, k] - rho[i-1, j, k]) * uy[i, j] + 0.5 * (rho[i, j+1, k] - rho[i, j-1, k]) * vy[i, j] + fac * (rho[i, j, kp1] - rho[i, j, km1]) * wy[i, j, k])

            rdz[i, j, k] = fac * (rho[i, j, kp1] - rho[i, j, km1]) * wz[i, j, k] * DLinv
            dudz[i, j, k] = UL * fac * (u[i, j, kp1, n] - u[i, j, km1, n]) * wz[i, j, k] * DLinv
            dvdz[i, j, k] = UL * fac * (v[i, j, kp1, n] - v[i, j, km1, n]) * wz[i, j, k] * DLinv

         end
      end
   end

   #     name the output file

   local outname = string("strain_", lpad(step, 5, "0"), ".cdf")
   local count = [NI, NJ, NK]


   local dims = [NcDim("x", NI), NcDim("y", NJ), NcDim("sigma", NK)]

   local varlist = [
      NcVar("xc", [dims[1]], t=Float32),
      NcVar("yc", [dims[2]], t=Float32),
      NcVar("zc", dims, t=Float32),
      NcVar("udx", dims, t=Float32),
      NcVar("udy", dims, t=Float32),
      NcVar("udz", dims, t=Float32),
      NcVar("vdx", dims, t=Float32),
      NcVar("vdy", dims, t=Float32),
      NcVar("vdz", dims, t=Float32),
      NcVar("wdx", dims, t=Float32),
      NcVar("wdy", dims, t=Float32),
      NcVar("wdz", dims, t=Float32),
      NcVar("Nsq", dims, t=Float32)#,
      #NcVar("vor1", dims, t=Float32),
      #NcVar("vor2", dims, t=Float32),
      #NcVar("vor3", dims, t=Float32),
      #NcVar("pv1", dims, t=Float32),
      #NcVar("pv2", dims, t=Float32),
      #NcVar("pv3", dims, t=Float32)
   ]

   NetCDF.create(string(dirout, "/", outname), varlist, mode=NC_CLASSIC_MODEL)

   local ncfile = NetCDF.open(string(dirout, "/", outname), mode=NC_WRITE)

   NetCDF.putvar(ncfile["xc"], OffsetArrays.no_offset_view(xc), start=[start[1]], count=[count[1]])
   NetCDF.putvar(ncfile["yc"], OffsetArrays.no_offset_view(yc), start=[start[2]], count=[count[2]])
   NetCDF.putvar(ncfile["zc"], OffsetArrays.no_offset_view(zc), start=start, count=count)
   NetCDF.putvar(ncfile["udx"], OffsetArrays.no_offset_view(dudx), start=start, count=count)
   NetCDF.putvar(ncfile["udy"], OffsetArrays.no_offset_view(dudy), start=start, count=count)
   NetCDF.putvar(ncfile["udz"], OffsetArrays.no_offset_view(dudz), start=start, count=count)
   NetCDF.putvar(ncfile["vdx"], OffsetArrays.no_offset_view(dvdx), start=start, count=count)
   NetCDF.putvar(ncfile["vdy"], OffsetArrays.no_offset_view(dvdy), start=start, count=count)
   NetCDF.putvar(ncfile["vdz"], OffsetArrays.no_offset_view(dvdz), start=start, count=count)
   NetCDF.putvar(ncfile["wdx"], OffsetArrays.no_offset_view(dwdx), start=start, count=count)
   NetCDF.putvar(ncfile["wdy"], OffsetArrays.no_offset_view(dwdy), start=start, count=count)
   NetCDF.putvar(ncfile["wdz"], OffsetArrays.no_offset_view(dwdz), start=start, count=count)
   #NetCDF.putvar(ncfile["vor1"],  OffsetArrays.no_offset_view(vor1), start= start, count= count)
   #NetCDF.putvar(ncfile["vor2"],  OffsetArrays.no_offset_view(vor2), start= start, count= count)
   #NetCDF.putvar(ncfile["vor3"],  OffsetArrays.no_offset_view(vor3), start= start, count= count)
   #NetCDF.putvar(ncfile["pv1"],  OffsetArrays.no_offset_view(pv1), start= start, count= count)
   #NetCDF.putvar(ncfile["pv2"],  OffsetArrays.no_offset_view(pv2), start= start, count= count)
   #NetCDF.putvar(ncfile["pv3"],  OffsetArrays.no_offset_view(pv3), start= start, count= count)


   NetCDF.sync(ncfile)

end