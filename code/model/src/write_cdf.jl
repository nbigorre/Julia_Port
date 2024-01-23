function write_cdf(step, n)

  #  USE header,ONLY : NI,NJ,NK,ntr,nconsume,dirout,out1d_int,out2d_int,out3d_int,rc_kind,pvt,pv1,pv2,pv3
  #  IMPLICIT NONE
  #
  #  INTEGER :: step,counter_2d,counter_3d,counter_1d,ksurf,islice,jslice,imooring,jmooring,n
  #  REAL :: sigma,z

  # out1d_int  - frequency of 1d output
  # out2d_int  - frequency of 2d output
  # out3d_int  - frequency of 3d output

  # 1D output
  if (mod(step, out1d_int) == 0)
    local counter_1d = div(step, out1d_int) + 1
    local imooring = div(NI, 2)
    local jmooring = div(NJ, 2)
    write_cdf_1D_mooring(imooring, jmooring, counter_1d)
    #@ccall "./PSOM_LIB.so".write_cdf_1d_mooring_(Ref(imooring)::Ref{Int}, Ref(jmooring)::Ref{Int}, Ref(counter_1d)::Ref{Int})::Cvoid
  end

  # 2D output
  if (mod(step, out2d_int) == 0)
    local counter_2d = div(step, out2d_int) + 1

    diag_pv(n)
    #@ccall "./PSOM_LIB.so".diag_pv_(Ref(n)::Ref{Int})::Cvoid
    pvt .= pv1 + pv2 + pv3

    local ksurf = NK
    write_cdf_2D_sigma(ksurf, counter_2d, n)
    #@ccall "./PSOM_LIB.so".write_cdf_2d_sigma_(Ref(ksurf)::Ref{Int}, Ref(counter_2d)::Ref{Int}, Ref(n)::Ref{Int})::Cvoid

    local islice = 5
    write_cdf_2D_x(islice, counter_2d, n)
    #@ccall "./PSOM_LIB.so".write_cdf_2d_x_(Ref(islice)::Ref{Int}, Ref(counter_2d)::Ref{Int}, Ref(n)::Ref{Int})::Cvoid
  end

  # 3D output
  if (mod(step, out3d_int) == 0)
    local counter_3d = div(step, out3d_int) + 1
    write_cdf_3D(step, n)
    #@ccall "./PSOM_LIB.so".write_cdf_3d_(Ref(step)::Ref{Int}, Ref(n)::Ref{Int})::Cvoid
    write_cdf_3D_strain(step, n)
    #@ccall "./PSOM_LIB.so".write_cdf_3d_strain_(Ref(step)::Ref{Int}, Ref(n)::Ref{Int})::Cvoid

  end

end