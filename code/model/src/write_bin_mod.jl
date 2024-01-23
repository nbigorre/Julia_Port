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

function w_pickup(filename)
   #use header, only : h,u,v,w,uf,vf,wf,T,S,ntr,Tr,p,gradhn,hxn,hyn,rc_kind
   #character(len=*) :: filename
   open(filename, "w") do file
      write(file, h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn)
   end
   open(string(filename,".meta"), "w") do file
      write(file, string("h ", size(h,1), " ",size(h,2)," 0 0\n"))
      write(file, string("u ",size(u,1)," ",size(u,2)," ",size(u,3)," ",size(u,4), "\n"))
      write(file, string("v ",size(v,1)," ",size(v,2)," ",size(v,3)," ",size(v,4), "\n"))
      write(file, string("w ",size(w,1)," ",size(w,2)," ",size(w,3)," ",size(w,4), "\n"))
      write(file, string("uf ",size(uf,1)," ",size(uf,2)," ",size(uf,3)," ","0", "\n"))
      write(file, string("vf ",size(vf,1)," ",size(vf,2)," ",size(vf,3)," ","0", "\n"))
      write(file, string("wf ",size(wf,1)," ",size(wf,2)," ",size(wf,3)," ","0", "\n"))
      write(file, string("T ",size(T,1)," ",size(T,2)," ",size(T,3)," ",size(T,4), "\n"))
      write(file, string("S ",size(S,1)," ",size(S,2)," ",size(S,3)," ",size(S,4), "\n"))
      write(file, string("Tr ",ntr," ",size(S,1)," ",size(S,2)," ",size(S,3)," ",size(S,4), "\n"))
      write(file, string("p ",size(p,1)," ",size(p,2)," ",size(p,3), "\n"))
      write(file, string("gradhn ",size(gradhn,1)," ",size(gradhn,2)," ",size(gradhn,3), "\n"))
      write(file, string("hxn ",size(hxn,1)," ",size(hxn,2)," ",size(hxn,3), "\n"))
      write(file, string("hyn ",size(hyn,1)," ",size(hyn,2)," ",size(hyn,3), "\n"))
   end
end

function r_pickup(step)
   #use header, only : h,u,v,w,uf,vf,wf,T,S,Tr,p,gradhn,hxn,hyn,rc_kind,dirout
   #integer :: step
   #character(len=10) :: stepchar

   open(filename, "r") do file
      h = read(file, typeof(h))
      u = read(file, typeof(u))
      v = read(file, typeof(v))
      w = read(file, typeof(w))
      uf = read(file, typeof(uf))
      vf = read(file, typeof(vf))
      wf = read(file, typeof(wf))
      T = read(file, typeof(T))
      S = read(file, typeof(S))
      Tr = read(file, typeof(Tr))
      p = read(file, typeof(p))
      gradhn = read(file, typeof(gradhn))
      hxn = read(file, typeof(hxn))
      hyn = read(file, typeof(hyn))
   end

   println("# pickup at step ", step)
end
