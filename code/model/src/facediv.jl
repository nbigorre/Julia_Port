function facediv(dtimel)
  local maxdiv= 0e0 
  for k in 1:NK 
    for j in 1:NJ 
      for i in 1:NI 
        ufdx= (uf[i,j,k]-uf[i-1,j,k]) 
        vfdy= (vf[i,j,k]-vf[i,j-1,k]) 
        wfdz= (wf[i,j,k]-wf[i,j,k-1]) 
        div= abs(ufdx+ vfdy + wfdz) 
        if (div > maxdiv) 
              maxdiv=div;
        end
     end
   end
 end
 return maxdiv
end