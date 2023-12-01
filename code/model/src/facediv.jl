function facediv(dtimel)
  local maxdiv= 0e0 
  for k in 1:NK 
    for j in 1:NJ 
      for i in 1:NI 
        ufdx= (uf[i+1,j,k]-uf[i,j,k]) 
        vfdy= (vf[i,j+1,k]-vf[i,j,k]) 
        wfdz= (wf[i,j,k+1]-wf[i,j,k]) 
        div= abs(ufdx+ vfdy + wfdz) 
        if (div > maxdiv) 
              maxdiv=div;
        end
     end
   end
 end
 return maxdiv
end