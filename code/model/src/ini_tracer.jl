function tracerinit(stepl) 
#     initializes tracer fields                                         
#     TRACER 1                                                          
#     =========                                                         
      local it=1 
      Tr .= 0e0
      Tr[it,:,:,1,0] .=1e0
end
