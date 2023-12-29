function pcorrect(pcorr) 
   @views @inbounds @. p[0:NI+1, 0:NJ+1, 0:NK+1] .+= pcorr[0:NI+1, 0:NJ+1, 0:NK+1]
end                                     
