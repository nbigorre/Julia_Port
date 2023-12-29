

function evalrho_sT(rhonew, nl2)

   for k in 0:NK+1
      for j in 0:NJ+1
         for i in 0:NI+1
            local se = s[i, j, k, nl2]
            local Te = T[i, j, k, nl2]
            rhonew[i, j, k] = potdens(se, Te)
         end
      end
   end

end