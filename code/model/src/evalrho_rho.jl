

function evalrho_rho(rhonew, n)

   for k in 0:NK+1
      for j in 0:NJ+1
         for i in 0:NI+1
            rhonew[i, j, k] = s[i, j, k, n]
         end
      end
   end

end