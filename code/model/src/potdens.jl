

function potdens(sl, Tl)
     local s15 = sl * sqrt(sl)

     local rw = 999.842594e0 +  6.793952e-2 * Tl - 9.095290e-3 * Tl^2 + 1.001685e-4 * Tl^3 -1.120083e-6 * Tl^4 + 6.56332e-9 * Tl^5

     return rw + (8.24493e-1 - 4.0899e-3 * Tl + 7.6438e-5 * Tl^2 - 8.2467e-7 * Tl^3 + 5.3875e-9 * Tl^4 ) * sl + (-5.72466e-3 + 1.0227e-4 * Tl - 1.6546e-6 * Tl^2) * s15 + 4.8314e-4*sl*sl
end