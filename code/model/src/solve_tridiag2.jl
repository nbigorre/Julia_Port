function solve_tridiag(a, b, c, d, x, n)
        local cp = zeros(rc_kind, n)
        local dp = zeros(rc_kind, n)
        cp[1] = c[1] / b[1]
        dp[1] = d[1] / b[1]

        for i in 2:n
                local m = b[i] - cp[i-1] * a[i]
                cp[i] = c[i]/m
                dp[i] = (d[i]-dp[i-1]*a[i])/m
        end

        x[n] = dp[n]

        for i in n-1:-1:1
                x[i] = dp[i]-cp[i]*x[i+1]
        end
end
