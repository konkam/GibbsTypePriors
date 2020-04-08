function Cnk(n, k, σ)
    # = Cnk(n, k, σ, 0)
    σ_arb = RR(σ)
    return 1//fac(RR(k)) * sum([ (-1)^i * binom(k, i, RR_in) * risingfac(-i * σ_arb,n) for i in 0:k])
end


@memoize function noncentral_generalised_factorial_coefficient(n, k, s::T, r)::T where T
    @assert n >= 0
    @assert k >= 0
    if k==0
        if n==0
            return 1
        else
            return risingfac(r, n)
        end
    else
        if k>n
            return 0
        else
            return (s * k + r - n + 1) * noncentral_generalised_factorial_coefficient(n - 1, k, s, r) + s * noncentral_generalised_factorial_coefficient(n - 1, k - 1, s, r)
        end
    end
end


@memoize function noncentral_generalised_factorial_coefficient(n, k, s::arb, r)
    #Getting ERROR: MethodError: Cannot `convert` an object of type Int64 to an object of type arb, which I am not able to fix. This is a specialised function to fix this error
    @assert n >= 0
    @assert k >= 0
    if k==0
        if n==0
            return RR(1)
        else
            return risingfac(r, n)
        end
    else
        if k>n
            return RR(0)
        else
            return (s * k + r - n + 1) * noncentral_generalised_factorial_coefficient(n - 1, k, s, r) + s * noncentral_generalised_factorial_coefficient(n - 1, k - 1, s, r)
        end
    end
end

function Cnk_rec(n, k, σ)
    # factor_k = factorial(k)
#   sum((-1)^(1:k) * gmp::chooseZ(n = k, k = 1:k) * Rmpfr::pochMpfr(-(1:k) * Gama, n) / factor_k)
  (-1)^(n - k) * noncentral_generalised_factorial_coefficient(n, k, σ, 0)
end
