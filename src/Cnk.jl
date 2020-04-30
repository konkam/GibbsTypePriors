function Cnk(n, k, σ; prec = 5000)
    RF = RealField(prec)
    σ_arb = RF(σ)
    return 1//fac(RF(k)) * sum([ (-1)^i * binom(k, i, RF) * risingfac(-i * σ_arb,n) for i in 0:k])
end


@memoize function noncentral_generalised_factorial_coefficient(n, k, s::T, r; prec = 5000)::T where T
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

function noncentral_generalised_factorial_coefficient(n, k, s::arb, r; prec = 5000)::arb
    RF = RealField(prec)
    return noncentral_generalised_factorial_coefficient_in(n, k, s, r, RF)
end

@memoize function noncentral_generalised_factorial_coefficient_in(n, k, s::arb, r, RF)::arb
    #Getting ERROR: MethodError: Cannot `convert` an object of type Int64 to an object of type arb, which I am not able to fix. This is a specialised function to fix this error
    @assert n >= 0
    @assert k >= 0
    if k==0
        if n==0
            return RF(1)
        else
            return risingfac(RF(r), n)
        end
    else
        if k>n
            return RF(0)
        else
            return (s * k + r - n + 1) * noncentral_generalised_factorial_coefficient_in(n - 1, k, s, r, RF) + s * noncentral_generalised_factorial_coefficient_in(n - 1, k - 1, s, r, RF)
        end
    end
end

function Cnk_rec(n, k, σ; prec = 5000)
    # factor_k = factorial(k)
#   sum((-1)^(1:k) * gmp::chooseZ(n = k, k = 1:k) * Rmpfr::pochMpfr(-(1:k) * Gama, n) / factor_k)
  (-1)^(n - k) * noncentral_generalised_factorial_coefficient(n, k, σ, 0; prec = 5000)
end
function Cnk_rec(n, k, σ::arb)
    #Same comment as  noncentral_generalised_factorial_coefficient(n, k, s::arb, r::arb)::arb

  (-1)^(n - k) * noncentral_generalised_factorial_coefficient(n, k, σ, arb_0)
end

function Cnk_robust(n, k, σ; verbose = false)
    res = Cnk(n, k, σ)
    if has_reasonable_precision(res)
        return res
    else
        rec_res = Cnk_rec(n, k, RR(σ))
        if verbose
            println("switching to recursive formula for Cnk")
        end
        if has_reasonable_precision(rec_res)
            return rec_res
        else
            error("Unable to compute Cnk with reasonable precision (see function has_reasonable_precision() )")
        end
    end
end
