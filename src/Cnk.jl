"""
    Cnk(n, k, σ)

Calculation of the Cnk using the direct formula with arbitrary precision. This is unstable for large values of k, more so for small values of σ.

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk(6, 5, 0.5)
[0.23437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 6.66e-1506]
```
"""
function Cnk(n, k, σ::Float64)
    # = Cnk(n, k, σ, 0)
    σ_arb = RR(σ)
    return Cnk(n, k, σ_arb)
end
function Cnk(n, k, σ::arb)
    return 1//fac(RR(k)) * sum([ (-1)^i * binom(k, i, RR) * risingfac(-i * σ,n) for i in 0:k])
end

"""
    noncentral_generalised_factorial_coefficient(n, k, s, r)

Calculation of the noncentral generalised factorial coefficient (Charalambides, 2005) using a recursive relation. This function uses memoization to speed up the recursion.
It is called by Cnk_rec but follows the sign convention of (Charalambides, 2005), which differs from that of Cnk.

References:
- Charalambides, Charalambos A. Combinatorial methods in discrete distributions. Vol. 600. John Wiley & Sons, 2005.

# Examples
```julia-repl
julia> GibbsTypePriors.noncentral_generalised_factorial_coefficient(10, 5, 0.2, 0.)
-40.19355648000003
```
"""
@memoize function noncentral_generalised_factorial_coefficient(n, k, s::T, r)::T where T<:Number
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


# @memoize function noncentral_generalised_factorial_coefficient(n, k, s::arb, r::arb)::arb
#     #Getting ERROR: MethodError: Cannot `convert` an object of type Int64 to an object of type arb, which I am not able to fix. Indeed, there is no default convert method, the bits of precision must be chosen. This is a specialised function to fix this error
#     @assert n >= 0
#     @assert k >= 0
#     if k==0
#         if n==0
#             return arb_1
#         else
#             return risingfac(r, n)
#         end
#     else
#         if k>n
#             return arb_0
#         else
#             return (s * k + r - n + arb_1) * noncentral_generalised_factorial_coefficient(n - 1, k, s, r) + s * noncentral_generalised_factorial_coefficient(n - 1, k - 1, s, r)
#         end
#     end
# end

@memoize function noncentral_generalised_factorial_coefficient(n, k, s::arb, r::arb; arb_0::arb=arb_0, arb_1::arb=arb_1)::arb
    #Getting ERROR: MethodError: Cannot `convert` an object of type Int64 to an object of type arb, which I am not able to fix. This is a specialised function to fix this error
    @assert n >= 0
    @assert k >= 0
    if k==0
        if n==0
            return arb_1
        else
            return risingfac(r, n)
        end
    else
        if k>n
            return arb_0
        else
            return (s * k + r - n + 1) * noncentral_generalised_factorial_coefficient(n - 1, k, s, r; arb_0=arb_0, arb_1=arb_1) + s * noncentral_generalised_factorial_coefficient(n - 1, k - 1, s, r; arb_0=arb_0, arb_1=arb_1)
        end
    end
end

"""
    Cnk_rec(n, k, σ)

Calculation of the Cnk using the recursive formula with arbitrary precision. This seems stable for large values of k, whatever the range of σ, but slower than Cnk(n, k, σ).

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk_rec(6, 5, 0.5)
0.234375
```
"""
function Cnk_rec(n, k, σ)
    # factor_k = factorial(k)
#   sum((-1)^(1:k) * gmp::chooseZ(n = k, k = 1:k) * Rmpfr::pochMpfr(-(1:k) * Gama, n) / factor_k)
  (-1)^(n - k) * noncentral_generalised_factorial_coefficient(n, k, σ, 0)
end
function Cnk_rec(n, k, σ::arb)
    #Same comment as  noncentral_generalised_factorial_coefficient(n, k, s::arb, r::arb)::arb

  (-1)^(n - k) * noncentral_generalised_factorial_coefficient(n, k, σ, arb_0)
end

"""
    Cnk_robust(n, k, σ; verbose = false)

Calculation of the Cnk using either the direct formula with arbitrary precision or the recursive formula. Accuracy is checked for the fast direct computation, and if it is too low (roughly lower than double precision), the function switches to recursive computation of the Cnk.

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk_robust(6, 5, 0.5)
[0.23437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 +/- 6.66e-1506]
```
"""
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
