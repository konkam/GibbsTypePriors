"""
    Cnk(n, k, σ, prec)

Calculation of the Cnk using the direct formula with arbitrary precision, using a specified number of bits (prec). This is unstable for large values of k, more so for small values of σ.

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk(6, 5, 0.5, 200)
[0.234375000000000000000000000000000000000000000000000000000000 +/- 5.85e-61]
```
"""
function Cnk(n, k, σ, prec)
    RF = RealField(prec)
    σ_arb = RF(σ)
    return 1//fac(RF(k)) * sum([ (-1)^i * binom(k, i, RF) * risingfac(-i * σ_arb,n) for i in 0:k])
end

function noncentral_generalised_factorial_coefficient_prec(n, k, s, r, prec)::arb
    RF = RealField(prec)
    arb_1 = RF(1)
    arb_0 = RF(0)
    arb_r = RF(r)
    arb_s = RF(s)
    return noncentral_generalised_factorial_coefficient(n, k, arb_s, arb_r; arb_0=arb_0, arb_1=arb_1)
end

function Cnk_rec(n, k, σ, prec)
    # factor_k = factorial(k)
#   sum((-1)^(1:k) * gmp::chooseZ(n = k, k = 1:k) * Rmpfr::pochMpfr(-(1:k) * Gama, n) / factor_k)
  (-1)^(n - k) * noncentral_generalised_factorial_coefficient_prec(n, k, σ, 0, prec)
end