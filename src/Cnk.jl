function Cnk(n, k, σ)
    # = Cnk(n, k, σ, 0)
    σ_arb = RR(σ)
    return 1//fac(RR(k)) * sum([ (-1)^i * binom(k, i, RR) * risingfac(-i * σ_arb,n) for i in 0:k])
end
