using SpecialFunctions, Nemo, StatsFuns
import Nemo.binom, Nemo.gamma

prec = 1000 ## Increase for better precision
RR = RealField(prec)
CC = ComplexField(prec)
binom(n::Int64, k::Int64, r::Nemo.ArbField) = binom(convert(UInt64, n),convert(UInt64, k), r)
gamma(x::Int64) = Nemo.gamma(RR(x))

function Cnk(n, k, σ)
    # = Cnk(n, k, σ, 0)
    σ_arb = RR(σ)
    return 1//fac(RR(k)) * sum([ (-1)^i * binom(k, i, RR) * risingfac(-i * σ_arb,n) for i in 0:k])
end
