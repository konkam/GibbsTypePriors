using SpecialFunctions, Nemo, StatsFuns
import Nemo.binom, Nemo.gamma

prec = 5000 ## Increase for better precision
RR = RealField(prec)
CC = ComplexField(prec)
binom(n::Int64, k::Int64, r::Nemo.ArbField) = binom(convert(UInt64, n),convert(UInt64, k), r)
gamma(x::Int64) = Nemo.gamma(RR(x))
