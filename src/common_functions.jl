using SpecialFunctions, Nemo, StatsFuns, Memoization
import Nemo.binom, Nemo.gamma

"Session-wide precision"
const prec = 5000 ## Increase for better precision

"Converter from real to arbitrary precision real (Arb)"
const RR = RealField(prec)

"Converter from real to arbitrary precision complex (Acb)"
const CC = ComplexField(prec)

binom(n::Int64, k::Int64, r::Nemo.ArbField) = binom(convert(UInt64, n),convert(UInt64, k), r)

# Making sure the factorial function uses the Nemo version of the gamma and returns an Arb
gamma(x::Int64) = Nemo.gamma(RR(x))

# Useful typed constants
"1 with required arbitrary precision"
const arb_1 = RR(1)
"0 with required arbitrary precision"
const arb_0 = RR(0)

"""
  unsigned_Stirling1(n::Integer, k::Integer)

Computation of unsigned Stirling numbers of the first kind, using Memoization and a recursive formula.

# Examples
```julia-repl
julia> GibbsTypePriors.unsigned_Stirling1(10, 5)
269325
```
"""
@memoize function unsigned_Stirling1(n::Integer, k::Integer)
  # special cases
  if k<0 || n<0
    throw(DomainError())
  end
  if k>n
    return big(0)
  end
  if n==0  # and, by logic, k==0
    return big(1)
  end
  if k==0  # and, by logic, n>0
    return big(0)
  end
  # end of special cases, invoke recursion
  return unsigned_Stirling1(n-1,k-1) + (n-1)*unsigned_Stirling1(n-1,k)
end

"""
  has_reasonable_precision(arb_num)

Check whether arb_num has at least roughly the double precision (maybe replace by 64 by 53).

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk_rec(6, 5, 0.5)
0.234375
```
"""
function has_reasonable_precision(arb_num)
    return accuracy_bits(arb_num) ≥ 53
end

import Nemo.risingfac
"""
  risingfac(r, n)

Computes the rising factorial.

# Examples
```julia-repl
julia> GibbsTypePriors.Cnk_rec(6, 5, 0.5)
0.234375
```
"""
function risingfac(r, n)
  if r == 0
    return 0
  else
    return prod(r + i for i in 0:(n-1))
  end
end
function risingfac(r::arb, n)
  if r == 0
    return arb_0
  else
    return prod(r + i for i in 0:(n-1))
  end
end
