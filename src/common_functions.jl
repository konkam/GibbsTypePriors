using SpecialFunctions, Nemo, StatsFuns, Memoization
import Nemo.binom, Nemo.gamma

prec = 5000 ## Increase for better precision
RR = RealField(prec)
CC = ComplexField(prec)
binom(n::Int64, k::Int64, r::Nemo.ArbField) = binom(convert(UInt64, n),convert(UInt64, k), r)
gamma(x::Int64) = Nemo.gamma(RR(x))

# const arb_1 = RR(1)
# const arb_0 = RR(0)

@memoize function unsigned_Stirling1(n::Integer,k::Integer)
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

function has_reasonable_precision(arb_num)
    return accuracy_bits(arb_num) ≥ 64
end

import Nemo.risingfac
function risingfac(r, n)
  if r == 0
    return 0
  else
    return prod(r + i for i in 0:(n-1))
  end
end
function risingfac(r::arb, n)
  if r == 0
    return 0*r #To keep type stability
  else
    return prod(r + i for i in 0:(n-1))
  end
end
