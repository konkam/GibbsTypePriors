using SpecialFunctions, Nemo, StatsFuns, Memoization
import Nemo.binom, Nemo.gamma

prec = 5000 ## Increase for better precision
RR = RealField(prec)
CC = ComplexField(prec)
binom(n::Int64, k::Int64, r::Nemo.ArbField) = binom(convert(UInt64, n),convert(UInt64, k), r)
gamma(x::Int64) = Nemo.gamma(RR(x))

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
