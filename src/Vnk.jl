<<<<<<< HEAD
function Vnk_NGG(n, k, β, σ; prec = 5000)
    RF = RealField(prec)
    CF = ComplexField(prec)
    β_arb = RF(β)
    σ_arb = RF(σ)
    return exp(β_arb) * σ_arb^(k-1) // RF(gamma(n)) * sum([binom(n-1, i, RF) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CF(σ_arb), CF(β_arb))) for i in 0:(n-1)])
end

function Vnk_2PD(n, k, θ, σ; prec = 5000)
    RF = RealField(prec)
    θ_arb = RF(θ)
    σ_arb = RF(σ)
=======
# function Vnk_NGG_old(n, k, β, σ)
#     β_arb = RR(β)
#     σ_arb = RR(σ)
#     return exp(β_arb) * σ_arb^(k-1) // RR(gamma(n)) * sum([binom(n-1, i, RR) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CC(σ_arb), CC(β_arb))) for i in 0:(n-1)])
# end
function Vnk_NGG(n, k, β, σ)
    β_arb = RR(β)
    σ_arb = RR(σ)
    n_m_1_Flint = FlintZZ(n-1)
    return exp(β_arb) * σ_arb^(k-1) // RR(gamma(n)) * sum(binomial(n_m_1_Flint, FlintZZ(i)) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CC(σ_arb), CC(β_arb))) for i in 0:(n-1))
end
#
# function Vnk_NGG2(n, k, β, σ)
#     gamma_n::arb = gamma(n)
#     return exp(β) * σ^(k-1) // gamma_n * sum([binom(n-1, i, RR) * (-1)^i * β^(i//σ) * real(gamma(k-i//CC(σ), CC(β))) for i in 0:(n-1)])
# end
# function Vnk_NGG3(n, k, β, σ)
#     gamma_n::arb = gamma(n)
#     return exp(β) * σ^(k-1) // gamma_n * sum(binom(n-1, i, RR) * (-1)^i * β^(i//σ) * real(gamma(k-i//CC(σ), CC(β))) for i in 0:(n-1))
# end
# function Vnk_NGG4(n, k, β, σ)
#     gamma_n::arb = gamma(n)
#     n_m_1_Flint = FlintZZ(n-1)
#     return exp(β) * σ^(k-1) // gamma_n * sum(binomial(n_m_1_Flint, FlintZZ(i)) * (-1)^i * β^(i//σ) * real(gamma(k-i//CC(σ), CC(β))) for i in 0:(n-1))
# end
function Vnk_2PD(n, k, θ, σ)
    θ_arb::arb = RR(θ)
    σ_arb::arb = RR(σ)
>>>>>>> master
    if k==1
        num::arb = arb_1
    else
        num = prod([θ_arb + i * σ_arb for i in 1:(k-1)])
    end
    return num // risingfac(1+θ_arb,n-1)
end
