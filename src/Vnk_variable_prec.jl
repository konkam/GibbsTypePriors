function Vnk_NGG(n, k, β, σ, prec)
    RF = RealField(prec)
    CF = ComplexField(prec)
    β_arb = RF(β)
    σ_arb = RF(σ)
    return exp(β_arb) * σ_arb^(k-1) // RF(gamma(n)) * sum([binom(n-1, i, RF) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CF(σ_arb), CF(β_arb))) for i in 0:(n-1)])
end

function Vnk_2PD(n, k, θ, σ, prec)
    RF = RealField(prec)
    θ_arb = RF(θ)
    σ_arb = RF(σ)
    if k==1
        num::arb = RF(1)
    else
        num = prod([θ_arb + i * σ_arb for i in 1:(k-1)])
    end
    return num // risingfac(1+θ_arb,n-1)
end