function Vnk_NGG(n, k, β, σ)
    β_arb = RR(β)
    σ_arb = RR(σ)
    return exp(β_arb) * σ_arb^(k-1) // RR(gamma(n)) * sum([binom(n-1, i, RR_in) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CC(σ_arb), CC(β_arb))) for i in 0:(n-1)])
end

function Vnk_2PD(n, k, θ, σ)
    θ_arb = RR(θ)
    σ_arb = RR(σ)
    if k==1
        num = 1
    else
        num =  prod([θ_arb + i * σ_arb for i in 1:(k-1)])
    end
    return num // risingfac(1+θ_arb,n-1)
end
