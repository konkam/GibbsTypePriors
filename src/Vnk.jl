function Vnk_NGG(n, k, β, σ)
    β_arb = RR(β)
    σ_arb = RR(σ)
    return exp(β_arb) * σ_arb^(k-1) // RR(gamma(n)) * sum([binom(n-1, i, RR) * (-1)^i * β_arb ^(i//σ_arb) * real(gamma(k-i//CC(σ_arb), CC(β_arb))) for i in 0:(n-1)])
end
