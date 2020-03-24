function Pkn_NGG_arb(k, n,  β, σ)
    σ_arb = RR(σ)
    Vnk_exact(n, k, β, σ) // σ_arb^k * Cnk(n, k, σ)
end
