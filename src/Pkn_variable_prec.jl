function Pkn_NGG_arb(k, n, β, σ, prec)
    RF = RealField(prec)
    if k>n || k==0
        return RF(0)
    else
        σ_arb = RF(σ)
        Vnk_NGG(n, k, β, σ, prec) // σ_arb^k * Cnk(n, k, σ, prec)
    end
end