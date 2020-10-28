@test_nowarn GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2)
@test GibbsTypePriors.log_βnk(1.3, 100, 25, 1.2) ≈ 2.1851379297320825
@test Float64(GibbsTypePriors.logxk(10, 5, 1.3, 0.9))  ≈ 0.7300787963127003
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(90, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(100, 100, 1.2, 0.6)
@test GibbsTypePriors.Pkn_NGG_approx(100, 100, 1.2, 0.6, 100, -1) == -1
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(100, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(1, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(1, 1000, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(20, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(20, 1.2, 0.6; verbose = true)
@test_nowarn GibbsTypePriors.Pkn_NGG(1, 20, 1.2, 0.6)

@test_nowarn GibbsTypePriors.Vnk_2PD(7, 5, 0.2, 0.4)

@test_nowarn GibbsTypePriors.Pkn_PY(7, 5, 0.2, 0.4)
@test_nowarn GibbsTypePriors.Pkn_Dirichlet(7, 10, 1.5)

@test_nowarn GibbsTypePriors.Pkn_Dirichlet(7, 10, 1.5)

# c=0.1, k=0.1, sigma= 0.7
Pkn_NGG(n, β, σ)
Pkn_NGG_pred_approx(n, β, σ)

Pkn_Dirichlet(7, 10, 1.5)
Pkn_Dirichlet_Mult(k, n, H, θ)


Pkn_PYM(k, n, H, θ, σ)


Pkn_NGGM_direct(k, n, H, β, σ)
Pkn_NGGM_precomp(k,n, H, β, σ, Pk_NGG)

# c=0.1, k=0.1, sigma= 0.7
#beta = a*k^(sigma)/sigma
##comparison finite-dim nrmi
pk = Pkn_NGG_mult(100, 50, 0.028503747356698292, 0.7)
plot(1:40, pk)

pk= Pkn_Dirichlet_Mult.(1:100, 100, 50, 11.0)
plot(1:40, pk[1:40])


pk_50 = Pkn_NGG_mult(100, 50, 0.028503747356698292, 0.7)
pk_250 = Pkn_NGG_mult(100, 250, 0.028503747356698292, 0.7)
pk = Pkn_NGG(100, 0.028503747356698292, 0.7)
plot(1:40, pk_50[1:40])
plot!(1:40, pk_250[1:40])
plot!(1:40, pk[1:40])


pk_50= Pkn_Dirichlet_Mult.(1:100, 100, 50, 5.87)
pk_250= Pkn_Dirichlet_Mult.(1:100, 100, 250, 5.87)
pk_DP= Pkn_Dirichlet.(1:100, 100, 5.87)
plot(1:40, pk_50[1:40])
plot!(1:40, pk_250[1:40])
plot!(1:40, pk_DP[1:40])




#c=0.5, k=1, sigma= 0.5
#beta = a*k^(sigma)/sigma
#
pk_50 = Pkn_NGG_mult(100, 50,1.0, 0.5)
pk_250 = Pkn_NGG_mult(100, 250,1.0, 0.5)
pk = Pkn_NGG(100,1.0, 0.5)
pk_FK = prior_Kn(0.5, 1.0, 0.5, 100, 5000)
plot(1:40, pk_50[1:40])
plot!(1:40, pk_250[1:40])
plot!(1:40, pk_FK[1].p_k[1:40]) ##not validated
plot!(1:40, pk[1:40])
