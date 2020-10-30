@test_nowarn GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2)
@test GibbsTypePriors.log_βnk(1.3, 100, 25, 1.2) ≈ 2.1851379297320825
@test Float64(GibbsTypePriors.logxk(10, 5, 1.3, 0.9))  ≈ 0.6969596408073483
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
