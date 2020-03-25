@test_nowarn GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2)
@test GibbsTypePriors.log_βnk(1.3, 100, 25, 1.2) ≈ 2.1851379297320825
@test Float64(GibbsTypePriors.logxk(10, 5, 1.3, 0.9))  ≈ 0.7300787963127003
@test_nowarn GibbsTypePriors.Pkn_approx(90, 100,  1.2, GibbsTypePriors.RR(0.6))
@test_nowarn GibbsTypePriors.Pkn_approx(100, 100,  1.2, GibbsTypePriors.RR(0.6))
@test_nowarn GibbsTypePriors.Pkn_robust(100, 100,  1.2, GibbsTypePriors.RR(0.6))
@test_nowarn GibbsTypePriors.Pkn_robust(1, 100,  1.2, GibbsTypePriors.RR(0.6))
@test_nowarn GibbsTypePriors.Pkn_robust(1, 1000,  1.2, GibbsTypePriors.RR(0.6))
@test_nowarn GibbsTypePriors.Pkn_approx(20, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_robust(20, 1.2, 0.6; verbose = true)
