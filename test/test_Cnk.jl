@test_nowarn GibbsTypePriors.Cnk(6, 5, 0.5)
@test_nowarn GibbsTypePriors.Cnk(6, 5, GibbsTypePriors.RR(0.5))

@test Float64(GibbsTypePriors.Cnk(6, 5, 0.5)) == Float64(GibbsTypePriors.Cnk_rec(6, 5, 0.5))
@test Float64(GibbsTypePriors.Cnk(600, 599, 0.5)) ≈ Float64(GibbsTypePriors.Cnk_rec(600, 599, 0.5)) atol=10^(-8)
@test Float64(GibbsTypePriors.Cnk(6, 5, 0.5)) == Float64(GibbsTypePriors.Cnk_rec(6, 5, GibbsTypePriors.RR(0.5)))
@test Float64(GibbsTypePriors.Cnk(600, 599, 0.5)) ≈ Float64(GibbsTypePriors.Cnk_rec(600, 599, GibbsTypePriors.RR(0.5))) atol=10^(-8)

@test_nowarn GibbsTypePriors.Cnk_robust(6, 5, 0.5)
@test_nowarn GibbsTypePriors.Cnk_robust(1000, 999, 0.01; verbose = true)

@test_nowarn GibbsTypePriors.noncentral_generalised_factorial_coefficient(10, 5, 0.2, 0.)

@test_nowarn GibbsTypePriors.Cnk_rec(6, 5, 0.5)