@test_nowarn GibbsTypePriors.Cnk(6, 5, 0.5)
@test_nowarn GibbsTypePriors.Cnk(6, 5, GibbsTypePriors.RR(0.5))

@test Float64(GibbsTypePriors.Cnk(6, 5, 0.5)) == Float64(GibbsTypePriors.Cnk_rec(6, 5, 0.5))
@test Float64(GibbsTypePriors.Cnk(600, 599, 0.5)) == Float64(GibbsTypePriors.Cnk_rec(600, 599, 0.5))
