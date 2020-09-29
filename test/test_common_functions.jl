@test 1==1

@test_nowarn GibbsTypePriors.binom(6, 5, GibbsTypePriors.RR)
@test_nowarn GibbsTypePriors.gamma(10)
@test_nowarn GibbsTypePriors.unsigned_Stirling1(10, 5)
@test GibbsTypePriors.has_reasonable_precision(GibbsTypePriors.RR(6)) == true
@test_nowarn GibbsTypePriors.risingfac(0.5, 10)
@test_nowarn GibbsTypePriors.risingfac(GibbsTypePriors.RR(.6), 10)