@test Nemo.accuracy_bits(GibbsTypePriors.Cnk(6, 5, 0.5, 200)) ≤ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Cnk(6, 5, 0.5, 2000)) ≥ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Cnk_rec(600, 5, 0.5, 200)) ≤ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Cnk_rec(600, 5, 0.5, 2000)) ≥ 200