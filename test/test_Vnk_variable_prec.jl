@test Nemo.accuracy_bits(GibbsTypePriors.Vnk_NGG(950, 50, 0.5, 0.2, 200)) ≤ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Vnk_NGG(950, 50, 0.5, 0.2, 2000)) ≥ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Vnk_2PD(7, 6, 0.5, 0.01, 200)) ≤ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Vnk_2PD(7, 6, 0.5, 0.01, 2000)) ≥ 200