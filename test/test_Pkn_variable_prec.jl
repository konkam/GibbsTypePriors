@test Nemo.accuracy_bits(GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2, 200)) ≤ 200
@test Nemo.accuracy_bits(GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2, 2000)) ≥ 200