@test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60.)
@test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60., 50)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4, 20)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_stable(10, 0.4, 20)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_stable(10, 0.4)

@test_nowarn GibbsTypePriors.E_Dirichlet(500, 60.)
@test_nowarn GibbsTypePriors.E_2PD(10, 0., 0.4)
@test_nowarn GibbsTypePriors.E_PY(10, 0., 0.4)
@test_nowarn GibbsTypePriors.E_stable(10, 0.4)

## Dirichlet test
E_Dirichlet_multinomial(100, 50,5.87)
E_Dirichlet_multinomial(100, 200,5.87)
E_Dirichlet(100, 5.87)


E_Dirichlet_multinomial(100, 50,11.)
E_Dirichlet_multinomial(100, 200,11.)
E_Dirichlet(100, 11.)

## PY test
convert(Float64,E_PY(100, 1., 0.5))  ≈ convert(Float64,E_2PD(100, 1.0, 0.5))
convert(Float64,E_PY(100, 1., 0.5)) ≈ convert(Float64, E_PY_exact(100, 1.0, 0.5))
V_PY_exact((n::N, θ::T, σ::T))

## NGG test

E_NGG(100, 0.028503747356698292, 0.7)
E_NGG_multinomial(100,100, 0.028503747356698292, 0.7)
E_NGG_multinomial(100,200, 0.028503747356698292, 0.7)



E_NGG(100, 1.0, 0.5)
E_NGG_multinomial(100,100, 1.0, 0.5)
E_NGG_multinomial(100,200, 1.0, 0.5)
