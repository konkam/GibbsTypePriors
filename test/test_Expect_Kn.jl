@test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60., 500)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4, 20)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_stable(10, 0.4, 20)

@test_nowarn GibbsTypePriors.E_Dirichlet(500, 60., 500)
@test_nowarn GibbsTypePriors.E_2PD(10, 0., 0.4, 20)
@test_nowarn GibbsTypePriors.E_PY(10, 0., 0.4, 20)
@test_nowarn GibbsTypePriors.E_stable(10, 0.4, 20)
