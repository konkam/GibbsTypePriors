@test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60.)
@test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60., 50)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4, 20)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_stable(10, 0.4, 20)
@test_nowarn GibbsTypePriors.expected_number_of_cluster_stable(10, 0.4)
@test_nowarn GibbsTypePriors.expected_number_of_clusters_PD_direct(100, 10., 0.1)
@test_nowarn GibbsTypePriors.expected_number_of_clusters_PD_direct(100, 0.01, 0.9)
@test_nowarn GibbsTypePriors.variance_number_of_clusters_PD_direct(100, 0.01, 0.9)
@test_nowarn GibbsTypePriors.variance_number_of_clusters_PD_direct(1000, 0.01, 0.1)


@test_nowarn GibbsTypePriors.E_Dirichlet(500, 60.)
@test_nowarn GibbsTypePriors.E_2PD(10, 0., 0.4)
@test_nowarn GibbsTypePriors.E_PY(10, 0., 0.4)
@test_nowarn GibbsTypePriors.E_stable(10, 0.4)



@testset "Expectation DP" begin
    @test_nowarn GibbsTypePriors.E_Dirichlet(500, 60.)
    @test_nowarn GibbsTypePriors.E_Dirichlet(100, 0.01)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60.)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet_direct(1000, 60.)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet_direct(100, 0.01)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial_direct(100, 100, 0.01)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial_direct(1000, 1000, 0.01)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial(100, 100, 0.01)
        @testset "Comparison between methods" begin
        @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60., 150) < @test_nowarn GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60.)
        @test Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet(500, 60.)) ≈  Float64(GibbsTypePriors.E_Dirichlet(500, 60.))
        @test Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet_direct(100, 6.)) ≈  Float64(GibbsTypePriors.E_Dirichlet(100, 6.))
        @test Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial_direct(500,500, 6.)) ==  Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial(500,500,  6.))
        @test Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial_direct(500,100, 6.)) <  Float64(GibbsTypePriors.E_Dirichlet(500,  6.))
        @test Float64(GibbsTypePriors.expected_number_of_clusters_Dirichlet_Multinomial(500,300, 6.)) <  Float64(GibbsTypePriors.E_Dirichlet(500,  6.))
    end
end


@testset "Expectation PY" begin
    @test_nowarn GibbsTypePriors.E_PY(500, 60., 0.1)
    @test_nowarn GibbsTypePriors.E_PY(100, 60., 0.9)
    @test_nowarn GibbsTypePriors.expected_number_of_clusters_PD_direct(100, 60., 0.9)
    @test_nowarn GibbsTypePriors.E_PY(500, 0.1, 0.99)
    @testset "Comparison between methods" begin
        @test Float64(GibbsTypePriors.E_PY(100, 60., 0.1)) ≈  Float64(GibbsTypePriors.expected_number_of_clusters_PD_direct(100, 60., 0.1))
        @test Float64(GibbsTypePriors.expected_number_of_cluster_2PD(10, 0., 0.4)) == Float64(GibbsTypePriors.E_PY(10, 0., 0.4))
        @test Float64(GibbsTypePriors.E_PY(100, 6., 0.9)) ≈  Float64(GibbsTypePriors.expected_number_of_clusters_PD_direct(100, 6., 0.9))
        @test Float64(GibbsTypePriors.E_2PD(100, 6., 0.9, 50)) < Float64(GibbsTypePriors.E_PY(100, 6.0, 0.9))
    end
end


@testset "Variance for Gibbs-type priors" begin
    @testset "Variance DP" begin
        @test_nowarn GibbsTypePriors.V_Dirichlet(100, 60.)
        @test_nowarn GibbsTypePriors.variance_number_of_clusters_Dirichlet(100, 0.6)
        @test_nowarn GibbsTypePriors.variance_number_of_clusters_Dirichlet_direct(1000, 60.)
        @test_nowarn GibbsTypePriors.variance_number_of_clusters_Dirichlet_direct(1000, 0.001)
        @testset "Comparison" begin
            @test Float64(GibbsTypePriors.variance_number_of_clusters_Dirichlet(100, 60.)) ≈ Float64(GibbsTypePriors.V_Dirichlet(100, 60.))
            @test Float64(GibbsTypePriors.variance_number_of_clusters_Dirichlet_direct(100, 60.)) ≈ Float64(GibbsTypePriors.V_Dirichlet(100, 60.))

        end
    end
    @testset "Variance PY" begin
        @test_nowarn GibbsTypePriors.V_PY(500, 60., 0.1)
        @test_nowarn GibbsTypePriors.V_PY(100, 60., 0.9)
        @test_nowarn GibbsTypePriors.V_PY(1000, 60., 0.9)
        @test_nowarn GibbsTypePriors.V_PY(100, 60., 0.5)
        @test_nowarn GibbsTypePriors.variance_number_of_clusters_PD(20, 5., 0.9)
        @test_nowarn GibbsTypePriors.variance_number_of_clusters_PD_direct(20, 5., 0.9)
        @testset "Comparison between methods" begin
            Float64(GibbsTypePriors.variance_number_of_clusters_PD(100, 0.5, 0.5)) ≈ Float64(GibbsTypePriors.V_PY(100, 0.5, 0.5))
        end
    end
end
