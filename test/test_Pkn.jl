@test_nowarn GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2)
@test GibbsTypePriors.log_βnk(1.3, 100, 25, 1.2) ≈ 2.1851379297320825
@test Float64(GibbsTypePriors.logxk(10, 5, 1.3, 0.9))  ≈ 0.6969596408073483
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(90, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(100, 100, 1.2, 0.6)
@test GibbsTypePriors.Pkn_NGG_approx(100, 100, 1.2, 0.6, 100, -1) == -1
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(100, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(1, 100, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(1, 1000, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_approx(20, 1.2, 0.6)
@test_nowarn GibbsTypePriors.Pkn_NGG_robust(20, 1.2, 0.6; verbose = true)
@test_nowarn GibbsTypePriors.Pkn_NGG(1, 20, 1.2, 0.6)

@test_nowarn GibbsTypePriors.Vnk_2PD(7, 5, 0.2, 0.4)

@test_nowarn GibbsTypePriors.Pkn_PY(5, 0.2, 0.4)
@test_nowarn GibbsTypePriors.Pkn_Dirichlet(7, 10, 1.5)

@test GibbsTypePriors.log_βnk(1.3, 100, 25, 1.2) ≈ 2.1851379297320825
@test Float64(GibbsTypePriors.log_βnk(1., 10, 1, 1.2))== log(10)
@test GibbsTypePriors.log_βnk(1., 10, 2, 0.5) ≈ (log(10) -2*log(2))
@test GibbsTypePriors.βnk(1., 10, 2, 0.5) ≈ 10/4
@test_nowarn Float64(GibbsTypePriors.logxk(10., 1,1.0, 0.5)) ==  Float64(log(0.5*1 + GibbsTypePriors.βnk(1.0, 10, 1, 0.5)) + log(GibbsTypePriors.Cnk(11,2,0.5)) - log(0.5) - log(GibbsTypePriors.Cnk(11,1,0.5)))
@test_nowarn GibbsTypePriors.logxk_py(10, 1, 1.0, 0.1)
@test_nowarn Float64(GibbsTypePriors.logxk_py(10., 1,1.0, 0.5)) ==  Float64(log(0.5*1 + 1.0) + log(GibbsTypePriors.Cnk(11,2,0.5)) - log(0.5) - log(GibbsTypePriors.Cnk(11,1,0.5)))


@testset "NGG prior Kn distribution" begin
    @test_nowarn GibbsTypePriors.Pkn_NGG_arb(5, 5, 0.2, 0.4)
    @test_nowarn GibbsTypePriors.Pkn_NGG_robust_in(7, 5, 0.2, 0.4)
    @test_nowarn GibbsTypePriors.Pkn_NGG_pred_approx(100, 1.0, 0.5)
    @test  Float64(GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2)) ≈  Float64(GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2))
    #@test Nemo.isequal(GibbsTypePriors.Pkn_NGG_robust(1, 100, 1.2, 0.5),GibbsTypePriors.Pkn_NGG_robust_in(1, 100, 1.2, 0.5))
    #@test Nemo.isequal(GibbsTypePriors.Pkn_NGG_robust_in(50, 100, 0.5, 0.2),GibbsTypePriors.Pkn_NGG_arb(50, 100, 0.5, 0.2))
    @test  Float64(sum(GibbsTypePriors.Pkn_NGG_robust(100, 0.5, 0.2))) ≈  1.0
    @test  Float64(sum(GibbsTypePriors.Pkn_NGG(100, 0.5, 0.2))) ≈  1.0
    @test  Float64(sum(GibbsTypePriors.Pkn_NGG_approx(100, 0.5, 0.2))) ≈  1.0
    @test  GibbsTypePriors.Pkn_NGG_pred_approx(10, 1.0, 0.5) == convert(Array{Float64,1},GibbsTypePriors.Pkn_NGG_approx(10, 1.0, 0.5))
    @test  GibbsTypePriors.Pkn_NGG_pred_approx(50, 10.0, 0.01) == convert(Array{Float64,1},GibbsTypePriors.Pkn_NGG_approx(50, 10.0, 0.01))
    @test_nowarn GibbsTypePriors.Pkn_PY_pred_approx(100, 1.0, 0.1)
    @test  Float64(sum(GibbsTypePriors.Pkn_PY_pred_approx(100, 1.0, 0.1))) ≈  1.0
    @test  Float64(sum(GibbsTypePriors.Pkn_NGG_pred_approx(50, 10.0, 0.01))) ≈  1.0
    @test  Float64(sum(GibbsTypePriors.Pkn_NGG_approx(10, 1.0, 0.5))) ≈  1.0
    @testset "NGG  multinomial" begin
        @test_nowarn GibbsTypePriors.Pkn_NGGM_arb(1,100,100, 1.0, 0.5)
        @test_nowarn GibbsTypePriors.Pkn_NGGM_arb(99,100,100, 100.0, 0.5)
        @test_nowarn GibbsTypePriors.Pkn_NGGM(50,100,100, 100.0, 0.5)
        @test  Float64(sum(GibbsTypePriors.Pkn_NGGM_arb.(1:10,10,10, 100.0, 0.5))) ≈  1.0
        @test  sum(GibbsTypePriors.Pkn_NGGM.(1:10,10,10, 100.0, 0.5)) ==  1.0
    end
end


@testset "PY approximation for prior Kn distribution" begin
    @test_nowarn GibbsTypePriors.Pkn_PY(50, 0.2, 0.4)
    @test_nowarn GibbsTypePriors.Pkn_PY( 500, 0.2, 0.4)
    @test  Float64(sum(GibbsTypePriors.Pkn_PY( 100, 0.2, 0.4))) ≈ 1.0
    @test_nowarn GibbsTypePriors.Pkn_PY_pred_approx(100, 1.0, 0.1)
    @test_nowarn  GibbsTypePriors.Pkn_2PD_arb(1,100, 0.2, 0.4)
    @test Float64( GibbsTypePriors.Pkn_2PD_arb(1,100, 0.2, 0.4)) ==  GibbsTypePriors.Pkn_2PD(1,100, 0.2, 0.4)
    @test Float64( sum(GibbsTypePriors.Pkn_2PD_arb.(1:100,100, 0.2, 0.4))) ≈  1.
    @test Float64( sum(GibbsTypePriors.Pkn_2PD.(1:100,100, 0.2, 0.4))) ≈  1.
    @test Float64( sum(GibbsTypePriors.Pkn_PY(100, 0.2, 0.4))) ≈  1.
    @test GibbsTypePriors.Pkn_PY(100, 0.2, 0.4)== convert(Array{Float64,1},GibbsTypePriors.Pkn_2PD_arb.(1:100,100, 0.2, 0.4))
    @test  Float64(sum(GibbsTypePriors.Pkn_PY_pred_approx(100, 1.0, 0.1))) ≈  1.
    @test_nowarn GibbsTypePriors.Pkn_PDM_arb(1,100,100,1.0, 0.1)
    @test_nowarn GibbsTypePriors.Pkn_PDM_arb(1,100,200,1.0, 0.1)
    @test_nowarn GibbsTypePriors.Pkn_PDM_arb(1,100,200,100.0, 0.9)
    @test Float64(sum(GibbsTypePriors.Pkn_PDM_arb.(1:100,100,200,100.0, 0.9))) ==1
    @test_nowarn GibbsTypePriors.Pkn_PYM(1,100,200,100.0, 0.9)
    @test sum(GibbsTypePriors.Pkn_PYM.(1:100,100,200,100.0, 0.9)) ≈  1
end



@testset "DP approximation for prior Kn distribution" begin
    @test_nowarn GibbsTypePriors.Pkn_Dirichlet(7, 50, 0.2)
    @test_nowarn GibbsTypePriors.Pkn_Dirichlet_arb(7, 50, 0.2)
    @test  Float64( GibbsTypePriors.Pkn_Dirichlet_arb(7, 50, 0.2)) ==  GibbsTypePriors.Pkn_Dirichlet(7, 50, 0.2)
    @test  convert(Array{Float64,1}, GibbsTypePriors.Pkn_Dirichlet_arb.(1:50, 50, 0.2)) ==  GibbsTypePriors.Pkn_Dirichlet.(1:50, 50, 0.2)
    @test_nowarn GibbsTypePriors.Pkn_Dirichlet_mult_arb(7, 50,50,  0.2)
    @test GibbsTypePriors.Pkn_Dirichlet_mult_arb(60, 50,100,  0.2) == 0
    @test  sum(GibbsTypePriors.Pkn_Dirichlet_mult.(1:100,100, 100, 0.4)) ≈ 1.0
    @test  GibbsTypePriors.Pkn_Dirichlet_mult(50,100, 100, 0.4) <  GibbsTypePriors.Pkn_Dirichlet_mult(50,100, 400, 0.4)
    @test  GibbsTypePriors.Pkn_Dirichlet_mult(50,100, 100, 0.4) <  GibbsTypePriors.Pkn_Dirichlet(50,100, 0.4)
end
