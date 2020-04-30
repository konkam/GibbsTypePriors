using Nemo, GibbsTypePriors, JLD
include("common_functions_for_tests.jl")


grid_n = [collect(10:500:5000); 5000] |> unique
grid_σ = [collect(0.1:0.3:0.9); 0.9]
grid_β = [0.5, 1., 2.]
param_grid = expand_grid(grid_n, grid_σ, grid_β) |>
    df -> rename(df, [:n, :σ, :β])

Nit = size(param_grid,1)

# Vnk_numerical_accuracy = GibbsTypePriors.Vnk_NGG.(param_grid[!,:n], 1, param_grid[!,:β], param_grid[!,:σ]) |> x -> accuracy_bits.(x)

using Base.Threads
println("Computing Vn1 using $(nthreads()) threads")
Vnk_numerical_accuracy = Array{Int64}(undef, Nit)
@threads for i in 1:Nit
          Vnk_numerical_accuracy[i] = GibbsTypePriors.Vnk_NGG.(param_grid[i,:n], 1, param_grid[i,:β], param_grid[i,:σ]) |> x -> accuracy_bits.(x)
       end

save("test/graphical_tests/saves_for_graphical_tests/accuracy_Vn1.jld", "Vnk_numerical_accuracy", Vnk_numerical_accuracy)
