using Nemo, GibbsTypePriors, JLD, DataFramesMeta

include("common_functions_for_tests.jl")

grid_n = [collect(10:100:4000); 4000]
grid_s = 0.1:0.2:0.9

accuracy_Cnk_sigma = expand_grid(grid_n, grid_s) |>
    df -> rename(df, [:n, :s]) |>
    df -> vcat(@transform(df, k = :n), @transform(df, k = repeat([1], inner = length(:n)))) |>
    df -> @transform(df, Accuracy = GibbsTypePriors.Cnk.(:n, :k, :s) |> x -> accuracy_bits.(x))
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma.jld", "accuracy_Cnk_sigma", accuracy_Cnk_sigma)
