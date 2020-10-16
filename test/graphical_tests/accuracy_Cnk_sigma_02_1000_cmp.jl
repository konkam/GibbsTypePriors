using Nemo, GibbsTypePriors, JLD

grid_k = [collect(1:25:1000); 1000]

Cnk_numerical_accuracy = GibbsTypePriors.Cnk.(1000, grid_k, 0.2) |> x -> accuracy_bits.(x)
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma_02_1000.jld", "Cnk_numerical_accuracy", Cnk_numerical_accuracy)
