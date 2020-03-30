using Nemo, GibbsTypePriors, JLD

grid_k = [collect(1:25:1000); 1000]

Vnk_numerical_accuracy = GibbsTypePriors.Vnk_NGG.(1000, grid_k, 1.2, 0.6) |> x -> accuracy_bits.(x)
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld", "Vnk_numerical_accuracy", Vnk_numerical_accuracy)
