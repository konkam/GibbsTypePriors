using Nemo, GibbsTypePriors, JLD

grid_k = [collect(1:25:1000); 1000]

Pkn_numerical_accuracy = GibbsTypePriors.Pkn_NGG_arb.(grid_k, 1000, 1.2, 0.6) |> x -> accuracy_bits.(x)
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld", "Pkn_numerical_accuracy", Pkn_numerical_accuracy)
