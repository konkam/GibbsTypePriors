using Nemo, GibbsTypePriors, JLD

grid_k = [collect(1:25:10000); 10000]

Pkn_numerical_accuracy = GibbsTypePriors.Pkn_NGG_arb.(grid_k, 10000, 1.2, 0.6) |> x -> accuracy_bits.(x)
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld", "Pkn_numerical_accuracy", Pkn_numerical_accuracy)
