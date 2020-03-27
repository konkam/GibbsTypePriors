using Nemo, GibbsTypePriors, JLD

grid_k = [collect(1:250:10000); 10000]


Cnk_numerical_accuracy = GibbsTypePriors.Cnk.(10000, grid_k, 0.6) |> x -> accuracy_bits.(x)
save("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld", "Cnk_numerical_accuracy", Cnk_numerical_accuracy)
