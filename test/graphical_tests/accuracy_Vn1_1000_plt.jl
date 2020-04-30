using RCall, JLD
include("common_functions_for_tests.jl")


grid_n = [collect(10:100:1000); 1000] |> unique
grid_σ = [collect(0.1:0.3:0.9); 0.9]
grid_β = [0.5, 1., 2.]
param_grid = expand_grid(grid_n, grid_σ, grid_β) |>
    df -> rename(df, [:n, :σ, :β])

Nit = size(param_grid,1)

# Vnk_numerical_accuracy = GibbsTypePriors.Vnk_NGG.(param_grid[!,:n], 1, param_grid[!,:β], param_grid[!,:σ]) |> x -> accuracy_bits.(x)

Vnk_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Vn1.jld", "Vnk_numerical_accuracy")
R"library(tidyverse)"
R"p = tibble(Vnk = $Vnk_numerical_accuracy, n = $(param_grid[!,:n]), s = $(param_grid[!,:σ]), b = $(param_grid[!,:β])) %>%
    ggplot(aes(x = n, y = Vnk, colour = s, group = s)) +
    theme_bw() +
    facet_wrap(~b) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red') +
    viridis::scale_colour_viridis()"

R"pdf('test/graphical_tests/figures_graphical_tests/accuracy_Vn1.pdf', height = 5, width = 5)
    plot(p)
dev.off()"
