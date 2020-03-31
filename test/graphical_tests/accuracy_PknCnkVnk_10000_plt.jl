using RCall, JLD

grid_k = [collect(1:25:10000); 10000]


Vnk_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_10000.jld", "Vnk_numerical_accuracy") |> x -> Float64.(x)
Pkn_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_10000.jld", "Pkn_numerical_accuracy") |> x -> Float64.(x)
Cnk_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_10000.jld", "Cnk_numerical_accuracy") |> x -> Float64.(x)

R"library(tidyverse)"
R"p = tibble(Pkn = $Pkn_numerical_accuracy, Cnk = $Cnk_numerical_accuracy, Vnk = $Vnk_numerical_accuracy, k = $(collect(grid_k)) ) %>%
    gather(Var, Accuracy, -k) %>%
    mutate(Var = factor(Var, levels = c('Vnk', 'Cnk', 'Pkn'))) %>%
    ggplot(aes(x = k, y = Accuracy, colour = Var, group = Var)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red') +
    scale_colour_discrete(name='')"

R"pdf('test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_10000.pdf', height = 5, width = 5)
    plot(p)
    dev.off()"
