using RCall, JLD

grid_k = [collect(1:25:1000); 1000]


Vnk_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_1000.jld", "Vnk_numerical_accuracy") |> x -> Float64.(x)
Pkn_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_1000.jld", "Pkn_numerical_accuracy") |> x -> Float64.(x)
Cnk_numerical_accuracy = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_1000.jld", "Cnk_numerical_accuracy") |> x -> Float64.(x)

Vnk_numerical_accuracy_sigma_02 = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Vnk_sigma_02_1000.jld", "Vnk_numerical_accuracy") |> x -> Float64.(x)
Pkn_numerical_accuracy_sigma_02 = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Pkn_sigma_02_1000.jld", "Pkn_numerical_accuracy") |> x -> Float64.(x)
Cnk_numerical_accuracy_sigma_02 = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma_02_1000.jld", "Cnk_numerical_accuracy") |> x -> Float64.(x)

R"library(tidyverse)"
R"p = tibble(Pkn = $Pkn_numerical_accuracy, Cnk = $Cnk_numerical_accuracy, Vnk = $Vnk_numerical_accuracy, k = $(collect(grid_k)), s = 0.6 ) %>%
    bind_rows(tibble(Pkn = $Pkn_numerical_accuracy_sigma_02, Cnk = $Cnk_numerical_accuracy_sigma_02, Vnk = $Vnk_numerical_accuracy_sigma_02, k = $(collect(grid_k)), s = 0.2 )) %>%
    gather(Var, Accuracy, -(k:s)) %>%
    mutate(Var = factor(Var, levels = c('Vnk', 'Cnk', 'Pkn'))) %>%
    ggplot(aes(x = k, y = Accuracy, colour = Var, group = Var)) +
    theme_bw() +
    facet_wrap(~s, labeller = label_both) +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red') +
    scale_colour_discrete(name='')"

R"pdf('test/graphical_tests/figures_graphical_tests/accuracy_PknCnkVnk_1000.pdf', height = 5, width = 5)
    plot(p)
    dev.off()"
