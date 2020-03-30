using Nemo, GibbsTypePriors, RCall, JLD

to_plot = Float64.(GibbsTypePriors.Pkn_robust(20, 1.2, 0.6; verbose = true))

R"library(tidyverse)"
R"tibble(Pkn = $to_plot) %>%
    rowid_to_column(var = 'k') %>%
    ggplot(aes(x = k, y = Pkn)) +
    theme_bw() +
    geom_point() +
    geom_line()"


to_plot = Float64.(GibbsTypePriors.Pkn_robust(1000, 1.2, 0.6; verbose = true))

R"tibble(Pkn = $to_plot) %>%
    rowid_to_column(var = 'k') %>%
    ggplot(aes(x = k, y = Pkn)) +
    theme_bw() +
    geom_point() +
    geom_line()"


to_plot = GibbsTypePriors.Pkn_NGG_arb.(1:1000, 1000, 1.2, 0.6) |> x -> accuracy_bits.(x)

save("Pkn_numerical_accuracy.jld", "Pkn_numerical_accuracy", to_plot)


R"p = tibble(Accuracy = $to_plot) %>%
    rowid_to_column(var = 'k') %>%
    ggplot(aes(x = k, y = Accuracy)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red')"

R"pdf('Pkn_numerical_accuracy.pdf')
plot(p)
dev.off()"

Cnk_numerical_accuracy = GibbsTypePriors.Cnk.(1000, 1:1000, 0.6) |> x -> accuracy_bits.(x)

save("Cnk_numerical_accuracy.jld", "Cnk_numerical_accuracy", Cnk_numerical_accuracy)

Vnk_numerical_accuracy = GibbsTypePriors.Vnk_NGG.(1000, 1:1000, 1.2, 0.6) |> x -> accuracy_bits.(x)

save("Vnk_numerical_accuracy.jld", "Vnk_numerical_accuracy", Vnk_numerical_accuracy)

Pkn_numerical_accuracy = load("Pkn_numerical_accuracy.jld", "Pkn_numerical_accuracy")
Cnk_numerical_accuracy = load("Cnk_numerical_accuracy.jld", "Cnk_numerical_accuracy")
Vnk_numerical_accuracy = load("Vnk_numerical_accuracy.jld", "Vnk_numerical_accuracy")

R"p = tibble(Pkn = $Pkn_numerical_accuracy, Cnk = $Cnk_numerical_accuracy, Vnk = $Vnk_numerical_accuracy ) %>%
    rowid_to_column(var = 'k') %>%
    gather(Var, Accuracy, -k) %>%
    mutate(Var = factor(Var, levels = c('Vnk', 'Cnk', 'Pkn'))) %>%
    ggplot(aes(x = k, y = Accuracy, colour = Var, group = Var)) +
    theme_bw() +
    geom_point() +
    geom_line() +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red') +
    scale_colour_discrete(name='')"

R"pdf('Pkn_Cnk_Vnk_numerical_accuracy.pdf')
    plot(p)
    dev.off()"
