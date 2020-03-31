using Nemo, GibbsTypePriors, RCall, JLD
R"library(tidyverse)"

to_plot = Float64.(GibbsTypePriors.Pkn_NGG_robust(1000, 1.2, 0.6; verbose = true))

R"p = tibble(Pkn = $to_plot) %>%
    rowid_to_column(var = 'k') %>%
    ggplot(aes(x = k, y = Pkn)) +
    theme_bw() +
    geom_point() +
    geom_line()"

R"pdf('test/graphical_tests/figures_graphical_tests/dist_Pkn_1000.pdf', height = 5, width = 5)
    plot(p)
    dev.off()"
