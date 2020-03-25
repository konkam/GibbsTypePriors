using Nemo, GibbsTypePriors, RCall

to_plot = Float64.(GibbsTypePriors.Pkn_robust(20, 1.2, 0.6; verbose = true))

R"library(tidyverse)"
R"tibble(Pkn = $to_plot) %>%
    rowid_to_column(var = 'k') %>%
    ggplot(aes(x = k, y = Pkn)) +
    theme_bw() +
    geom_point() +
    geom_line()"


    to_plot = Float64.(GibbsTypePriors.Pkn_robust(1000, 1.2, 0.6; verbose = true))

    R"library(tidyverse)"
    R"tibble(Pkn = $to_plot) %>%
        rowid_to_column(var = 'k') %>%
        ggplot(aes(x = k, y = Pkn)) +
        theme_bw() +
        geom_point() +
        geom_line()"


to_plot = GibbsTypePriors.Pkn_NGG_arb.(1:1000, 1000, 1.2, 0.6) |> x -> accuracy_bits.(x)
