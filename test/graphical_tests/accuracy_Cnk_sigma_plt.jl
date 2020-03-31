using RCall, JLD, DataFrames, DataFramesMeta, GibbsTypePriors, Nemo

accuracy_Cnk_sigma = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_sigma.jld", "accuracy_Cnk_sigma") |>
    df -> @transform(df, Accuracy = min.(:Accuracy, GibbsTypePriors.RR(22)/GibbsTypePriors.RR(7.) |> accuracy_bits))

accuracy_Cnk_rec_sigma = load("test/graphical_tests/saves_for_graphical_tests/accuracy_Cnk_rec_sigma.jld", "accuracy_Cnk_sigma") |>
        df -> @transform(df, Accuracy = min.(:Accuracy, GibbsTypePriors.RR(22)/GibbsTypePriors.RR(7.) |> accuracy_bits))  |>
        df -> @transform(df, s = Float64.(:s))


R"library(tidyverse)"

R"p = as_tibble($accuracy_Cnk_sigma) %>%
    mutate(type = ifelse(test = n == k, yes = 'Large k', no = 'Small k'), fun = 'Iterative') %>%
    bind_rows(as_tibble($accuracy_Cnk_rec_sigma) %>%
        mutate(type = ifelse(test = n == k, yes = 'Large k', no = 'Small k'), fun = 'Recursive')) %>%
    ggplot(aes(x = n, y = Accuracy, colour = s, group = interaction(s, type))) +
    theme_bw() +
    facet_wrap(~fun) +
    geom_point(aes(shape = type)) +
    geom_line(aes(linetype = type)) +
    scale_linetype(name='') +
    scale_shape(name='') +
    viridis::scale_colour_viridis(name=expression(sigma)) +
    geom_hline(yintercept = 64, linetype = 'dotted', colour = 'red')
    "

R"pdf('test/graphical_tests/figures_graphical_tests/accuracy_Cnk_sigma.pdf', height = 5, width = 5)
    plot(p)
dev.off()"
