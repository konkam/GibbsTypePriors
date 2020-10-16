using Nemo, GibbsTypePriors, RCall, JLD

grid_k = 1:100

n, β, σ = 100, 1.2, 0.6


exact_Pkn = GibbsTypePriors.Pkn_NGG_arb.(grid_k, n, β, σ)
approx_Pkn = GibbsTypePriors.Pkn_NGG_approx.(grid_k, n, β, σ, n, exact_Pkn[n])

R"library(tidyverse)"
R"p2 = tibble(k = $grid_k,
         exact_Pkn = $(Float64.(exact_Pkn)),
         approx_Pkn = $(Float64.(approx_Pkn))) %>%
         #mutate(abs_rel_error = abs(approx_Pkn-exact_Pkn)/exact_Pkn) %>%
         mutate(Error = approx_Pkn-exact_Pkn) %>%
         #ggplot(aes(x = k, y = abs_rel_error)) +
         ggplot(aes(x = k, y = Error)) +
         theme_bw() +
         geom_point() +
         geom_line()"

R"p1 = tibble(k = $grid_k,
        exact_Pkn = $(Float64.(exact_Pkn)),
        approx_Pkn = $(Float64.(approx_Pkn))) %>%
        gather(type, value, -k) %>%
        mutate(type = gsub('_Pkn', '', type)) %>%
        ggplot(aes(x = k, y = value, colour = type, group = type)) +
        theme_bw() +
        geom_point() +
        geom_line() +
        theme(legend.position = 'left') +
        scale_colour_discrete(name='') +
        ylab('Pkn')"

R"p = gridExtra::grid.arrange(p1, p2, nrow = 1)"

R"pdf('test/graphical_tests/figures_graphical_tests/Pkn_NGG_approx_quality.pdf', height = 5, width = 7)
        plot(p)
dev.off()"
