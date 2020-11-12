

n, β, σ, M = 100, 10, 0.5, 200
using GibbsTypePriors
Pkn_FK = Pkn_NGG_FK(n, β, σ, M; runs=2*10^2)
Pkn_FK_fast = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)

Pkn_exact = Pkn_NGG(n, β, σ)

to_plot = DataFrame(k = 1:n, Pkn_FK = Pkn_FK, Pkn_FK_fast = Pkn_FK_fast, Pkn_exact = Pkn_exact)

R"$to_plot %>%
    as_tibble %>%
    gather(type, p_k, Pkn_FKkn_exact) %>%
    ggplot(aes(x = k, y = p_k, colour = type)) +
    theme_minimal() +
    geom_point() +
    geom_line()"
