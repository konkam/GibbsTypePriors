## SB Pkn_NGG
R"library(tidyverse)
library(copula)"


using Distributions

function R_q(q,b,k_n, alpha)
    t1 = exp(- ((b/k_n)^alpha)*((1 + q^(1/alpha))^alpha - q))
    t2 = (1 + q^(1/alpha))^(alpha -1)
    return t1*t2
end

function sample_xi_n(alpha, b, theta, k_n, n)
    while r_n
end

function approximation_prior_distribution(beta,sigma,N,Nt)
    #i = findall(x->x==sigma, sigma_arr)
    Pkn_numeric_ = Pkn_NGG_numeric.(1:N, N, beta, sigma)
    Pkn_order2_ = Pkn_NGG_pred_approx(N, beta, sigma)
    Pkn_NGGMult = Pkn_NGGM_precomp.(1:N,N,Nt,beta,sigma, [Pkn_numeric_])
    #kappa = (beta*1*sigma)^(1/sigma)
    alpha = (beta*sigma)
    Pkn_NGG_FK_ = prior_Kn(alpha,1.0, sigma, N,3000; runs=10^4)[1].p_k
