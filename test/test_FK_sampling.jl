##NGG_FK test

include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_FK_sampling.jl")


## test different parametrizations
p_ngg_5 = TruncatedNGG(5.0, 1.0, 0.2, 5000,0.01)
p_ngg_1 = TruncatedNGG(1.0, ((5.0)^0.2)/0.2, 0.2, 5000,0.01)
p_test = TruncatedNGG_test(0.25, 1.0, 0.25, 10000,0.01)


p_ngg_5 = prior_Kn(0.25, 1.0, 0.25, 100, 5000)
p_ngg_1 = prior_Kn(1.0,  0.00390625, 0.25, 100, 5000)

using Plots

plot(1:100, p_ngg_5[1].p_k)
plot!(1:100, p_ngg_1[1].p_k)

plot(1:100, p_ngg_1[1].p_k)

p_ngg_1 = prior_Kn(2.5,  1.0, 0.25, 100, 5000,runs=10^3)

p_ngg_1_ = prior_Kn(0.25,  1.0, 0.25, 100, 5000,runs=10^4)

p_ngg_1_2 = prior_Kn(0.75,  1.0, 0.75, 100, 5000,runs=10^3)
p_ngg_1_3 = prior_Kn(7.5,  1.0, 0.75, 100, 5000,runs=10^3)

plot(1:100, p_ngg_1[1].p_k)
plot(1:100, p_ngg_1_[1].p_k)

plot!(1:100, p_ngg_1_2[1].p_k)

plot!(1:100, p_ngg_1_2[1].p_k)


plot(1:100, p_ngg_1_3[1].p_k)


## Exp  ≈ 15
p_ngg_test1 = prior_Kn(0.5,  1, 0.5, 100, 5000,runs=10^3)
plot(1:100, p_ngg_test1[1].p_k)

p_ngg_test2 = prior_Kn(0.1,  1, 0.7, 100, 50,runs=10^3)
plot(1:100, p_ngg_test2[1].p_k)




sigma_array = [0.25, 0.75]
beta_array = [(1*0.25)^(1/0.25), (1*0.75)^(1/0.75)]
n = 100
H = 3000
Kn_priors = exp_K_n_prior_NGG(sigma_array,beta_array, H, n, 100)

sigma_array = [0.25, 0.75]
beta_array = [(10*0.25)^(1/0.25), (10*0.75)^(1/0.75)]
n = 1000
H = 3000
Kn_priors = exp_K_n_prior_NGG(sigma_array,beta_array, H, n, 100)


sigma_array = [0.25, 0.75]
beta_array = [(1*0.25)^(1/0.25), (1*0.75)^(1/0.75)]
n = 1000
H = 3000
#Kn_priors = exp_K_n_prior_NGG(sigma_array,beta_array, H, n, 10000)


sigma_array = [0.25, 0.75]
beta_array = [1.0, 10.0]
n = 1000
H = 3000
#Kn_priors = exp_K_n_prior_NGG(sigma_array,beta_array, H, n, 10000)
