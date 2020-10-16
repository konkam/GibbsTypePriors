##NGG_FK test

include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_FK_sampling.jl")



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
