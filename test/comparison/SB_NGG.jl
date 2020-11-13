## SB Pkn_NGG

include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/MvInv.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_FK_sampling.jl")

using DataFrames, DataFramesMeta, RCall
R"library(tidyverse)
library(latex2exp)
"
using Distributions

R"library(tidyverse)
library(copula)"


function R_q(q, b, k_n, alpha)
    t1 = exp(- ((b/k_n)^alpha)*((1 + q^(1/alpha))^alpha - q))
    t2 =  ((q^(1/alpha))/(1+q^(1/alpha)))^(1-alpha)
    return t1*t2
end

function sample_xi_n(alpha, b, theta, k_n, n)
   u =  Uniform(0,1)
   g = Gamma(theta/alpha +n, 1/((b/k_n)^(alpha)))
   r_n = rand(u,1)[1]
   dz_n  = rand(g,1)[1]
   while r_n > R_q(dz_n, b, k_n,alpha)
       u =  Uniform(0,1)
       r_n = rand(u,1)[1]
       dz_n  = rand(g,1)[1]
       if r_n<= dz_n
        # println("true")
       end
         #return (b/k_n)*((dz_n)^(1/alpha))
   end
   return (b/k_n)*((dz_n)^(1/alpha))
end

function density_xi_n(x, alpha, b, theta, k_n, n)
  term1 = exp(-((b/k_n + x)^alpha))
  term2 = (b/k_n + x)^(alpha - 1)
  term3 = x^(theta + (n-1)*alpha)
  return term1*term2*term3
end

## test xi_n density

x_t = collect(range(0,1.2*10, length=100))
xi = Array{Float64}(undef, 10000)
for i in 1:10000
  xi[i] = sample_xi_n(0.75, 1, 0.75, 1,1)
end

## TEST for xi_n
using Plots
histogram(xi,bins=250)
x_t = collect(range(0,1.2*10^2, length=1000))
plot!(x_t,2000*map(x ->density_xi_n(x,0.75, 1.0, 0.75, 1, 1),x_t))

#histogram(sim.^(1/0.25),bins=250, normalize=true)
plot!(x_t,map(x ->density_xi_n(x,0.25, 1.0, 0.25, 1, 1),x_t)/20)


plot!(x_t,2000*map(x ->density_xi_n(x,0.75, 1.0, 0.75, 1, 1),x_t))

function SB_NGG(alpha_, b, theta,N)
  u_n =Array{Float64}(undef, N-1)
  k_n = 1
  for i in (1:(N-1))
    if i == 1
      k_n  = 1
    else
      k_n =prod(vec(ones(1,i-1)).-u_n[1:(i-1)])
    end
    #println(i)
    #f(x) = (x^((1/alpha_)-1))*exp(-((b/k_n + (x)^(1/alpha_))^alpha_))*((b/k_n + (x)^(1/alpha_))^(alpha_ - 1))* x^(theta/alpha_ + (i-1))
    #support = (0.0, Inf)
    # Build the sampler and simulate 10,000 samples
    #sampler = RejectionSampler(f, support)
    #sim = run_sampler!(sampler, 1)[1]
    #xi_n = (sim)^(1/alpha_)
    xi_n = sample_xi_n(alpha_, b, theta, k_n,i)
    println(xi_n)
    R"zn = rgamma(1,1- $alpha_, $b/$k_n + $xi_n)"
    #g = Gamma(1- alpha_, 1/(b/k_n + xi_n))
    #z_n  = rand(g,1)[1]
    z_n  = @rget zn
    #h_n  = b/k_n + xi_n
    R" xn <- retstable($alpha_, 1, h =$b/$k_n + $xi_n,method = 'MH')"
    x_n = @rget xn
    u_n[i] = z_n/ (z_n + x_n)
 end

  p = Array{Float64}(undef, N)
  p[1] = u_n[1]
  for l in 2:(N-1)
    p[l] = prod(vec(ones(1,l - 1)).-u_n[1:(l - 1)])*u_n[l]
  end
  p[N]=prod(vec(ones(1,(N-1))).- u_n)
  return p
end

n_tr=10
p_k = SB_NGG(.25, 1,0.25,n_tr)
plot(1:,p_k)



function Pkn_NGG_SB(n, β, σ, M; runs=10^4)
  theta = 1
  b = (σ*theta*β)^(1/σ)
  alpha_ = σ
    array_clust_num = Array{Int64}(undef, runs)
    for i in 1:runs
        weights_NGG_SB = SB_NGG(alpha_,b, theta, M)
        c = wsample(1:M, weights_NGG_SB, n, replace=true)
        n_c =  length(unique(c))
        array_clust_num[i] = n_c
    end
    return proportions(array_clust_num, n)
end



function Pkn_NGG_SB_1(n, β, σ, M; runs=10^4)
  theta = (σ*β)
  b = 1
  alpha_ = σ
    array_clust_num = Array{Int64}(undef, runs)
    for i in 1:runs
        weights_NGG_SB = SB_NGG(alpha_,b, theta, M)
        c = wsample(1:M, weights_NGG_SB, n, replace=true)
        n_c =  length(unique(c))
        array_clust_num[i] = n_c
    end
    return proportions(array_clust_num, n)
end

pks = Pkn_NGG_SB_1(20,1.0,0.25,50; runs=2*100)
pk_fk = Pkn_NGG_FK(20,1.0,0.75,50; runs=2*100)
p_true = Pkn_NGG_numeric.(1:10, 10, 1.0, 0.25)
#pks_1 = Pkn_NGG_SB_1(100,1.0,0.25,100; runs=2*100)

#E_1 =  pks|> ar -> map(*, ar, 1:100) |> sum
#E_2 =  pks_1|> ar -> map(*, ar, 1:100) |> sum
#println([E_1, E_2])

plot(collect(range(1,10, length=10)), pks)
#plot!(collect(range(1,10, length=10)), pk_fk)
plot!(collect(range(1,10, length=10)), p_true)
#plot!(collect(range(1,100, length=100)), pk_fk)

#E_1 =  pks_1|> ar -> map(*, ar, 1:100) |> sum
#E_2 =  pk_fk|> ar -> map(*, ar, 1:100) |> sum
#println([E_1, E_2])


#using Optim, OptimTestProblems


#function quantilepk(pk)
#  n = length(pk)
#  cdfk = cumsum(pk)
#  quantile_grid = (1:n)/(n+1)
#  return [argmin(cdfk .< q) for q in quantile_grid]
#end

#function gamma_qq(pk, α, β)
#  n = length(pk)
#  quantile_grid = (1:n)/(n+1)
#  return sum((quantile.(Gamma(α, 1/β), quantile_grid) - quantilepk(pk)).^2)
#end


 #function smooth_pk(Pkn)
#  lower = [1., 0.]
#  upper = [10000, 10000]
#  initial_x = [2.0, 2.0]
#  objfun_pkn(x) = gamma_qq(Pkn, x[1], x[2])
#  res = optimize(objfun_pkn, lower, upper, initial_x, SAMIN(), Optim.Options(iterations=10^6))
#  α_, β_ = res.minimizer
#  return Pkn__smooth = pdf.(Gamma(α_, 1/β_), eachindex(Pkn))
#end

#pks_smooth = smooth_pk(pks)
#Pkn_fk__smooth = smooth_pk(pk_fk)
#pks_smooth_1 =  smooth_pk(pks_1)


#plot(collect(range(1,100, length=100)), pks)
#plot(collect(range(1,100, length=100)), pks_smooth)
#plot!(collect(range(1,100, length=100)), Pkn_fk__smooth)
#plot!(collect(range(1,100, length=100)), pks_smooth_1)


#using AdaptiveRejectionSampling
#μ, σ, theta = 1.0, 0.75, 0.75
#kn = 1
#n = 1
#f(x) = (x^((1/σ)-1))*exp(-((μ/kn + (x)^(1/σ))^σ))*((μ/kn + (x)^(1/σ))^(σ - 1))* x^(theta/σ + (n-1))

#support = (0.0, Inf)

# Build the sampler and simulate 10,000 samples
#sampler = RejectionSampler(f, support)
#@time sim = run_sampler!(sampler, 100000);


pks = Pkn_NGG_SB_1(100,1.0,0.25,250; runs=2*100)

R"
pk_sb_1_25=$pks
save(pk_sb_1_25,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_100_1.Rdata')
"

pk_sb_1_75 = Pkn_NGG_SB_1(100,1.0,0.75,250; runs=2*100)
R"
pk_sb_1_75=$pk_sb_1_75
save(pk_sb_1_75,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_100_1.Rdata')
"

pk_sb_10_25 = Pkn_NGG_SB_1(100,10.0,0.25,250; runs=2*100)
R"
pk_sb_10_25=$pk_sb_10_25
save(pk_sb_10_25,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_100_1.Rdata')
"
pk_sb_10_75 = Pkn_NGG_SB_1(100,10.0,0.75,250; runs=2*100)
R"
pk_sb_10_75=$pk_sb_10_25
save(pk_sb_10_75,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_100_1.Rdata')
"

pk_sb_1_25_1000 = Pkn_NGG_SB_1(1000,1.0,0.25,250; runs=2*100)
R"
pk_sb_1_25_1000=$pk_sb_1_25_1000
save(pk_sb_1_25_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_1000_1.Rdata')
"
pk_sb_1_75_1000 = Pkn_NGG_SB_1(1000,1.0,0.75,250; runs=2*100)
R"
pk_sb_1_75_1000=$pk_sb_1_75_1000
save(pk_sb_1_75_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_1000_1.Rdata')
"
pk_sb_10_25_1000 = Pkn_NGG_SB_1(1000,10.0,0.25,250; runs=2*100)
R"
pk_sb_10_25_1000=$pk_sb_10_25_1000
save(pk_sb_10_25_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_1000_1.Rdata')
"
pk_sb_10_75_1000 = Pkn_NGG_SB_1(1000,10.0,0.75,250; runs=2*100)
R"
pk_sb_10_75_1000=$pk_sb_10_75_1000
save(pk_sb_10_75_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_1000_1.Rdata')
"


plot(collect(range(1,1000, length=1000)), pk_sb_10_25_1000)
plot!(collect(range(1,1000, length=1000)), pk_sb_1_75_1000)
E_1 =  pk_sb_10_25_1000|> ar -> map(*, ar, 1:1000) |> sum
E_2 =  pk_sb_10_25_1000|> ar -> map(*, ar, 1:1000) |> sum

pk_sb_10_25_1000 = Pkn_NGG_SB_1(1000,10.0,0.25,250; runs=2*100)



#####################################################################


pk_sb_1_25_100_1k = Pkn_NGG_SB(100,1.0,0.25,1000; runs=2*100)

R"
pk_sb_1_25_100_1k=$pk_sb_1_25_100_1k
save(pk_sb_1_25_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_100_1k.Rdata')
"

pk_sb_1_75_100_1k = Pkn_NGG_SB_1(100,1.0,0.75,1000; runs=2*100)

R"
pk_sb_1_75_100_1k=$pk_sb_1_75_100_1k
save(pk_sb_1_75_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_100_1k_.Rdata')
"

pk_sb_10_25_100_1k = Pkn_NGG_SB_1(100,10.0,0.25,1000; runs=2*100)

R"
pk_sb_10_25_100_1k=$pk_sb_10_25_100_1k
save(pk_sb_10_25_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_100_1k.Rdata')
"

pk_sb_10_75_100_1k = Pkn_NGG_SB_1(100,10.0,0.75,1000; runs=2*100)

R"
pk_sb_10_75_100_1k=$pk_sb_10_75_100_1k
save(pk_sb_10_75_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_100_1k.Rdata')
"





pk_sb_1_25_1k_1k = Pkn_NGG_SB_1(1000,1.0,0.25,1000; runs=2*100)

R"
pk_sb_1_25_1k_1k=$pk_sb_1_25_1k_1k
save(pk_sb_1_25_1k_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_1k_1k.Rdata')
"


pk_sb_1_75_1k_1k = Pkn_NGG_SB_1(1000,1.0,0.75,1000; runs=2*100)

R"
pk_sb_1_75_1k_1k=$pk_sb_1_75_1k_1k
save(pk_sb_1_75_1k_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_1k_1k.Rdata')
"

pk_sb_10_25_1k_1k = Pkn_NGG_SB_1(1000,10.0,0.25,1000; runs=2*100)

R"
pk_sb_10_25_1k_1k=$pk_sb_10_25_1k_1k
save(pk_sb_10_25_1k_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_1k_1k.Rdata')
"


pk_sb_10_75_1k_1k = Pkn_NGG_SB_1(1000,10.0,0.75,1000; runs=2*100)

R"
pk_sb_10_75_1k_1k=$pk_sb_10_75_1k_1k
save(pk_sb_10_75_1k_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_1k_1k.Rdata')
"
