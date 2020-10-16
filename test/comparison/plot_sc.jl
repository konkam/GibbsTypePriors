
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_FK_sampling.jl")

using DataFrames, DataFramesMeta, RCall
R"library(tidyverse)
library(latex2exp)
"

function approximation_prior_distribution(beta,sigma,N,Nt)
    #i = findall(x->x==sigma, sigma_arr)
    Pkn_numeric_ = Pkn_NGG_numeric.(1:N, N, beta, sigma)
    Pkn_order2_ = Pkn_NGG_pred_approx(N, beta, sigma)
    Pkn_NGGMult = Pkn_NGGM_precomp.(1:N,N,Nt,beta,sigma, [Pkn_numeric_])
    kappa = (beta*1*sigma)^(1/sigma)
    Pkn_NGG_FK_ = prior_Kn(1,kappa, sigma, N,3000; runs=10^4)[1].p_k
    df= DataFrame(Pkn_numeric = Pkn_numeric_,
                  Pkn_order2 = Pkn_order2_,
                  Pkn_NGGM = Pkn_NGGMult,
                  Pkn_NGG_FK= Pkn_NGG_FK_)
       return df
end

function plot_draw_prior_distribution(df,N,beta,sigma)
               R"p = ggplot(data.frame(k = 1:$N,
                                Pkn_numeric = $(df.Pkn_numeric[1:N]),
                                Pkn_order2 = $(df.Pkn_order2[1:N]),
                                Pkn_NGGM = $(df.Pkn_NGGM[1:N]),
                                Pkn_NGG_FK = $(df.Pkn_NGG_FK[1:N])
                            ) %>%
                    gather(Process_type, density, Pkn_numeric:Pkn_NGG_FK),
               aes(x=k, y = density, colour = Process_type)) +
                geom_line() + xlab('k') + ylab('') + ggtitle(TeX(sprintf('$\\tau =%2.f$, $\\sigma = %.2f$',$beta,$sigma)))+
               ggthemes::scale_colour_ptol() + theme_minimal()+
               theme(plot.title = element_text(hjust = 0.5))"
    return R"p"
end


#n=100
#R"setwd('~/Documents/GitHub/GibbsTypePriors/plots')
#load(file = 'FK_NGG_0.250.75beta_0.681420222312052Ns_100N_tr_3000.Rda')"
#@rget priors_list
n=100
beta= 1.0
ntr=100
sigma_vec= [0.25,0.75]
DF_all_100 = map(x ->approximation_prior_distribution(beta,x,n,ntr),sigma_vec)

P_all_approx_100 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
             P_all_approx_100[i]= plot_draw_prior_distribution(DF_all_100[i],n,beta,sigma_vec[i])

end

P_all_approx_100


R"
m1=as.list($P_all_approx_100)
save(m1,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/P_100_df_b1.Rdata')
"

n=100
beta= 10.0
ntr=100
sigma_vec= [0.25,0.75]
DF_all_100_10 = map(x ->approximation_prior_distribution(beta,x,n,ntr),sigma_vec)
N_plot = [30, 50]
P_all_approx_100_10 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
             P_all_approx_100_10[i]= plot_draw_prior_distribution(DF_all_100_10[i],n,beta,sigma_vec[i])

end

P_all_approx_100_10


R"
m2=as.list($P_all_approx_100_10)
save(m2,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/P_100_df_b10.Rdata')
"


R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_100)
m2 = as.list($P_all_approx_100_10)
df_100 = $DF_all_100
df_100_10 =  $DF_all_100_10
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow,ncol = 1,rel_heights = c(10, 1))
ggsave(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/Plots_sigma_all_approximation_100.pdf', width= 6, height = 4,p)
save(m1,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/P_100_df_b1.Rdata')
save(m2,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/P_100_df_b10.Rdata')
save(df_100,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/DF_100_df_b1.Rdata')
save(df_100_10,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/DF_100_df_b10.Rdata')
p"




#####################################################################################
#n=1000

#R"load(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/plots/FK_NGG_0.250.75beta_0.681420222312052Ns_1000N_tr_3000.Rda')"
#@rget priors_list

n=1000
beta= 1.0
ntr=500
sigma_vec= [0.25,0.75]
DF_all_1000 = map(x ->approximation_prior_distribution(beta,x,n,ntr),sigma_vec)
N_plot = [300, 500]
P_all_approx_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
             P_all_approx_1000[i]= plot_draw_prior_distribution(DF_all_1000[i],N_plot[i],beta,sigma_vec[i])
end

R"
m1=as.list($P_all_approx_1000)
df_1000 = $DF_all_1000
save(df_1000,file ='DF_1000_df_b1.Rdata')
save(m1,file ='P_1000_df_b1.Rdata')
"
 P_all_approx_1000

#R"load(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/plots/FK_NGG_0.250.75beta_14.6807536543832Ns_1000N_tr_3000.Rda')"
#@rget priors_list
n=1000
beta= 10.0
ntr=500
sigma_vec= [0.25,0.75]
DF_all_1000_10 = map(x ->approximation_prior_distribution(beta,x,n,ntr),sigma_vec)
N_plot = [300, 500]
P_all_approx_1000_10 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
             P_all_approx_1000_10[i]= plot_draw_prior_distribution(DF_all_1000_10[i],N_plot[i],beta,sigma_vec[i])

end

#P_all_approx_100
R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_1000)
m2 = as.list($P_all_approx_1000_10)
df_1000 = $DF_all_1000
df_1000_10 =  $DF_all_1000_10
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))
ggsave(file = 'Plots_approximation_1000.pdf', width= 6, height = 4,p)
#save(m1,file ='P_1000_df_b1.Rdata')
#save(m2,file ='P_1000_df_b10.Rdata')
#save(df_1000,file ='DF_1000_df_b1.Rdata')
#save(df_1000_10,file ='DF_1000_df_b10.Rdata')
p"

using BenchmarkTools
@btime Pkn_NGG_full_approximation(10, 1.0, 0.1, logxk)
@btime Pkn_NGG_numeric.(1:10, 10, 1.0, 0.1)
@btime Pkn_NGGM.(1:10,10,10,1.0,0.1)
a= Pkn_NGG_numeric.(1:10, 10, 1.0, 0.1)
@btime Pkn_NGGM_precomp.(1:10,10,10,1.0,0.1,[a])
