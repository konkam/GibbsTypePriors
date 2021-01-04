
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
library(viridis)
"

using Optim, OptimTestProblems



function approximation_prior_distribution(beta,sigma,N,Nt,sigma_arr)
    i = findall(x->x==sigma, sigma_arr)
    Pkn_numeric_ = Pkn_NGG_numeric.(1:N, N, beta, sigma)
    Pkn_order2_ = Pkn_NGG_pred_approx(N, beta, sigma)
    Pkn_NGGMult = Pkn_NGGM_precomp.(1:N,N,Nt[i],beta,sigma, [Pkn_numeric_])
    #Pkn_FK = Pkn_NGG_FK(n, β, σ, M; runs=2*10^2)
    #prior_Kn(alpha,1.0, sigma, N,3000; runs=10^4)[1].p_k
    df= DataFrame(Pkn_numeric = Pkn_numeric_,
                  Pkn_order2 = Pkn_order2_,
                  Pkn_NGGM = Pkn_NGGMult)
       return df
end


function approximation_prior_distribution_precomp(beta,sigma,N,Nt,sigma_arr,df)
    i = findall(x->x==sigma, sigma_arr)
    Pkn_numeric_ = df[i[1]].Pkn_numeric
    Pkn_order2_ = df[i[1]].Pkn_order2
    Pkn_NGGMult = Pkn_NGGM_precomp.(1:N,N,Nt[i],beta,sigma, [Pkn_numeric_])
    df= DataFrame(Pkn_numeric = Pkn_numeric_,
                  Pkn_order2 = Pkn_order2_,
                  Pkn_NGGM = Pkn_NGGMult)
       return df
end

# scale_color_viridis(discrete=TRUE)
#ggthemes::scale_colour_ptol()
function plot_draw_prior_distribution(df,N,beta,sigma,y_l,x_lab,n,m)
               R"p = ggplot(data.frame(k = 1:$N,
                                Pkn_5_numeric = $(df.Pkn_numeric[1:N]),
                                Pkn_4_order2 = $(df.Pkn_order2[1:N]),
                                Pkn_3_NGG_sM = $(df.Pkn_NGGM[1:N]),
                                Pkn_2_NGG_FK = $(df.Pkn_FK[1:N]),
                                Pkn_1_NGG_SB = $(df.Pkn_SB[1:N])
                            ) %>%
                    gather(Process_type, density, Pkn_5_numeric:Pkn_1_NGG_SB),
               aes(x=k, y = density, colour = Process_type)) +
                geom_line(aes(linetype =Process_type) ) + xlab($x_lab) + scale_linetype_manual(values=c('solid','solid','solid','solid','dashed')) +
                ylab('') + ggtitle(TeX(sprintf('$\\tau =%2.f$, $\\sigma = %.2f$,$\\n = %3.f$',$beta,$sigma,$n)))+
              scale_color_viridis(discrete=TRUE, direction=-1)+ theme_minimal()+ylim(0,$y_l)+xlim(1,$N)+scale_x_continuous(limits = c(1, $N), expand = c(0, 0),breaks= c(1,seq(0,$N,length=5)[2:5]))+
               theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10), plot.margin = unit(c(1,$m, 0, 0),'pt'))"
    return R"p"
end


function quantilepk(pk)
  n = length(pk)
  cdfk = cumsum(pk)
  quantile_grid = (1:n)/(n+1)
  return [argmin(cdfk .< q) for q in quantile_grid]
end

function gamma_qq(pk, α, β)
  n = length(pk)
  quantile_grid = (1:n)/(n+1)
  return sum((quantile.(Gamma(α, 1/β), quantile_grid) - quantilepk(pk)).^2)
end


 function smooth_pk(Pkn)
  lower = [1., 0.]
  upper = [10000, 10000]
  initial_x = [2.0, 2.0]
  objfun_pkn(x) = gamma_qq(Pkn, x[1], x[2])
  res = optimize(objfun_pkn, lower, upper, initial_x, SAMIN(), Optim.Options(iterations=10^6))
  α_, β_ = res.minimizer
  return Pkn__smooth = pdf.(Gamma(α_, 1/β_), eachindex(Pkn))
end



R"
load(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/DF_100_df_b1_250.Rdata')
DF_all_100_1000 = df_100_250
"
@rget DF_all_100_1000

n=100
β= 1.0
ntr=[1000,1000]
sigma_vec= [0.25,0.75]
DF_all_100_1000 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_100_1000),sigma_vec)

R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_100_1k.Rdata')
pk_sb_1_25_100_1k = pk_sb_1_25_100_1k"
@rget pk_sb_1_25_100_1k

R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_100_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_100_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_100_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_100_1k.Rdata')
pk_sb_1_25_100_1k = pk_sb_1_25_100_1k
pk_sb_10_25_100_1k = pk_sb_10_25_100_1k
pk_sb_1_75_100_1k = pk_sb_1_75_100_1k
pk_sb_10_75_100_1k_ = pk_sb_10_75_100_1k"
@rget pk_sb_1_25_100_1k
@rget pk_sb_10_25_100_1k
@rget pk_sb_1_75_100_1k_
@rget pk_sb_1_75_100_1k


R"
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_100_1k_.Rdata')
pk_sb_1_75_100_1k_ = pk_sb_1_75_100_1k"
@rget pk_sb_1_75_100_1k_

Pkn_SB_1_25_smooth = smooth_pk(pk_sb_1_25_100_1k)
Pkn_SB_1_75_smooth = smooth_pk(pk_sb_1_75_100_1k_)


Pkn_NGG_FK_025_1_100 = Pkn_NGG_FK(n, 1.0, 0.25, 1000; runs=2*10^2)
Pkn_NGG_FK_075_1_100 = Pkn_NGG_FK(n, 1.0, 0.75, 1000; runs=2*10^2)
Pkn_FK_1_25_smooth = smooth_pk(Pkn_NGG_FK_025_1_100)
Pkn_FK_1_75_smooth = smooth_pk(Pkn_NGG_FK_075_1_100)

R"
Pkn_NGG_FK_025_1_100_1000 = $Pkn_NGG_FK_025_1_100
save(Pkn_NGG_FK_025_1_100_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_025_1_100_1k.Rdata')"
R"
Pkn_NGG_FK_075_1_100_1000 = $Pkn_NGG_FK_075_1_100
save(Pkn_NGG_FK_075_1_100_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_075_1_100_1k.Rdata')"


DF_all_100_1000[1].Pkn_FK = Pkn_FK_1_25_smooth
DF_all_100_1000[1].Pkn_SB = Pkn_SB_1_25_smooth

DF_all_100_1000[2].Pkn_FK = Pkn_FK_1_75_smooth
DF_all_100_1000[2].Pkn_SB = Pkn_SB_1_75_smooth

R"
df_100_1000 = $DF_all_100_1000
save(df_100_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_100_df_b1_1000.Rdata')"

R"load('/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_100_df_b1_1000.Rdata')"
DF_all_100_1000 =@rget df_100_1000

N_plot = [40, 100]
y_l = [0.3,0.1]
P_all_approx_100_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
x_lab=[" ","k"]
for i in (1:length(sigma_vec))
             P_all_approx_100_1000[i]= plot_draw_prior_distribution(DF_all_100_1000[i],N_plot[i],β,sigma_vec[i],y_l[i],x_lab[i],n,0)

end

R"
m1=as.list($P_all_approx_100_1000)
save(m1,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/P_100_df_b1_H1000.Rdata')



#R"
#load(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/DF_100_df_b10_250.Rdata')
#DF_all_100_10_1000 = df_100_10_250
#"
#@rget DF_all_100_10_1000

n=100
β= 10.0
ntr=[1000,1000]
sigma_vec= [0.25,0.75]
DF_all_100_10_1000 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_100_10_1000),sigma_vec)



Pkn_NGG_FK_025_10_100_1k = Pkn_NGG_FK(n, 10.0, 0.25, 1000; runs=2*10^2)
Pkn_NGG_FK_075_10_100_1k = Pkn_NGG_FK(n, 10.0, 0.75, 1000; runs=2*10^2)

R"
Pkn_NGG_FK_025_10_100_1k = $Pkn_NGG_FK_025_10_100_1k
save(Pkn_NGG_FK_025_10_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_025_10_100_1k.Rdata')"
R"
Pkn_NGG_FK_075_10_100_1k = $Pkn_NGG_FK_075_10_100_1k
save(Pkn_NGG_FK_075_10_100_1k,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_075_10_100_1k.Rdata')"


Pkn_FK_10_25_smooth_1k = smooth_pk(Pkn_NGG_FK_025_10_100_1k)
Pkn_FK_10_75_smooth_1k  = smooth_pk(Pkn_NGG_FK_075_10_100_1k)


Pkn_SB_10_25_smooth_1k  = smooth_pk(pk_sb_10_25_100_1k)
Pkn_SB_10_75_smooth_1k  = smooth_pk(pk_sb_10_75_100_1k)


DF_all_100_10_1000[1].Pkn_FK = Pkn_FK_10_25_smooth_1k
DF_all_100_10_1000[1].Pkn_SB = Pkn_SB_10_25_smooth_1k

DF_all_100_10_1000[2].Pkn_FK = Pkn_FK_10_75_smooth_1k
DF_all_100_10_1000[2].Pkn_SB = Pkn_SB_10_75_smooth_1k

R"
df_100_10_1000 = $DF_all_100_10_1000
save(df_100_10_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_100_df_b10_1000.Rdata')"

R"load('/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_100_df_b10_1000.Rdata')"
DF_all_100_10_1000 =@rget df_100_10_1000

N_plot = [40, 100]
y_l = [0.2,0.1]
P_all_approx_100_10_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
x_lab=[" ","k"]
for i in (1:length(sigma_vec))
             P_all_approx_100_10_1000[i]= plot_draw_prior_distribution(DF_all_100_10_1000[i],N_plot[i],β,sigma_vec[i],y_l[i],x_lab[i],n,10)

end



R"
m2=as.list($P_all_approx_100_10_1000)
save(m2,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/P_100_df_b10.Rdata')
"


R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_100_1000)
m2 = as.list($P_all_approx_100_10_1000)
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow,ncol = 1,rel_heights = c(10, 1))
ggsave(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Plots_approximation_100_H1000_v2.pdf', width= 6, height = 4,p)
p"




#####################################################################################
#n=1000

R"#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparisonp_ngg_1_25_n1000.Rdata')
#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparisonp_ngg_1_75_n1000.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b1.Rdata')
#p_ngg_1_25_pk = p_ngg_1_25$pk
#p_ngg_1_75_pk = p_ngg_1_75$pk"
#@rget p_ngg_1_25_pk
#@rget p_ngg_1_75_pk
DF_all_1000_1 = @rget df_1000

n=1000
β= 1.0
ntr=[1000,1000]
sigma_vec= [0.25,0.75]
DF_all_1000_1000 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_1000_1),sigma_vec)
#R"
#df_1000 = $DF_all_1000
#save(df_1000,file ='/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b1.Rdata')"
Pkn_NGG_FK_025_1_1000 = Pkn_NGG_FK(n, 1.0, 0.25, 1000; runs=2*10^2)
Pkn_NGG_FK_075_1_1000 = Pkn_NGG_FK(n, 1.0, 0.75, 1000; runs=2*10^2)
Pkn_NGG_FK_025_10_1000 = Pkn_NGG_FK(n, 10.0, 0.25, 1000; runs=2*10^2)
Pkn_NGG_FK_075_10_1000 = Pkn_NGG_FK(n, 10.0, 0.75, 1000; runs=2*10^2)


R"
Pkn_NGG_FK_025_1_1000 = $Pkn_NGG_FK_025_1_1000
save(Pkn_NGG_FK_025_1_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_025_1_1000_1k.Rdata')"
R"
Pkn_NGG_FK_075_1_1000 = $Pkn_NGG_FK_075_1_1000
save(Pkn_NGG_FK_075_1_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_075_1_1000_1k.Rdata')"
R"
Pkn_NGG_FK_025_10_1000 = $Pkn_NGG_FK_025_10_1000
save(Pkn_NGG_FK_025_10_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_025_10_1000_1k.Rdata')"
R"
Pkn_NGG_FK_075_10_1000 = $Pkn_NGG_FK_075_10_1000
save(Pkn_NGG_FK_075_10_1000,file ='/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Pkn_NGG_FK_075_10_1000_1k.Rdata')"

R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_1k_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_1k_1k.Rdata')
pk_sb_1_25_1k_1k = pk_sb_1_25_1k_1k
pk_sb_1_75_1k_1k = pk_sb_1_75_1k_1k"
@rget pk_sb_1_25_1k_1k
@rget pk_sb_1_75_1k_1k


Pkn_FK_1_25_1k_smooth = smooth_pk(Pkn_NGG_FK_025_1_1000)
Pkn_FK_1_75_1k_smooth = smooth_pk(Pkn_NGG_FK_075_1_1000)
Pkn_SB_1_25_1k_smooth = smooth_pk(pk_sb_1_25_1k_1k)
Pkn_SB_1_75_1k_smooth = smooth_pk(pk_sb_1_75_1k_1k)


DF_all_1000_1000[1].Pkn_FK = Pkn_FK_1_25_1k_smooth
DF_all_1000_1000[1].Pkn_SB = Pkn_SB_1_25_1k_smooth

DF_all_1000_1000[2].Pkn_FK = Pkn_FK_1_75_1k_smooth
DF_all_1000_1000[2].Pkn_SB = Pkn_SB_1_75_1k_smooth

R"
df_1000_1_1000 = $DF_all_1000_1000
save(df_1000_1_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_1000_df_b1_1000.Rdata')"

R"load('~/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_1000_df_b1_1000.Rdata')"
DF_all_1000_1000 = @rget df_1000_1_1000

N_plot = [100, 600]
y_l = [0.15,0.03]
x_lab=[" ","k"]

P_all_approx_1000_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
            P_all_approx_1000_1000[i]= plot_draw_prior_distribution(DF_all_1000_1000[i],N_plot[i],β,sigma_vec[i], y_l[i],x_lab[i],n,0)
end

# P_all_approx_100_10[i]= plot_draw_prior_distribution(DF_all_100_10[i],N_plot[i],β,sigma_vec[i],y_l[i],x_lab[i],n,10)

R"
m1=as.list($P_all_approx_1000_1000)
save(m1,file ='P_1000_df_b1_1000.Rdata')
"
#P_all_approx_1000


 R"
 p=$(P_all_approx_1000_1000[2])
 pdf(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Plot_all_1000_1_075_1000.pdf')
 plot(p)
 dev.off()"
#R"load(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/plots/FK_NGG_0.250.75beta_14.6807536543832Ns_1000N_tr_3000.Rda')"
#@rget priors_list

R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_1k_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_1k_1k.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b10.Rdata')
pk_sb_10_25_1k_1k = pk_sb_10_25_1k_1k
pk_sb_10_75_1k_1k = pk_sb_10_75_1k_1k
#DF_all_1000_10 = df_1000_10"
@rget pk_sb_10_25_1k_1k
@rget pk_sb_10_75_1k_1k
#@rget DF_all_1000_10
n=1000
β= 10.0
ntr=[1000,1000]
sigma_vec= [0.25,0.75]
#DF_all_1000_10 = map(x ->approximation_prior_distribution(β,x,n,ntr,sigma_vec),sigma_vec)
DF_all_1000_10_1000 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_1000_10),sigma_vec)

#E_1 =  Pkn_SB_10_25_1k__smooth|> ar -> map(*, ar, 1:1000) |> sum


Pkn_FK_10_25_1k_smooth = smooth_pk(Pkn_NGG_FK_025_10_1000)
Pkn_FK_10_75_1k_smooth = smooth_pk(Pkn_NGG_FK_075_10_1000)
Pkn_SB_10_25_1k__smooth = smooth_pk(pk_sb_10_25_1k_1k)
Pkn_SB_10_75_1k_smooth = smooth_pk(pk_sb_10_75_1k_1k)


DF_all_1000_10_1000[1].Pkn_FK = Pkn_FK_10_25_1k_smooth
DF_all_1000_10_1000[1].Pkn_SB = Pkn_SB_10_25_1k__smooth

DF_all_1000_10_1000[2].Pkn_SB = Pkn_SB_10_75_1k_smooth
DF_all_1000_10_1000[2].Pkn_FK = Pkn_FK_10_75_1k_smooth
R"
df_1000_10_1000 = $DF_all_1000_10_1000
save(df_1000_10_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_1000_df_b10_n1000.Rdata')"

R"load('~/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/DF_1000_df_b10_n1000.Rdata')"
DF_all_1000_10_1000 = @rget df_1000_10_1000

N_plot = [100, 600]
y_l = [0.15,0.03]
x_lab=[" ","k"]

P_all_approx_1000_10_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
            P_all_approx_1000_10_1000[i]= plot_draw_prior_distribution(DF_all_1000_10_1000[i],N_plot[i],β,sigma_vec[i], y_l[i],x_lab[i],n,10)
end



R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_1000_1000)
m2 = as.list($P_all_approx_1000_10_1000)
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))
ggsave(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/H1000/Plots_approximation_1000_1000_v2.pdf', width= 6, height = 4,p)
#save(m1,file ='P_1000_df_b1.Rdata')
#save(m2,file ='P_1000_df_b10.Rdata')
p"



##### TESTS

using BenchmarkTools


@btime Pkn_NGG_full_approximation(100, 1.0, 0.1, logxk)
@btime Pkn_NGG_numeric.(1:10, 10, 1.0, 0.1)
@btime Pkn_NGGM.(1:10,10,10,1.0,0.1)
a= Pkn_NGG_numeric.(1:10, 10, 1.0, 0.1)
@btime Pkn_NGGM_precomp.(1:10,10,10,1.0,0.1,[a])



using BenchmarkTools
n, β, σ, M = 100, 1., 0.05, 1000
@btime Pkn_FK = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)
n, β, σ, M = 100, 1., 0.05, 250
@btime Pkn_FK = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)

@btime Pkn_NGG_full_approximation(100, 1.0, 0.1, logxk)

@btime a = Pkn_NGG_numeric.(1:1000, 1000, 1.0, 0.25)
@btime Pkn_NGGM_precomp.(1:1000, 1000, 1.0, 0.25,[a])


n, β, σ, M = 100, 10., 0.25, 250
using BenchmarkTools
@btime p1 = Pkn_NGG_FK(n, β, σ, M; runs=10^2)
@btime p2 = Pkn_NGG_FK(n, β, σ, M; runs=10^2)

n, β, σ, M = 100, 1., 0.3, 250
@btime p2 = Pkn_NGG_FK(n, β, σ, M; runs=10^2)

####
@btime Pkn_NGG_FK_fast(n, β, σ, M; runs=1)

@btime Pkn_NGG_FK_fast(n, β, σ, M; runs=1)



n, β, σ, M = 100, 10, 0.25, 200


@time Pkn_FK = Pkn_NGG_FK(n, β, σ, M; runs=2*10^2)

@time Pkn_FK_fast = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)

@time Pkn_exact = Pkn_NGG(n, β, σ)


Pkn_FK_2 =  Pkn_NGG_FK(n, β, σ, M; runs=2*10^2)
Pkn_FK_3= Pkn_NGG_FK(n, β, σ, M; runs=2*10^2)

function gamma_l2(pk, α, β)
  n = length(pk)
  return sum((pdf.(Gamma(α, 1/β), 1:n) - pk).^2)
end

using Optim, OptimTestProblems

lower = [1., 0.]
upper = [10000, 10000]
initial_x = [2.0, 2.0]
objfun_FK(x) =gamma_l2(Pkn_FK, x[1], x[2])
res = optimize(objfun_FK, lower, upper, initial_x, SAMIN(), Optim.Options(iterations=10^6))
α_FK, β_FK = res.minimizer
Pkn_FK_smooth = pdf.(Gamma(α_FK, 1/β_FK), eachindex(Pkn_FK))


to_plot = DataFrame(k = 1:n, Pkn_FK = Pkn_FK, Pkn_FK_fast = Pkn_FK_fast,
Pkn_exact = Pkn_exact,
Pkn_FK_smooth = Pkn_FK_smooth)

p = R"$to_plot %>%
    as_tibble %>%
    gather(type, p_k, -k) %>%
    ggplot(aes(x = k, y = p_k, colour = type)) +
    theme_minimal() +
    geom_point() +
    geom_line()"

R"png('FK_smooth_L2_fit.png')
  plot($p)
  dev.off()"


  function quantilepk(pk)
    n = length(pk)
    cdfk = cumsum(pk)
    quantile_grid = (1:n)/(n+1)
    return [argmin(cdfk .< q) for q in quantile_grid]
  end

  function gamma_qq(pk, α, β)
    n = length(pk)
    quantile_grid = (1:n)/(n+1)
    return sum((quantile.(Gamma(α, 1/β), quantile_grid) - quantilepk(pk)).^2)
  end


  lower = [1., 0.]
  upper = [10000, 10000]
  initial_x = [2.0, 2.0]
  objfun_FK(x) = gamma_qq(Pkn_FK, x[1], x[2])
  res = optimize(objfun_FK, lower, upper, initial_x, SAMIN(), Optim.Options(iterations=10^6))
  α_FK, β_FK = res.minimizer
  Pkn_FK_smooth = pdf.(Gamma(α_FK, 1/β_FK), eachindex(Pkn_FK))

  lower = [1., 0.]
  upper = [10000, 10000]
  initial_x = [2.0, 2.0]
  objfun_SB(x) = gamma_qq(p_ngg_10_25_pk_100, x[1], x[2])
  res = optimize(objfun_SB, lower, upper, initial_x, SAMIN(), Optim.Options(iterations=10^6))
  α_SB, β_SB = res.minimizer
  Pkn_SB_smooth = pdf.(Gamma(α_SB, 1/β_SB), eachindex(p_ngg_10_25_pk_100))



  to_plot = DataFrame(k = 1:n, Pkn_FK = Pkn_FK, Pkn_FK_fast = Pkn_FK_fast, Pkn_exact = Pkn_exact, Pkn_FK_smooth = Pkn_FK_smooth,
  Pkn_SB = Pkn_SB_smooth)

  to_plot2 = DataFrame(k = 1:n, Pkn_exact = Pkn_exact, Pkn_FK_smooth = Pkn_FK_smooth, Pk_SB = p_ngg_10_25_pk_100, Pk_SB_sm= Pkn_SB_smooth, PK_FK2 = Pkn_FK_2)

  p = R"$to_plot %>%
      as_tibble %>%
      gather(type, p_k, -k) %>%
      ggplot(aes(x = k, y = p_k, colour = type)) +
      theme_minimal() +
      geom_point() +
      geom_line()"

  R"png('FK_smooth_qq_fit.png')
    plot($p)
    dev.off()"


plot(collect(range(1,100, length=100)),Pkn_FK_smooth)
plot!(collect(range(1,100, length=100)),Pkn_SB_smooth)
plot!(collect(range(1,100, length=100)),p_ngg_10_25_pk_100)


n=112

using Distributions
Alpha = rand(Gamma(1.939519, 1/0.3114096),1000)
mean(Alpha)
DF_DPM = map(x ->Pkn_Dirichlet_Mult.(1:112,112,112,x),Alpha)
DF_DPM_16= hcat(DF_DPM...)'
Ar =mean(DF_DPM_16,dims=1)
L= Array{Float64}(undef, 112)
for l in (1:112)
  L[l]= Ar[1,l]
end
L |> ar -> map(*, ar, 1:112) |> sum

to_plot = DataFrame(k = 1:n,
Pkn_DPM_1 = vcat(Pkn_Dirichlet_Mult.(1:16, 112,16, 112.0),zeros(96)),
Pkn_DPM_2 = vcat(Pkn_Dirichlet_Mult.(1:112, 112,112, 6.23)),
Pkn_DPM_3 = L)

#Pkn_DPM_3 = Pkn_Dirichlet.(1:112, 112, 6.23))

p = R"$to_plot %>%
    as_tibble %>%
    gather(type, p_k, -k) %>%
    ggplot(aes(x = k, y = p_k, colour = type)) + scale_color_manual(values=c('#B63679FF', '#FB8861FF','#31688EFF'), labels = c(expression(DPM[1]),expression(DPM[2]), expression(DP[c])), name ='Model') +

    theme_minimal() +
    geom_point() +
    ylab('Probability')+
    theme_bw()  +
    theme(axis.text.x = element_text(angle = 0, hjust = 1,size = 10), strip.text = element_text(size = 10),legend.position = 'right', plot.title = element_text(hjust = 0.5),legend.text.align = 0)+
    geom_line()"

R"pdf('DP_mult_3.pdf', width= 6, height = 4)
  plot($p)
  dev.off()"
