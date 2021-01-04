
include("~/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("~//Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("~//Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("~//Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("~//Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("~/Documents/GitHub/GibbsTypePriors/src/MvInv.jl")
include("~//Documents/GitHub/GibbsTypePriors/test/comparison/NGG_FK_sampling.jl")

using DataFrames, DataFramesMeta, RCall
R"library(tidyverse)
library(latex2exp)
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


function plot_draw_prior_distribution(df,N,beta,sigma,y_l,x_lab,n,m)
               R"p = ggplot(data.frame(k = 1:$N,
                                Pkn_numeric = $(df.Pkn_numeric[1:N]),
                                Pkn_order2 = $(df.Pkn_order2[1:N]),
                                Pkn_NGG_sM = $(df.Pkn_NGGM[1:N]),
                                Pkn_NGG_FK = $(df.Pkn_FK[1:N]),
                                Pkn_NGG_SB = $(df.Pkn_SB[1:N])
                            ) %>%
                    gather(Process_type, density, Pkn_numeric:Pkn_NGG_SB),
               aes(x=k, y = density, colour = Process_type)) +
                geom_line(aes(linetype =Process_type) ) + xlab($x_lab) + scale_linetype_manual(values=c('solid','solid','solid','solid','dashed')) +
                ylab('') + ggtitle(TeX(sprintf('$\\tau =%2.f$, $\\sigma = %.2f$,$\\n = %3.f$',$beta,$sigma,$n)))+
               ggthemes::scale_colour_ptol() + theme_minimal()+ylim(0,$y_l)+xlim(1,$N)+scale_x_continuous(limits = c(1, $N), expand = c(0, 0),breaks= c(1,seq(0,$N,length=5)[2:5]))+
               theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10), plot.margin = unit(c(1,$m, 0, 0),'pt'))"
    return R"p"
end


function plot_draw_prior_distribution_v2(df,N,beta,sigma)
               R"p = ggplot(data.frame(k = 1:$N,
                                Pkn_numeric = $(df.Pkn_numeric[1:N]),
                                Pkn_order2 = $(df.Pkn_order2[1:N]),
                                Pkn_NGGM = $(df.Pkn_NGGM[1:N])
                            ) %>%
                    gather(Process_type, density, Pkn_numeric:Pkn_NGGM),
               aes(x=k, y = density, colour = Process_type)) +geom_line()
               "
               R"p<- p + geom_line(data = data.frame(k = 1:$N,
                                Pkn_FK = $(df.Pkn_FK[1:N]),
                                Pkn_SB = $(df.Pkn_SB[1:N])) %>% gather(Process_type, density, Pkn_FK:Pkn_SB),
                           aes(x=k, y = density, alpha = 0.)) + geom_smooth(se=FALSE, span =0.2)+
                xlab('k') + ylab('') + x_lim(1,N)+ ggtitle(TeX(sprintf('$\\tau =%2.f$, $\\sigma = %.2f$',$beta,$sigma)))+
               ggthemes::scale_colour_ptol() + theme_minimal()+
               theme(plot.title = element_text(hjust = 0.5))"
    return R"p"
end



function plot_draw_prior_distribution_v3(df,N,beta,sigma, smooth, y_max)
               R"p = ggplot(data.frame(k = 1:$N,
                   Pkn_FK = $(df.Pkn_FK[1:N]),
                   Pkn_SB = $(df.Pkn_SB[1:N])
                            ) %>%
                    gather(Process_type, density, Pkn_FK:Pkn_SB),
               aes(x=k, y = density, colour = Process_type)) + geom_line(se=FALSE, span =$smooth, size=0.5)
               "
               R"p<- p + geom_line(data = data.frame(k = 1:$N,
                   Pkn_numeric = $(df.Pkn_numeric[1:N]),
                   Pkn_order2 = $(df.Pkn_order2[1:N]),
                   Pkn_NGGM = $(df.Pkn_NGGM[1:N])) %>% gather(Process_type, density, Pkn_numeric:Pkn_NGGM),
                           aes(x=k, y = density, colour = Process_type)) +
                xlab('k') + ylab('') + xlim(c(1,$N))+ ggtitle(TeX(sprintf('$\\tau =%2.f$, $\\sigma = %.2f$',$beta,$sigma)))+
               ggthemes::scale_colour_ptol() + theme_minimal()+ylim(0,$y_max) +
               theme(plot.title = element_text(hjust = 0.5))"
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
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_100_df_b1_250.Rdata')
DF_all_100_10 = df_100_10_250
"
@rget DF_all_100_10

n=100
β= 1.0
ntr=[250,250]
sigma_vec= [0.25,0.75]
DF_all_100 = map(x ->approximation_prior_distribution(β,x,n,ntr,sigma_vec),sigma_vec)


R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_100_1.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_100_1.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_100_1.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_100_1.Rdata')
pk_sb_1_25_100 = pk_sb_1_25
pk_sb_1_75_100 = pk_sb_1_75
pk_sb_10_25_100 = pk_sb_10_25
pk_sb_10_75_100 = pk_sb_10_75"
@rget pk_sb_1_25_100
@rget pk_sb_1_75_100
@rget pk_sb_10_25_100
@rget pk_sb_10_75_100


Pkn_NGG_FK_025_1_100 = Pkn_NGG_FK(n, 1.0, 0.25, 250; runs=2*10^2)
Pkn_NGG_FK_075_1_100 = Pkn_NGG_FK(n, 1.0, 0.75, 250; runs=2*10^2)

Pkn_FK_1_25_smooth = smooth_pk(Pkn_NGG_FK_025_1_100)
Pkn_FK_1_75_smooth = smooth_pk(Pkn_NGG_FK_075_1_100)
Pkn_SB_1_25_smooth = smooth_pk(pk_sb_1_25_100)
Pkn_SB_1_75_smooth = smooth_pk(pk_sb_1_75_100)

DF_all_100[1].Pkn_FK = Pkn_FK_1_25_smooth
DF_all_100[1].Pkn_SB = Pkn_SB_1_25_smooth

DF_all_100[2].Pkn_FK = Pkn_FK_1_75_smooth
DF_all_100[2].Pkn_SB = Pkn_SB_1_75_smooth

R"
df_100_250 = $DF_all_100
save(df_100_250,file ='~//Documents/GitHub/GibbsTypePriors/test/comparison/DF_100_df_b1_250.Rdata')"


N_plot = [40, 100]
y_l = [0.2,0.1]
P_all_approx_100 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
x_lab=[" ","k"]
for i in (1:length(sigma_vec))
             P_all_approx_100[i]= plot_draw_prior_distribution(DF_all_100[i],N_plot[i],β,sigma_vec[i],y_l[i],x_lab[i],n,0)

end


R"
m1=as.list($P_all_approx_100)
save(m1,file ='~/Documents/GitHub/GibbsTypePriors/P_100_df_b1.Rdata')


n=100
β= 10.0
ntr=[250,250]
sigma_vec= [0.25,0.75]
DF_all_100_10 = map(x ->approximation_prior_distribution(β,x,n,ntr,sigma_vec),sigma_vec)

Pkn_NGG_FK_025_10_100 = Pkn_NGG_FK(n, 10.0, 0.25, 250; runs=2*10^2)
Pkn_NGG_FK_075_10_100 = Pkn_NGG_FK(n, 10.0, 0.75, 250; runs=2*10^2)


Pkn_FK_10_25_smooth = smooth_pk(Pkn_NGG_FK_025_10_100)
Pkn_FK_10_75_smooth = smooth_pk(Pkn_NGG_FK_075_10_100)
Pkn_SB_10_25_smooth = smooth_pk(pk_sb_10_25_100)
Pkn_SB_10_75_smooth = smooth_pk(pk_sb_10_75_100)


DF_all_100_10[1].Pkn_FK = Pkn_FK_10_25_smooth
DF_all_100_10[1].Pkn_SB = Pkn_SB_10_25_smooth

DF_all_100_10[2].Pkn_FK = Pkn_FK_10_75_smooth
DF_all_100_10[2].Pkn_SB = Pkn_SB_10_75_smooth

R"
df_100_10_250 = $DF_all_100_10
save(df_100_10_250,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_100_df_b10_250.Rdata')"



N_plot = [40, 100]
y_l = [0.2,0.1]
P_all_approx_100_10 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
x_lab=[" ","k"]
for i in (1:length(sigma_vec))
             P_all_approx_100_10[i]= plot_draw_prior_distribution(DF_all_100_10[i],N_plot[i],β,sigma_vec[i],y_l[i],x_lab[i],n,10)

end

P_all_approx_100_10


R"
m2=as.list($P_all_approx_100_10)
save(m2,file ='~/Documents/GitHub/GibbsTypePriors/P_100_df_b10.Rdata')
"


R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_100)
m2 = as.list($P_all_approx_100_10)
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow,ncol = 1,rel_heights = c(10, 1))
ggsave(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/Plots_sigma_all_approximation_100_n.pdf', width= 6, height = 4,p)
#save(m1,file ='~/Documents/GitHub/GibbsTypePriors/P_100_df_b1.Rdata')
#save(m2,file ='~/Documents/GitHub/GibbsTypePriors/P_100_df_b10.Rdata')
p"



#####################################################################################

R"#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparisonp_ngg_1_25_n1000.Rdata')
#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparisonp_ngg_1_75_n1000.Rdata')
#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b1.Rdata')
#p_ngg_1_25_pk = p_ngg_1_25$pk
#p_ngg_1_75_pk = p_ngg_1_75$pk"

n=1000
β= 1.0
ntr=[250,250]
sigma_vec= [0.25,0.75]
DF_all_1000 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_1000),sigma_vec)
#R"
#df_1000 = $DF_all_1000
#save(df_1000,file ='/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b1.Rdata')"
Pkn_NGG_FK_025_1_1000 = Pkn_NGG_FK(n, 1.0, 0.25, 250; runs=2*10^2)
Pkn_NGG_FK_075_1_1000 = Pkn_NGG_FK(n, 1.0, 0.75, 250; runs=2*10^2)
Pkn_NGG_FK_025_10_1000 = Pkn_NGG_FK(n, 10.0, 0.25, 250; runs=2*10^2)
Pkn_NGG_FK_075_10_1000 = Pkn_NGG_FK(n, 10.0, 0.75, 250; runs=2*10^2)



R"
Pkn_NGG_FK_025_1_1000 = $Pkn_NGG_FK_025_1_1000
save(Pkn_NGG_FK_025_1_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/Pkn_NGG_FK_025_1_1000.Rdata')"
R"
Pkn_NGG_FK_075_1_1000 = $Pkn_NGG_FK_075_1_1000
save(Pkn_NGG_FK_075_1_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/Pkn_NGG_FK_075_1_1000.Rdata')"
R"
Pkn_NGG_FK_025_10_1000 = $Pkn_NGG_FK_025_10_1000
save(Pkn_NGG_FK_025_10_1000,file ='~//Documents/GitHub/GibbsTypePriors/test/comparison/Pkn_NGG_FK_025_10_1000.Rdata')"
R"
Pkn_NGG_FK_075_10_1000 = $Pkn_NGG_FK_075_10_1000
save(Pkn_NGG_FK_075_10_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/Pkn_NGG_FK_075_10_1000.Rdata')"

R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_25_1000_1.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_1_75_1000_1.Rdata')
pk_sb_1_25_1000 = pk_sb_1_25_1000
pk_sb_1_75_1000 = pk_sb_1_75_1000"
@rget pk_sb_1_25_1000
@rget pk_sb_1_75_1000


Pkn_FK_1_25_1k_smooth = smooth_pk(Pkn_NGG_FK_025_1_1000)
Pkn_FK_1_75_1k_smooth = smooth_pk(Pkn_NGG_FK_075_1_1000)
Pkn_SB_1_25_1k_smooth = smooth_pk(pk_sb_1_25_1000)
Pkn_SB_1_75_1k_smooth = smooth_pk(pk_sb_1_75_1000)


DF_all_1000[1].Pkn_FK = Pkn_FK_1_25_1k_smooth
DF_all_1000[1].Pkn_SB = Pkn_SB_1_25_1k_smooth

DF_all_1000[2].Pkn_FK = Pkn_FK_1_75_1k_smooth
DF_all_1000[2].Pkn_SB = Pkn_SB_1_75_1k_smooth

R"
df_1000 = $DF_all_1000
save(df_1000,file ='~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b1.Rdata')"


N_plot = [100, 600]
y_l = [0.15,0.05]
x_lab=[" ","k"]

P_all_approx_1000 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
            P_all_approx_1000[i]= plot_draw_prior_distribution(DF_all_1000[i],N_plot[i],β,sigma_vec[i], y_l[i],x_lab[i],n,0)
end

R"
m1=as.list($P_all_approx_1000)
df_1000 = $DF_all_1000
save(df_1000,file ='DF_1000_df_b1.Rdata')
save(m1,file ='P_1000_df_b1.Rdata')
"
#
 R"
 p=$(P_all_approx_1000[2])
 pdf(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/Plot_all_1000_1_075.pdf')
 plot(p)
 dev.off()"


R"load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_25_1000_1.Rdata')
load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/NGG_sb_10_75_1000_1.Rdata')
#load(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/DF_1000_df_b10.Rdata')
pk_sb_10_25_1000 = pk_sb_10_25_1000
pk_sb_10_75_1000 = pk_sb_10_75_1000
#DF_all_1000_10 = df_1000_10"
@rget pk_sb_10_25_1000
@rget pk_sb_10_75_1000
@rget DF_all_1000_10
n=1000
β= 10.0
ntr=[250,250]
ntr_fk=500
c= [0.25,0.75]
DF_all_1000_10 = map(x ->approximation_prior_distribution_precomp(β,x,n,ntr,sigma_vec,DF_all_1000_10),sigma_vec)


Pkn_FK_10_25_1k_smooth = smooth_pk(Pkn_NGG_FK_025_10_1000)
Pkn_FK_10_75_1k_smooth = smooth_pk(Pkn_NGG_FK_075_10_1000)
Pkn_SB_10_25_1k__smooth = smooth_pk(pk_sb_10_25_1000)
Pkn_SB_10_75_1k_smooth = smooth_pk(pk_sb_10_75_1000)


DF_all_1000_10[1].Pkn_FK = Pkn_FK_10_25_1k_smooth
DF_all_1000_10[1].Pkn_SB = Pkn_SB_10_25_1k__smooth

DF_all_1000_10[2].Pkn_SB = Pkn_SB_10_75_1k_smooth
DF_all_1000_10[2].Pkn_FK = Pkn_FK_10_75_1k_smooth
R"
df_1000_10 = $DF_all_1000_10
save(df_1000_10,file ='DF_1000_df_b10_n1000.Rdata')"


N_plot = [100, 600]
y_l = [0.15,0.05]
x_lab=[" ","k"]

P_all_approx_1000_10 =Array{RObject{VecSxp}}(undef,length(sigma_vec))
for i in (1:length(sigma_vec))
            P_all_approx_1000_10[i]= plot_draw_prior_distribution(DF_all_1000_10[i],N_plot[i],β,sigma_vec[i], y_l[i],x_lab[i],n,10)
end


P_all_approx_1000_10
R"library(gridExtra)
library(cowplot)
m1=as.list($P_all_approx_1000)
m2 = as.list($P_all_approx_1000_10)
prow <- plot_grid(
  m1[[1]] + theme(legend.position='none'),
  m2[[1]] + theme(legend.position='none'),
  m1[[2]] + theme(legend.position='none'),
  m2[[2]]+ theme(legend.position='none'),
  nrow = 2
)
legend_b <- get_legend(m1[[1]]+theme(legend.position ='top'))
p <- plot_grid(prow, ncol = 1,rel_heights = c(10, 1))
ggsave(file = '~/Documents/GitHub/GibbsTypePriors/test/comparison/Plots_approximation_1000.pdf', width= 6, height = 4,p)
save(m1,file ='P_1000_df_b1.Rdata')
save(m2,file ='P_1000_df_b10.Rdata')
p"
