include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/MvInv.jl")

using DataFrames, DataFramesMeta, RCall

function EV_of_number_of_clusters_NGG(n,H,par_vec)
    pk_ngg = Pkn_NGG.(1:n, n, par_vec[1], par_vec[2])
    E_ngg =  pk_ngg|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg).^2
    V_ngg = pk_ngg |> ar -> map(*, ar, x) |> sum
    pk_nggm = Pkn_NGGM_precomp.(1:n,n,H,par_vec[1], par_vec[2], [pk_ngg])
    E_nggm =  pk_nggm|> ar -> map(*, ar, 1:n) |> sum
    x_nggm = ((1:n).-E_nggm).^2
    V_nggm = pk_nggm |> ar -> map(*, ar, x_nggm) |> sum
    return [E_ngg, V_ngg] , [E_nggm, V_nggm]
end


function EV_of_number_of_clusters_NGG_approx(n,par_vec)
    pk_ngg_approx = Pkn_NGG_pred_approx(n,par_vec[1], par_vec[2])
    E_ngg_approx =  pk_ngg_approx|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg_approx).^2
    V_ngg_approx = pk_ngg_approx |> ar -> map(*, ar, x) |> sum
    return E_ngg_approx, V_ngg_approx
end


function EV_of_number_of_clusters_NGG_FK(n,par_vec,M)
    pk_ngg_fk = Pkn_NGG_FK_fast(n,par_vec[1], par_vec[2],M)
    E_ngg_fk =  pk_ngg_fk|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg_fk).^2
    V_ngg_fk = pk_ngg_fk |> ar -> map(*, ar, x) |> sum
    return E_ngg_fk, V_ngg_fk
end



function EV_of_number_of_clusters_NGG_FK_slow(n,par_vec,M)
    println([par_vec[2],par_vec[1]])
    if ((par_vec[2] < 0.35)&&(par_vec[2] < 10))
        return -1,-1
    end
    pk_ngg_fk = Pkn_NGG_FK(n,par_vec[1], par_vec[2], 250; runs = 2*10^2)
    E_ngg_fk =  pk_ngg_fk|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg_fk).^2
    V_ngg_fk = pk_ngg_fk |> ar -> map(*, ar, x) |> sum
    return E_ngg_fk, V_ngg_fk
end


sigma = collect(range(0.05,0.99, length=10))
#alpha = collect(range(0.05,20, length=5))
alpha = collect(exp.(range(log(1), log(200), length =10)))
grid = collect(Iterators.product(alpha, sigma))
#grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))
n=100
grid_vec = vec(grid)
H = 250


EV = EV_of_number_of_clusters_NGG.(n,H,grid_vec)
EV_NGG= first.(EV)
EV_NGGM = last.(EV)

EV_NGG_approx = EV_of_number_of_clusters_NGG_approx.(n,grid_vec)

#R"
#EV_NGG = $EV_NGG
#save(EV_NGG,file ='EV_NGG.Rdata')"
#R"
R"
load(file ='/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/EV_NGG_approx.Rdata')
"
@rget EV_NGG_approx

R"
load(file ='/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/EV_NGG.Rdata')
"
@rget EV_NGG


R"
load(file ='/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/EV_NGGM.Rdata')
"
@rget EV_NGGM
#R"
#EV_NGG_approx = $(DataFrame(EV_NGG_approx))
#save(EV_NGG_approx,file ='EV_NGG_approx.Rdata')"

#R"
#EV_NGGM = $EV_NGGM
#save(EV_NGGM,file ='EV_NGGM.Rdata')"

#R"
#EV_NGG_FK = $(DataFrame(EV_NGG_FK))
#save(EV_NGG_FK,file ='EV_NGG_FK.Rdata')"


EV_NGG_FK = EV_of_number_of_clusters_NGG_FK.(n,grid_vec,250)


#EV_NGG_FK_2 = EV_of_number_of_clusters_NGG_FK_slow.(n,grid_vec,250)

R"library(tidyverse)
library(latex2exp)
library(gridExtra)
library(cowplot)
library(grid)"





DF_NGG = DataFrame(Exp = first.(EV_NGG),Var =last.(EV_NGG), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGG.std = sqrt.(DF_NGG.Var)
DF_NGG.beta_log = log.(DF_NGG.beta)
DF_NGG


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_NGG
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_ngg <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.8, linetype='dotted') +
  xlim(0, 100)+ ylim(0,20)+ labs(y='Std', x='Expectation')+scale_colour_identity()+
  #ggtitle('NGG')+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10), plot.margin = unit(c(0,0, 0, 0),'pt'))
"
R"p_ngg
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG.pdf', width= 4, height = 4)
plot(p_ngg)
dev.off()"


R"
df_color=data.frame(sigma= $(last.(grid_vec)), beta_log=  log($(first.(grid_vec))))
df_color$beta_scaled= range01(df_color$beta_log)
df_color$color= rgb(df_color$beta_scaled,df_color$sigma,0.5,1)
df_color$color_sigma= rgb(0,df_color$sigma,0.5,1)
df_color$color_beta= rgb(df_color$beta_scaled,0,0.5,1)
p <- df_color %>% ggplot(aes(x=beta_log, y = sigma, group=sigma)) + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash')
p_color <- p  + geom_line(data =df_color, aes(x=beta_log, y = sigma, group=beta_log, color = color_beta), alpha = 0.8, linetype='dotted') +
  labs(y='Std', x='Expectation')+
  #ggtitle('NGG')+
  theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
p_color<- p_color + geom_point(data =df_color, aes(x=beta_log, y = sigma, colour=color), size=4)
"



R"p_color = ggplot(df_color,aes(x=beta_log, y = sigma, colour = color)) +
geom_point(alpha =1, size=4) + labs(y='sigma', x='beta_scaled')+scale_colour_identity()+ theme_classic()
p_color<- p_color + geom_line(data =df_color, aes(x=beta_log, y = sigma, group=beta_log), alpha = 0.5,size=0.3)
p_color<- p_color + geom_line(data =df_color, aes(x=beta_log, y = sigma, group=sigma),, alpha = 0.5,size=0.3)
p_color<- p_color + theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none', plot.margin = unit(c(5,0, 0, 0),'pt'))+
ylab(TeX(sprintf('$\\sigma$')))+ xlab(TeX(sprintf('$\\beta$')))+ scale_x_continuous(breaks=c(0,2,4),labels=c(expression(e^0),expression(e^2),expression(e^4)))
"
R"p_color
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGGs_color_grid.pdf',width= 4, height = 4)
plot(p_color)
dev.off()
"
R"z = grid.arrange(p_ngg,p_color, nrow =1)
z
#ggsave(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG.pdf',z)
"


DF_NGGM = DataFrame(Exp = first.(EV_NGGM),Var =last.(EV_NGGM), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGGM.std = sqrt.(DF_NGGM.Var)
DF_NGGM.beta_log = log.(DF_NGGM.beta)
DF_NGGM


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}"
R"
df_r= $DF_NGGM
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
"
R"
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_nggm <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
         xlim(0, 100)+ ylim(0,20)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
 pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGGM.pdf',width= 4, height = 4)
 plot(p_nggm)
 dev.off()
"


DF_NGGM_approx = DataFrame(Exp = EV_NGG_approx[1],Var =EV_NGG_approx[2], beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGGM_approx.std = sqrt.(DF_NGGM_approx.Var)
DF_NGGM_approx.beta_log= log.(DF_NGGM_approx.beta)


R"
df_r= $DF_NGGM_approx
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
"
R"
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_nggm_approx <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
         xlim(0, 100)+ ylim(0,20)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
 pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_approx.pdf',width= 4, height = 4)
 plot(p_nggm_approx)
 dev.off()

"


R"z = grid.arrange(p_ngg,p_nggm, p_color, nrow =1)
z
z =plot_grid(p_ngg,p_nggm, p_color, align = 'h', nrow = 1, rel_heights = c(4/10, 4/10, 1/5), rel_width = c(3/8, 3/8, 1/4))
z
#ggsave(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_NGGM.pdf',z)
"
R"
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_NGGM.pdf')
pushViewport(viewport(layout = grid.layout(2,2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p_ngg, vp = vplayout(2, 1))
print(p_nggm_approx, vp = vplayout(1, 1))
print(p_nggm, vp = vplayout(1, 2))
print(p_color, vp = vplayout(2, 2))
#ggsave(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_NGGM.pdf')
dev.off()
"


DF_NGG_FK = DataFrame(Exp = first.(EV_NGG_FK),Var =last.(EV_NGG_FK), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGG_FK.std = sqrt.(DF_NGG_FK.Var)
DF_NGG_FK.beta_log= log.(DF_NGG_FK.beta)


R"
df_r= $DF_NGG_FK
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
"
R"
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_ngg_fk <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
         xlim(0, 100)+ ylim(0,20)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+scale_colour_identity()+theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
 pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_FK.pdf',width= 4, height = 4)
 plot(p_ngg_fk)
 dev.off()
"

















using BenchmarkTools
n, β, σ, M = 100, 1., 0.05, 1000
@btime Pkn_FK = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)
n, β, σ, M = 100, 1., 0.05, 250
@btime Pkn_FK = Pkn_NGG_FK_fast(n, β, σ, M; runs=2*10^2)

@btime a = Pkn_NGG_numeric.(1:1000, 1000, 1.0, 0.25)
@btime Pkn_NGGM_precomp.(1:1000, 1000, 1.0, 0.25,[a])



R"
library(ggalt)
px = ggplot(data.frame(k = 1:1000,
     Pkn_numeric = $Pkn_FK),
aes(x=k, y =  Pkn_numeric, colour = 'red')) +
geom_point()+  geom_smooth(se=FALSE, span=0.2)
px
"

##############################################################################
##############################################################################



function EV_of_number_of_clusters_PY(n,par_vec)
    #pk_py = Pkn_NGG.(1:n, n, par_vec[1], par_vec[2])
    E =  E_PY_exact(n, par_vec[1], par_vec[2])
    V  = V_PY_exact(n, par_vec[1], par_vec[2])
    return E,V
end


sigma = collect(range(0.05,0.99, length=10))
#alpha = collect(range(0.05,20, length=5))
alpha = collect(exp.(range(log(1), log(200), length =10)))
grid = collect(Iterators.product(alpha, sigma))
#grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))
n=100
grid_vec = vec(grid)
EV_PY = EV_of_number_of_clusters_PY.(n, grid_vec)


DF_PY= DataFrame(Exp = first.(EV_PY),Var =last.(EV_PY), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_PY.std = sqrt.(DF_PY.Var)
DF_PY.beta_log = log.(DF_PY.beta)

DF_PY


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_PY
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
print(df_r[1:15,])
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_py <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.8, linetype='dotted') +
  xlim(0, 100)+ ylim(0,20)+ labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')

"

R"p_py
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_PY_100.pdf')
plot(p_py)
dev.off()"





function EV_of_number_of_clusters_PYM(n,H,par_vec)
    pk_pym = Pkn_PYM.(1:n, n, H,  par_vec[1], par_vec[2])
    E_pym =  pk_pym|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_pym).^2
    V_pym = pk_pym |> ar -> map(*, ar, x) |> sum
    return E_pym, V_pym
end


H= 250
EV_PYM = EV_of_number_of_clusters_PYM.(n,H,grid_vec)


DF_PYM= DataFrame(Exp = first.(EV_PYM),Var =last.(EV_PYM), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_PYM.std = sqrt.(DF_PYM.Var)
DF_PYM.beta_log = log.(DF_PYM.beta)

DF_PYM


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_PYM
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
print(df_r[1:15,])
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_pym <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.8, linetype='dotted') +
  xlim(0, 100)+ ylim(0,20)+ labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()
"

R"p_pym
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_PYM_250_100.pdf')
plot(p_pym)
dev.off()"


R"
df_color=data.frame(sigma= $(last.(grid_vec)), beta_log=  log($(first.(grid_vec))))
df_color$beta_scaled= range01(df_color$beta_log)
df_color$color= rgb(df_color$beta_scaled,df_color$sigma,0.5,1)
df_color$color_sigma= rgb(0,df_color$sigma,0.5,1)
df_color$color_beta= rgb(df_color$beta_scaled,0,0.5,1)
p <- df_color %>% ggplot(aes(x=beta_log, y = sigma, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_color<- p +  geom_line(data =df_color, aes(x=beta_log, y = sigma, group=beta_log), alpha = 0.8, linetype='dotted') +theme_classic()+theme(legend.position='none')
p_color
"


R"p_color
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_PYM_PY_color.pdf')
plot(p_color)
dev.off()"
##### Expected number of clusters  Dirichlet process
