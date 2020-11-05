include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/bystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")

using DataFrames, DataFramesMeta, RCall

function EV_of_number_of_clusters_NGG(n,par_vec)
    pk_ngg = Pkn_NGG.(1:n, n, par_vec[1], par_vec[2])
    E_ngg =  pk_ngg|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg).^2
    V_ngg = pk_ngg |> ar -> map(*, ar, x) |> sum
    pk_nggm = Pkn_NGGM_precomp.(1:n,n,n,par_vec[1], par_vec[2], [pk_ngg])
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


sigma = collect(range(0.01,0.99, length=10))
#alpha = collect(range(0.05,20, length=5))
alpha = collect(exp.(range(log(0.05), log(200), length =10)))
grid = collect(Iterators.product(alpha, sigma))
#grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))
n=50
grid_vec = vec(grid)

EV = EV_of_number_of_clusters_NGG.(n,grid_vec)

EV_NGG= first.(EV)
EV_NGGM = last.(EV)
EV_NGG_approx = EV_of_number_of_clusters_NGG_approx.(n,grid_vec)

#using Plots
#plot(X, Y,seriestype = :scatter)

R"library(tidyverse)
library(latex2exp)
library(gridExtra)
library(cowplot)
library(grid)"


DF_NGG = DataFrame(Exp = first.(EV_NGG),Var =last.(EV_NGG), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGG.std = sqrt.(DF_NGG.Var)
DF_NGG


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_NGG
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
print(df_r[1:15,])
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_ngg <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.8, linetype='dotted') +
  xlim(0, 50)+ ylim(0,10)+ labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()
"
R"p_ngg"

R"
df_color=data.frame(sigma= $(last.(grid_vec)), beta_log=  log($(first.(grid_vec))))
df_color$beta_scaled= range01(df_color$beta_log)
print(head(df_color))"
R"df_color$color= rgb(df_color$beta_scaled,df_color$sigma,0.5,1)
df_color$color_sigma= rgb(0,df_color$sigma,0.5,1)
df_color$color_beta= rgb(df_color$beta_scaled,0,0.5,1)
p_color = ggplot(df_color,aes(x=beta_log, y = sigma, colour = color)) +
geom_point(alpha =1, size=3) +
labs(y='sigma', x='beta_scaled')+scale_colour_identity()+
theme_classic()"
R"p_color"
R"z = grid.arrange(p_ngg,p_color, nrow =1)
z
#ggsave(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG.pdf',z)
"


DF_NGGM = DataFrame(Exp = first.(EV_NGGM),Var =last.(EV_NGGM), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGGM.std = sqrt.(DF_NGGM.Var)
DF_NGGM.beta_log = log.(DF_NGGM.beta)
DF_NGGM


R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_NGGM
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
"
R"
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_nggm <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
         xlim(0, 50)+ ylim(0,10)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()
"


DF_NGGM_approx = DataFrame(Exp = first.(EV_NGG_approx),Var =last.(EV_NGG_approx), beta = first.(grid_vec), sigma = last.(grid_vec))
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
p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_nggm_approx <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
         xlim(0, 50)+ ylim(0,10)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()
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






##############################################################################
##############################################################################



function EV_of_number_of_clusters_PY(n,par_vec)
    #pk_py = Pkn_NGG.(1:n, n, par_vec[1], par_vec[2])
    E =  E_PY_exact(10, par_vec[1], par_vec[2])
    V  = V_PY_exact(10, par_vec[1], par_vec[2])
    return E,V
end

sigma = collect(range(0.01,0.99, length=200))
alpha = collect(range(0.05,200, length=200))

grid = collect(Iterators.product(alpha, sigma))
grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))
n= 100

grid_vec = vec(grid)
EV_PY = EV_of_number_of_clusters_PY.(n, grid_vec)


DF_PY= DataFrame(Exp = first.(EV_PY),Var =last.(EV_PY), alpha = first.(grid_vec), sigma = last.(grid_vec))
DF_PY.std = sqrt.(DF_PY.Var)
DF_PY

R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_PY
df_r$alpha_scaled= range01(df_r$alpha)
df_r$color = rgb(df_r$alpha_scaled, df_r$sigma,0.5,1)
p = ggplot(df_r,aes(x=Exp, y = std, colour = color)) +
geom_point(alpha =0.9, size=2) +
labs(y='Std', x='Expectation')+scale_colour_identity()+
theme_minimal()"
R"p"

R"
df_color=data.frame(sigma= $(last.(grid_vec)), beta=  $(first.(grid_vec)))
df_color$beta_scaled= range01(df_color$beta)
df_color$color= rgb(df_color$beta_scaled,df_color$sigma,0.5,1)
p_color = ggplot(df_color,aes(x=beta, y = sigma, colour = color)) +
geom_point(alpha =0.7, size=5) +
labs(y='sigma', x='beta_scaled')+scale_colour_identity()+
theme_minimal()"
R"p_color"
R"z = grid.arrange(p,p_color, nrow =1)
z
ggsave(file = 'Std_variance_PY.pdf',z)
"



function EV_of_number_of_clusters_PYM(n,H,par_vec)
    pk_pym = Pkn_PYM.(1:n, n, H,  par_vec[1], par_vec[2])
    E_pym =  pk_pym|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_pym).^2
    V_pym = pk_pym |> ar -> map(*, ar, x) |> sum
    return E_pym, V_pym
end



sigma = collect(range(0.01,0.99, length=5))
alpha = collect(range(0.05,200, length=5))

grid = collect(Iterators.product(alpha, sigma))
#grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))

grid_vec = vec(grid)

H= 10
EV_PYM = EV_of_number_of_clusters_PYM.(n,10,grid_vec)


##### Expected number of clusters  Dirichlet process
