include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/common_functions.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Cnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Vnk.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Expect_Kn.jl")
include("/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/src/Pkn.jl")

using DataFrames, DataFramesMeta, RCall

function EV_of_number_of_clusters_NGG(n,par_vec)
    pk_ngg = Pkn_NGG.(1:n, n, par_vec[1], par_vec[2])
    E_ngg =  pk_ngg|> ar -> map(*, ar, 1:n) |> sum
    x = ((1:n).-E_ngg).^2
    V_ngg = pk_ngg |> ar -> map(*, ar, x) |> sum
    pk_nggm = Pkn_NGGM_precomp.(1:n,n,n,par_vec[1], par_vec[2], [pk_ngg])
    E_nggm =  pk_nggm|> ar -> map(*, ar, 1:n) |> sum
    x_nggm = ((1:n).-E_nggm).^2
    V_nggm = pk_nggm |> ar -> map(*, ar, x) |> sum
    return [E_ngg, V_ngg] , [E_nggm, V_nggm]
end


sigma = collect(range(0.01,0.99, length=10))
alpha = collect(range(0.05,20, length=10))
grid = collect(Iterators.product(alpha, sigma))
#grid_borders = vcat(collect(Iterators.product(alpha[1], sigma)),collect(Iterators.product(alpha[nb], sigma)), collect(Iterators.product(alpha, sigma[1])),collect(Iterators.product(alpha, sigma[ns])))
n=50
grid_vec = vec(grid)

EV = EV_of_number_of_clusters_NGG.(n,grid_vec)

EV_NGG= first.(EV)
EV_NGGM = last.(EV)
#
#using Plots
#plot(X, Y,seriestype = :scatter)
using DataFrames, DataFramesMeta, RCall

R"library(tidyverse)
library(latex2exp)"


DF_NGG = DataFrame(Exp = first.(EV_NGG),Var =last.(EV_NGG), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGG.std = sqrt.(DF_NGG.Var)
DF_NGG



R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_NGG
df_r$beta_scaled= range01(df_r$beta)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
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
ggsave(file = 'Std_variance_NGG.pdf',z)
"


DF_NGGM = DataFrame(Exp = first.(EV_NGGM),Var =last.(EV_NGGM), beta = first.(grid_vec), sigma = last.(grid_vec))
DF_NGGM.std = sqrt.(DF_NGGM.Var)
DF_NGGM



R"
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r= $DF_NGGM
df_r$beta_scaled= range01(df_r$beta)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
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
ggsave(file = 'Std_variance_NGGM.pdf',z)
"


#############



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

using DataFrames, DataFramesMeta, RCall

R"library(tidyverse)
library(latex2exp)"

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
