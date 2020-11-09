setwd("~/Documents/GitHub/GibbsTypePriors")

#samplng NGG
M=100
alpha = 1
gama= 10 #kappa
sigma = 0.9 #gamma
Meps = 0.01
z= rexp(M, 1)
source("MvInv.R")

MvInv <-
  function(eps, u = 0.5, alpha = 1, beta = 1, gama = 1 / 2, N = 3001,M) # eps no longer required
  {
    x <- -log(seq(from = exp(-1e-05), to = exp(-10), length = N))
    f <- alpha / gamma(1 - gama) * x^(-(1 + gama)) * exp(-(u +
                                                             beta) * x)
    dx <- diff(x)
    h <- (f[-1] + f[-N]) / 2
    Mv <- c(rev(cumsum(rev(dx[-N] * h[-N]))), 0)
    
    #M <- ceiling(thresholdGG(
    #  alpha = alpha,
    #  kappa = beta + u,
    #  gama = gama
    #)) # upper bound defined via the grid
    M <- max(10, M) # We wish to make sure we at least use a few jumps
    W <- rexp(n = M)
    W <- cumsum(W)
    # x_which_min = function(w){ # I guess that this function could be defined outside of MvInV
    #   x[which.min(Mv > w)]
    # }
    # x_which_min = Vectorize(x_which_min, vectorize.args = "w")
    # v <- x_which_min(W)
    if (M < 25) {
      ## This version is faster because it has no loop, but it involves many passes over Mv which has 3001 elements
      return(fill_v1(M, Mv, W, x))
    }
    else {
      return(fill_v2(M, Mv, W, N, x)) # Faster for large values of N
    }
  }

#Js<- MvInv(eps = Meps, u = 0, alpha = alpha, beta = gama, gama = sigma, N = 50001,M) 
#Js_norm<- Js/sum(Js) 
#plot(1:length(Js),Js_norm)


StickBreakingNGG <- function(alpha,gama, sigma, N, Meps = 0.01) {
  Js<- MvInv(eps = Meps, u = 0, alpha = alpha, beta = gama, gama = sigma, N = 50001,M=N) 
  Js_norm<- Js/sum(Js) 
  return(Js_norm)
}

Prior_on_K_SB<- function(alpha, gama, sigma, N_s,N_tr, runs=10^4 ){
  array_nc<-c()
  i=1
  while (i<=runs){
    weights_NGG<- StickBreakingNGG(alpha,gama, sigma, N=N_tr)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_NGG )
    n_c<- length(unique(c))
    array_nc[i]<- n_c
    i=i+1
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:N_tr)[!1:N_tr%in%p_k$k], 
                  p_k=rep(0,length((1:N_tr)[!1:N_tr%in%p_k$k])))
  pk_padded= rbind(p_k, p_zeros)
  E_k=sum((1:N_tr) *(pk_padded$p_k))
  V_k = sum((((1:N_tr) - E_k)^2) *(pk_padded$p_k))
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k ))
}

P<- Prior_on_K_SB(alpha=1,gama=1, sigma=0.8,N_s= 1000,N_tr=1000, runs=10)


NGG_weights<- function(sigma_array, Ntr, Ns, runs, gama, alpha=1){
  priors_list <- list()
  for (i in 1:length(sigma_array)){
  NGG_prior <- Prior_on_K_SB(alpha=alpha,gama=gama, sigma=sigma_array[i],N_s= Ns,N_tr=Ntr,runs=runs)
  Pk_weights<- NGG_prior$pk
  priors_list[[i]]<- Pk_weights
  save(Pk_weights, file=  paste0("FK_NGG_",sigma_array[i],"beta_",gama,"Ns_",Ns,"N_tr_",Ntr,".Rda"))
  }
  save(priors_list, file=  paste0("FK_NGG_",sigma_array[1],sigma_array[length(sigma_array)],"beta_",gama,"Ns_",Ns,"N_tr_",Ntr,".Rda"))
  
}
NGG_weights(sigma_array=c(0.25,0.75), Ntr=300, Ns=1000, runs=10^5, gama=10, alpha=1)

################################################################################################