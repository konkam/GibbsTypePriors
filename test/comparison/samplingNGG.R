library(tidyverse)
setwd("~/Documents/GitHub/GibbsTypePriors/plots")

#samplng NGG
M=100
alpha = 1
beta= 10 #kappa
sigma = 0.9 #gamma
Meps = 0.01

fill_v1 <- function(M, Mv, W, x) {
  v <- rep(NA, M)
  for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
  return(v)
}


fill_v2 <- function(M, Mv, W, N, x) {
  v <- rep(NA, M)
  iMv <- N
  for (i in seq(M)) {
    while (iMv > 0 && Mv[iMv] < W[i]) {
      iMv <- iMv - 1
      # print(paste(iMv, Mv[iMv], W[i]))
    }
    v[i] <- x[iMv + 1] # This index shift is to keep consistency  with previous version of the function, not necessary.
  }
  # for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
  return(v)
}


MvInv <-
  function(eps, u = 0, alpha = 1, beta = 1, gama = 1 / 2, N = 3001,M) # eps no longer required
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

TruncatedNGG <- function(alpha,beta, sigma, Nw, Meps = 0.01) {
  Js<- MvInv(eps = Meps, u = 0, alpha = alpha, beta = beta, gama = sigma, N = 50001,M=Nw) 
  Js_norm<- Js/sum(Js) 
  return(Js_norm)
}

Prior_on_K_trunc<- function(alpha, beta, sigma, N_s,N_tr, runs=10^4 ){
  array_nc<-c()
  i=1
  while (i<=runs){
    weights_NGG<- TruncatedNGG(alpha,beta, sigma, N=N_tr)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_NGG)
    n_c<- length(unique(c))
    array_nc[i]<- n_c
    i=i+1
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:N_s)[!1:N_s%in%p_k$k], 
                  p_k=rep(0,length((1:N_s)[!1:N_s%in%p_k$k])))
  pk_padded= rbind(p_k, p_zeros)
  E_k=sum((1:N_s)*(pk_padded$p_k))
  V_k = sum((((1:N_s) - E_k)^2) *(pk_padded$p_k))
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k ))
}

P<- Prior_on_K_trunc(alpha=1,beta=(10*0.75)^(1/0.75), sigma=0.75,N_s= 1000,N_tr=3000, runs=1000)


NGG_prior_K_dist<- function(sigma_array, Ntr, Ns, runs, beta_array, alpha=1){
  priors_list <- list()
  for (i in 1:length(sigma_array)){
  NGG_prior <- Prior_on_K_trunc(alpha=alpha,beta=beta_array[i], sigma=sigma_array[i],N_s= Ns,N_tr=Ntr,runs=runs)
  Pk_weights<- NGG_prior$pk
  priors_list[[i]]<- Pk_weights
  }
  save(priors_list, file=  paste0("FK_NGG_",sigma_array[1],sigma_array[length(sigma_array)],"beta_",beta_array[i],"Ns_",Ns,"N_tr_",Ntr,".Rda"))
}

NGG_prior_K_dist(sigma_array=c(0.25,0.75), Ntr=3000, Ns=100, runs=10^5, beta_array=(1*c(0.25,0.75))^(1/c(0.25,0.75)), alpha=1)
NGG_prior_K_dist(sigma_array=c(0.25,0.75), Ntr=3000, Ns=100, runs=10^5, beta_array=(10*c(0.25,0.75))^(1/c(0.25,0.75)), alpha=1)


NGG_prior_K_dist(sigma_array=c(0.25,0.75), Ntr=3000, Ns=1000, runs=10^5, beta_array=(1*c(0.25,0.75))^(1/c(0.25,0.75)), alpha=1)
NGG_prior_K_dist(sigma_array=c(0.25,0.75), Ntr=3000, Ns=1000, runs=10^5, beta_array=(10*c(0.25,0.75))^(1/c(0.25,0.75)), alpha=1)

################################################################################################