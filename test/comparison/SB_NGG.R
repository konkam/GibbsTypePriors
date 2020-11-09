

## sampling NGG stick breaking
library(tidyverse)
library(copula)
R_q <- function(q, b, k_n,alpha){
  term1 = exp(-((b/k_n)^alpha)*((1 + q^(1/alpha))^alpha - q))
  #term2 = ((q^(1/alpha))/(1+q^(1/alpha)))^(1-alpha)
  term2 = (1+q^(1/alpha))^(alpha-1)
  #term3 = (b/kn)^(theta + alpha*(n-1) -1)
  R_q = term1*term2
  return(R_q)
}


sample_xi_n<- function(alpha, b, theta, k_n,n){
  repeat{
        r_n <- runif(1,0,1)
        xi_n <- rgamma(1,theta/alpha + n,(b/k_n)^(alpha))
        if (r_n <= R_q(xi_n, b, k_n, alpha)){ 
          return( (b/k_n)*((xi_n)^(1/alpha)))
          }
        }
}

#density

density_phi_n<- function(x, alpha, b, theta, k_n,n){
  term1 = exp(-((b/k_n)^alpha)*(1 + x^(1/alpha))^alpha - x)
  term2 = (1 + x^(1/alpha))^(alpha - 1)
  return (term1*term2)
}
density_xi_n<- function(x, alpha, b, theta, k_n,n){
  term1 = exp(-((b/k_n + x)^alpha))
  term2 = (b/k_n + x)^(alpha - 1)
  term3 = x^(theta + (n-1)*alpha)
  return (term1*term2*term3)
}

xi<- vector("numeric", length=10000)
for(i in 1:10000){
  xi[i]<- sample_xi_n(0.75, 1, 1, 0.5,1)
}
hist(xi)


##Test for the sampling of xi_n
x= seq(0,100,length=5000)
hist(xi, breaks=25,probability = TRUE)
#x= seq(0,100,length=1000)
#plot(x,density_phi_n(x,0.5, 1, 1, 0.5,1))
lines(x,density_xi_n(x,0.75,1,1,0.5,1)*2.5)


SB_NGG<- function(alpha, b, theta,N){
u_n <- vector("numeric",length=N-1)
for(i in 1:(N-1)){
  if (i==1) {
  k_n <- 1 
  } else {
    k_n = prod(1- u_n[1:i-1])
  }
  xi_n <- sample_xi_n(alpha, b, theta, k_n,i)
  z_n <- rgamma(1,1- alpha, b/k_n + xi_n)
  x_n <- retstable(alpha, 1, h = b/k_n + xi_n)
  u_n[i] = z_n/(z_n + x_n)
} 
p <- vector("numeric",length=N)
p[1] <- u_n[1]
for(l in 2:(N-1))p[l] <- prod(1 - u_n[1:(l - 1)])*u_n[l]
p[N] <- prod(1 - u_n)   
p
}

n_tr=1000
p_k <- SB_NGG(.75, 1,7.5,n_tr)
plot(1:n_tr,p_k )

Prior_on_K_SB_NGG<- function(alpha,  b, theta, N_s,N_tr, runs=10^4 ){
  array_nc<-c()
  i=1
  err=0
  while (i<=runs){
    weights_NGG<- SB_NGG(alpha, b, theta,N_tr)
    #plot(1:N_tr, weights_NGG)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_NGG )
    n_c<- length(unique(c))
    #n_c
    array_nc[i]<- n_c
    i=i+1
    err = err +weights_NGG[N_tr]
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:N_s)[!1:N_s%in%p_k$k], 
                  p_k=rep(0,length((1:N_s)[!1:N_s%in%p_k$k])))
  pk_padded= rbind(p_k, p_zeros) %>% arrange(k)
  E_k=sum((1:N_s) *(pk_padded$p_k))
  V_k = sum((((1:N_s) - E_k)^2) *(pk_padded$p_k))
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=sqrt(V_k), error = err/runs ))
}

p_ngg_10_25 <- Prior_on_K_SB_NGG(0.25, 1,2.5, 100,300, runs=10^3)
plot(1:100, p_ngg_10_25$pk, type="l")

p_ngg <- Prior_on_K_SB_NGG(0.25, 39.0625,1, 100,300, runs=10^3)
plot(1:100, p_ngg$pk, type="l")


p_ngg_1_025 <- Prior_on_K_SB_NGG(0.25, 1,0.25, 100,100, runs=10^3)
plot(1:100, p_ngg_1_025$pk, type="l")


p_ngg_1_075 <- Prior_on_K_SB_NGG(0.75, 1,0.75, 100,500, runs=10^3)
p_ngg
plot(1:100, p_ngg_1_075$pk, type="l")


p_ngg_10_075 <- Prior_on_K_SB_NGG(0.75, 1, 7.5, 100,2000, runs=10^2)
p_ngg
plot(1:100, p_ngg_10_075$pk, type="l")





p_ngg <- Prior_on_K_SB_NGG(0.5, 1, 0.5, 100,250, runs=10^3)
plot(1:100, p_ngg$pk, type="l")


p_ngg <- Prior_on_K_SB_NGG(0.7, 0.1, 0.1, 100,100, runs=10)
plot(1:100, p_ngg$pk, type="l")



