

## sampling NGG stick breaking
library(tidyverse)
library(copula)

R_q <- function(q, b, k_n,alpha){
  term1 = exp(-((b/k_n)^alpha)*((1 + q^(1/alpha))^alpha - q))
  #term2 = ((q^(1/alpha))/(1+q^(1/alpha)))^(1-alpha)
  term2 = (1+q^(1/alpha))^(alpha-1)
  #term3 = (b/k_n)^(theta + alpha*(n-1) -1)
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

xi<- vector("numeric", length=35000)
for(i in 1:35000){
  xi[i]<- sample_xi_n(0.75, 1, 1, 1,1)
}
hist(xi)


x= seq(0,100,length=200)
plot(density(xi))
#x= seq(0,100,length=1000)
#plot(x,density_phi_n(x,0.5, 1, 1, 0.5,1))
hist(xi, breaks=40,probability = TRUE)
lines(x,density_xi_n(x,0.75,1,1,1,1)*1.4)


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

n_tr=100
p_k <- SB_NGG(.75, 1,7.5,n_tr)
plot(1:n_tr,p_k )

n_tr=100
p_k <- SB_NGG(0.25, 1,2.5,n_tr)
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

#p_ngg_10_25 <- Prior_on_K_SB_NGG(0.25, 1,2.5, 100,100, runs=100)
#plot(1:100, p_ngg_10_25$pk, type="l")



########PY-SB

StickBreakingPY <- function(alpha,sigma, N) {
  a<- rep(1-sigma,N-1)
  b<- alpha+ sigma*(1:N-1)
  V <- rbeta(N-1 , a, b)
  p    <- vector("numeric",length=N)
  p[1] <- V[1]
  for(l in 2:(N - 1))p[l] <- prod(1 - V[1:(l - 1)])*V[l]
  p[N] <- prod(1 - V)   
  p
  return(p)
}

Prior_on_K_SB_PY<- function(alpha, sigma, N_s=100,N_tr=100, runs=10^4 ){
  array_nc<-c()
  i=1
  while (i<=runs){
    weights_PY<- StickBreakingPY(alpha,sigma,N_tr)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_PY )
    n_c<- length(unique(c))
    array_nc[i]<- n_c
    i=i+1
  }
  p_k = tibble(k=as.numeric(names(table(array_nc))), 
               p_k=as.vector(table(array_nc))/sum(table(array_nc)))
  p_zeros= tibble(k=(1:N_s)[!1:N_s%in%p_k$k], 
                  p_k=rep(0,length((1:N_s)[!1:N_s%in%p_k$k])))
  pk_padded= rbind(p_k, p_zeros)%>% arrange(k)
  E_k=sum((1:N_s) *(pk_padded$p_k))
  V_k = sum((((1:N_s) - E_k)^2) *(pk_padded$p_k))
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k ))
}

folderpath = "/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/comparison/"


pk_py_1_25 <- Prior_on_K_SB_PY(1.0, 0.25,100, 250, runs=2*10^2)
save(pk_py_1_25, file = paste0(folderpath,"pk_py_1_25.Rdata"))


pk_py_10_25 <- Prior_on_K_SB_PY(10.0, 0.25,100, 250, runs=2*10^2)
save(pk_py_10_25, file = paste0(folderpath,"pk_py_10_25.Rdata"))


pk_py_1_75 <- Prior_on_K_SB_PY(1.0, 0.75,100, 250, runs=2*10^2)
save(pk_py_1_75, file = paste0(folderpath,"pk_py_1_75.Rdata"))


pk_py_10_75 <- Prior_on_K_SB_PY(10.0, 0.75,100, 250, runs=2*10^2)
save(pk_py_10_75, file = paste0(folderpath,"pk_py_10_75.Rdata"))




pk_py_1_25 <- Prior_on_K_SB_PY(1.0, 0.25,1000, 250, runs=2*10^2)
save(pk_py_1_25, file = paste0(folderpath,"pk_py_1_25_1000.Rdata"))


pk_py_10_25 <- Prior_on_K_SB_PY(10.0, 0.25,1000, 250, runs=2*10^2)
save(pk_py_10_25, file = paste0(folderpath,"pk_py_10_25_1000.Rdata"))


pk_py_10_75 <- Prior_on_K_SB_PY(10.0, 0.75,1000, 250, runs=2*10^2)
save(pk_py_10_75, file = paste0(folderpath,"pk_py_10_75_1000.Rdata"))


pk_py_1_75 <- Prior_on_K_SB_PY(1.0, 0.75,1000, 250, runs=2*10^2)
save(pk_py_1_75, file = paste0(folderpath,"pk_py_1_75_1000.Rdata"))





grid<- expand.grid(beta = exp(seq(from = log(1),to=log(200), length=10)),
                   sigma = seq(from = 0.05, to=0.99, length= 10))


exp_var_SB_PY<- function(par_vec,n=100, M=250){
  pk = Prior_on_K_SB_PY(par_vec[1],par_vec[2],n,M, runs=100)
  return(list(E=pk$Ek, V= pk$Vk ))
}


df_100 = data.frame(matrix(ncol=2, nrow=100))
for (i in 1:100){
  par_vec = grid[i,]
  pk = Prior_on_K_SB_PY(par_vec[[1]],par_vec[[2]],100,250, runs=1000)
  df_100[i,1]=  pk$Ek
  df_100[i,2]=  pk$Vk
}
names(df_100)<-c("E","V")
df_100$Std = sqrt(df_100$V)
df_100$beta = unlist(grid[1])
df_100$sigma = unlist(grid[2])
df_100$beta_log = log(unlist(grid[1]))


range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r = df_100
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=E, y = Std, group=sigma))  + geom_line(alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_py_sb <- p  + geom_line(data =df_r, aes(x=E, y = Std, group=beta), alpha = 0.8, linetype='dotted') +
xlim(0, 100)+ ylim(0,15)+ labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic() + theme(plot.title = element_text(hjust = 0.5,size = 15), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), plot.margin = unit(c(0,0, 0, 0),'pt'))

pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_PY_SB_250_100.pdf',width= 4, height = 4)
plot(p_py_sb)
dev.off()


p <- df_r %>% ggplot(aes(x=Exp, y = std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_nggm_approx <- p  + geom_line(data =df_r, aes(x=Exp, y = std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
  xlim(0, 100)+ ylim(0,20)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_approx.pdf',width= 4, height = 4)
plot(p_nggm_approx)
dev.off()


