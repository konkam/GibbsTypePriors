

## sampling NGG stick breaking
library(tidyverse)
library(copula)
R_q <- function(q, b, k_n,alpha){
  term1 = exp(-((b/k_n)^alpha)*((1 + q^(1/alpha))^alpha - q))
  term2 = ((q^(1/alpha))/(1+q^(1/alpha)))^(1-alpha)
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

density_phi_n<- function(x, alpha, b, theta, k_n,n){
  term1 = exp(-((b/k_n + x^(1/alpha))^alpha))
  term2 = (b/k_n + x^(1/alpha))^(alpha - 1)
  term3 = x^(theta/alpha + (n-1))
  return (term1*term2*(x^((1/alpha)-1)))
}

log_phi_n<- function(x, alpha, b, theta, k_n,ns){
  term1 = -((b/k_n + x^(1/alpha))^alpha)
  term2 = (alpha - 1)*log(b/k_n + x^(1/alpha))
  term3 = (theta/alpha + (ns-1) + 1/alpha - 1)*log(x)
  return (term1 + term2)
}

density_xi_n<- function(x, alpha, b, theta, k_n,n){
  term1 = exp(-((b/k_n + x)^alpha))
  term2 = (b/k_n + x)^(alpha - 1)
  term3 = x^(theta + (n-1)*alpha)
  return (term1*term2*term3)
}


xi<- vector("numeric", length=100000)
for(i in 1:100000){
  xi[i]<- sample_xi_n(0.25, 0.00390625,1, 1,1)
}
hist(xi)


x= seq(0,100,length=200)
plot(density(xi))
#x= seq(0,100,length=1000)
#plot(x,density_phi_n(x,0.5, 1, 1, 0.5,1))
hist(xi, breaks=100, probability = TRUE)
lines(x,density_xi_n(x,0.75,1,0.75,1,1)*1.6)


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
  h = b/k_n + xi_n
  x_n <- retstable(alpha, 1, h = b/k_n + xi_n)
  u_n[i] = z_n/(z_n + x_n)
} 
p <- vector("numeric",length=N)
p[1] <- u_n[1]
for(l in 2:(N-1)){
  p[l] <- prod(1 - u_n[1:(l - 1)])*u_n[l]
}
p[N] <- prod(1 - u_n)   
p
}

n_tr=10
p_k <- SB_NGG(.25, 1,2.5,n_tr)
plot(1:n_tr,p_k )

Prior_on_K_SB_NGG<- function(beta, sigma, N_s,N_tr, runs=10^4 ){
  theta =1
  alpha = sigma
  b = (sigma*theta*beta)^(1/sigma)
  array_nc<-c()
  i=1
  err=0
  while (i<=runs){
    weights_NGG<- SB_NGG(alpha, b, theta,N_tr)
    #plot(1:N_tr, weights_NGG)
    c <- sample(1:N_tr,size=N_s, replace=TRUE, prob=weights_NGG )
    n_c<- length(unique(c))
    n_c
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
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k, error = err/runs ))
}


Prior_on_K_SB_NGG_1<- function(beta, sigma, N_s,N_tr, runs=10^4 ){
  #theta =1
  theta =sigma*beta
  alpha = sigma
  #b = (sigma*theta*beta)^(1/sigma)
  b=1
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
  return(list(pk=pk_padded$p_k, Ek= E_k,Vk=V_k, error = err/runs ))
}

folderpath =  "~/Documents/GitHub/GibbsTypePriors/test/comparison/"




p_ngg_1_25 <- Prior_on_K_SB_NGG(1.0, 0.25, 100,250, runs=2*10^2)
p_ngg_1_25_1 <- Prior_on_K_SB_NGG_1(1.0, 0.25, 100,250, runs=2*10^2)

plot(1:100, p_ngg_1_25_1$pk, type="l")
lines(1:100, p_ngg_1_25$pk, type="l")
#save(p_ngg_1_25_1, file = paste0(folderpath,"p_ngg_1_25.Rdata"))


#p_ngg_10_25 <- Prior_on_K_SB_NGG(10.0, 0.25, 100,250, runs=2*10^2)
p_ngg_10_25_1 <- Prior_on_K_SB_NGG_1(10.0, 0.25, 100,250, runs=2*10^2)

plot(1:100, p_ngg_10_25$pk, type="l")
lines(1:100, p_ngg_10_25_1$pk, type="l")
#save(p_ngg_10_25_1, file = paste0(folderpath,"p_ngg_10_25.Rdata"))

#p_ngg_1_5 <- Prior_on_K_SB_NGG(1.0, 0.5,100, 250, runs=2*10^3)
#p_ngg_1_5_1 <- Prior_on_K_SB_NGG_1(1.0, 0.5,100, 250, runs=2*10^3)
#plot(1:100, p_ngg_1_5$pk, type="l")
#lines(1:100, p_ngg_1_5_1$pk, type="l")


#p_ngg_1_1 <- Prior_on_K_SB_NGG(1.0, 0.1,100, 250, runs=2*10^2)
#p_ngg_1_1_1 <- Prior_on_K_SB_NGG_1(1.0, 0.1,100, 250, runs=2*10^2)
#plot(1:100, p_ngg_1_1$pk, type="l")
#lines(1:100, p_ngg_1_1_1$pk, type="l")


#p_ngg_10_1 <- Prior_on_K_SB_NGG(10.0, 0.7,100, 250, runs=2*10^2)
#p_ngg_10_1_1 <- Prior_on_K_SB_NGG_1(10.0, 0.7,100, 250, runs=2*10^2)
#plot(1:100, p_ngg_10_1$pk, type="l")
#lines(1:100, p_ngg_10_1_1$pk, type="l")


#p_ngg_1_75 <- Prior_on_K_SB_NGG(1.0, 0.75,100, 250, runs=2*10^2)
p_ngg_1_75_1 <- Prior_on_K_SB_NGG_1(1.0, 0.75,100, 250, runs=2*10^2)
plot(1:100, p_ngg_1_75_1$pk, type="l")
lines(1:100, p_ngg_1_75$pk, type="l")
#save(p_ngg_1_75_1, file = paste0(folderpath,"p_ngg_1_75_25.Rdata"))

p_ngg_10_75 <- Prior_on_K_SB_NGG(10.0, 0.75,100, 250, runs=2*10^2)
p_ngg_10_75_1 <- Prior_on_K_SB_NGG_1(10.0, 0.75,100, 250, runs=2*10^2)
plot(1:100, p_ngg_10_75_1$pk, type="l")
#save(p_ngg_10_75_1, file = paste0(folderpath,"p_ngg_10_75_25.Rdata"))


#p_ngg_1_25_1000 <- Prior_on_K_SB_NGG(1.0, 0.25, 100,1000, runs=2*10^2)
p_ngg_1_25_1000_1 <- Prior_on_K_SB_NGG_1(1.0, 0.25, 100,1000, runs=2*10^2)

plot(1:100, p_ngg_1_25_1000_1$pk, type="l")
folderpath =  "~/Documents/GitHub/GibbsTypePriors/test/comparison/"
#save(p_ngg_1_25_1000_1, file = paste0(folderpath,"p_ngg_1_25_1000.Rdata"))


p_ngg_10_25_1000_1 <- Prior_on_K_SB_NGG_1(10.0, 0.25, 100,1000, runs=2*10^2)
plot(1:100, p_ngg_10_25_1000_1$pk, type="l")
#save(p_ngg_10_25_1000_1, file = paste0(folderpath,"p_ngg_10_25_1000.Rdata"))


p_ngg_1_75_1000_1 <- Prior_on_K_SB_NGG_1(1.0, 0.75,100, 1000, runs=2*10^2)
plot(1:100, p_ngg_1_75_1000_1$pk, type="l")
#save(p_ngg_1_75_1000_1, file = paste0(folderpath,"p_ngg_1_75_1000.Rdata"))



p_ngg_10_75_1000_1 <- Prior_on_K_SB_NGG_1(10.0, 0.75, 100,1000, runs=2*10^2)
plot(1:100, p_ngg_10_75$pk, type="l")
#save(p_ngg_10_75_1000_1, file = paste0(folderpath,"p_ngg_10_75_1000_1.Rdata"))







####### n=1000
folderpath =  "~/Documents/GitHub/GibbsTypePriors/test/comparison/"
Nsamples = 1000

p_ngg_1_25_1000_1k <- Prior_on_K_SB_NGG_1(1.0, 0.25, 1000,1000, runs=2*10^2)
plot(1:1000, p_ngg_1_25_1000_1k$pk, type="l")
save(p_ngg_1_25_1000_1k, file = paste0(folderpath,"p_ngg_1_25_1k_1k.Rdata"))


p_ngg_10_25_1000_1k <- Prior_on_K_SB_NGG_1(10.0, 0.25, 1000,1000, runs=2*10^2)
plot(1:1000, p_ngg_10_25_1000_1k$pk, type="l")
save(p_ngg_10_25_1000_1k, file = paste0(folderpath,"p_ngg_10_25_1k_1k.Rdata"))

p_ngg_1_75_1000_1k <- Prior_on_K_SB_NGG_1(1.0, 0.75,1000, 1000, runs=2*10^2)
plot(1:1000, p_ngg_1_75_1000_1k$pk, type="l")
save(p_ngg_1_75_1000_1k, file = paste0(folderpath,"p_ngg_1_75_1k_1k.Rdata"))

p_ngg_10_75_1000_1k <- Prior_on_K_SB_NGG_1(10.0, 0.75, 1000,1000, runs=2*10^2)
plot(1:1000, p_ngg_10_75_1000_1k$pk, type="l")
save(p_ngg_10_75_1000_1k, file = paste0(folderpath,"p_ngg_10_75_1k_1k.Rdata"))









########################################################


grid<- expand.grid(beta = exp(seq(from = log(1),to=log(200), length=10)),
                      sigma = seq(from = 0.05, to=0.99, length= 10))


exp_var_SB_NGG<- function(par_vec,n=100, M=100){
 pk = Prior_on_K_SB_NGG_1(par_vec[1],par_vec[2],n,M, runs=100)
 return(list(E=pk$Ek, V= pk$Vk ))
}

df_100_0 = df_100

df_100 = data.frame(matrix(ncol=2, nrow=100))
for (i in 3:100){
  par_vec = grid[i,]
  pk = Prior_on_K_SB_NGG_1(par_vec[[1]],par_vec[[2]],100,250)
  df_100[i,1]=  pk$Ek
  df_100[i,2]=  pk$Vk
  print(i)
  }
names(df_100)<-c("E","V")
df_100$Std = sqrt(df_100$V)
df_100$beta = unlist(grid[1])
df_100$sigma = unlist(grid[2])
df_100$beta_log = log(unlist(grid[1]))


save(df_100, file = paste0(folderpath,"H1000/df_100_EV_NGG_SB.Rdata"))


########
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
df_r = df_100

for (i in 1:100){
  par_vec = grid[i,]
  df_r[i,3]= par_vec[[1]]
  df_r[i,4]=  par_vec[[2]]
  print(i)
}
names(df_r)<- c("Exp","V","beta","sigma")
df_r$beta_log = log(df_r$beta)
df_r$Std = sqrt(df_r$V)
df_r$beta_scaled= range01(df_r$beta_log)
df_r$color = rgb(df_r$beta_scaled, df_r$sigma,0.5,1)
df_r$color_sigma = rgb(0, df_r$sigma,0.5,1)
df_r$color_beta = rgb(df_r$beta_scaled, 0,0.5,1)
p <- df_r %>% ggplot(aes(x=Exp, y = Std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=3)
p_ngg_sb <- p  + geom_line(data =df_r, aes(x=Exp, y = Std, group=beta, color = color_beta), alpha = 0.8, linetype='dotted') +
  xlim(0, 100)+ ylim(0,20)+ labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()

p_ngg_sb
#pdf(file = '/Users/bystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_SB_100_100.pdf')
#plot(p_py_sb)
#dev.off()

p <- df_r %>% ggplot(aes(x=Exp, y = Std, group=sigma))  + geom_line(aes(color = color_sigma), alpha = 0.8,linetype = 'longdash') + geom_point(aes(colour=color), size=4)
p_sb_ngg <- p  + geom_line(data =df_r, aes(x=Exp, y = Std, group=beta, color = color_beta), alpha = 0.9,linetype = 'dotted') +
  xlim(0, 100)+ ylim(0,20)+labs(y='Std', x='Expectation')+scale_colour_identity()+theme_classic()+ theme(plot.title = element_text(hjust = 0.5,size = 10), axis.text.x = element_text(size=10),legend.position='none')
pdf(file = '/Users/dariabystrova/Documents/GitHub/GibbsTypePriors/test/expectation_variance/Std_variance_NGG_sb.pdf',width= 4, height = 4)
plot(p_sb_ngg)
dev.off()

#############################


p1  <- Prior_on_K_SB_NGG_1(1.0, 0.05,100, 50, runs=2*10^2)
plot(1:100, p1$pk, type="l")

p2  <- Prior_on_K_SB_NGG(1.0, 0.25,Nsamples, 250, runs=2*10^2)
plot(1:100, p1$pk, type="l")


save(p_ngg_1_75, file = paste0(folderpath,"p_ngg_1_75_n1000.Rdata"))


###########


p_ngg_1_25_1000_250 <- Prior_on_K_SB_NGG_1(1.0, 0.25, 1000,1000, runs=2*10^2)
plot(1:1000, p_ngg_1_25_1000_250$pk, type="l")
lines(1:1000, p_ngg_10_25_1000_250$pk, type="l")

p_ngg_10_25_1000_250 <- Prior_on_K_SB_NGG_1(10.0, 0.25, 1000,1000, runs=2*10^2)
plot(1:1000, p_ngg_1_25_1000_250$pk, type="l")