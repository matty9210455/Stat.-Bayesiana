############################################
##    A BAYESIAN HIERARCHICAL Model       ##
##    Ex. 8.4 from P. Hoff's book         ##
############################################
# Y_{i,j} ~ N(theta_j, sigma^2) student i in school j, i=1,...,n_j,
#                                                      j=1,...,J
# theta_j ~ N(mu,tau^2)

# Multilevel data
setwd("C:/alessandra/Documents/Didattica/StatisticaBayesiana/mat1516")
Y=read.table("school.mathscore.txt") 

dim(Y)
## Scores at the math test for 1993 students in 100 US schools 
## First column: index of the school of the student
## Second column: math score of the student


### group means of the math scores
m = length(unique(Y[,1]))
n<-sv<-ybar<-rep(NA,m) 
for(j in 1:m) 
{ 
  ybar[j]<- mean(Y[Y[,1]==unique(Y[,1])[j],2])   #medie per gruppo
  sv[j]<-var(Y[Y[,1]==unique(Y[,1])[j],2])  #varianze empiriche in ogni gruppo
  n[j]<-sum(Y[,1]==unique(Y[,1]))        #numerosit? in ogni gruppo (scuola)
}

#summary(ybar)

ybar
n
sv

############ Grafico dei dati ###########
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
hist(ybar,main="",xlab="sample mean")
plot(n,ybar,xlab="sample size",ylab="sample mean")
# LEFT: histogram of the empirical means per school ybar_j
# RIGHT: very extreme sample averages (very small or very large)  tend to be associated 
#       with schools with small sample size
# Dal grafico a destra, sembra che le scuole con ybar pi? grande (sample average) o pi? piccolo
#  corrispondano a numerosit? piccole
# Questo si spiega perch?, a parit? di media campionaria, la varianza campionaria ? 
# \sigma^2/n_i: pi? ? picolo n_i e maggiore sar? la variabilit? di ybar.


############################# PARAMETRI ##########################################
# PARAMETERS: (theta_1,...,theta_100, mu, tau^2,sigma^2)
# theta_j=mean in school j, mu= grand mean,
# tau^2 = variance between schools, sigma^2 = variance within school j (costant)
# la prior per theta_1,...,theta_100  ? di tipo gerarchico - vedi i lucidi

#################### PRIOR for(sigma^2, mu , tau^2)############################
### 1/sigma^2 \sim gamma(nu0/2, (nu0*sigma0^2)/2)
### 1/tau^2 \sim gamma(eta0/2, (eta0*tau0^2)/2)
### mu \sim N(mu0,gamma0^2)
#  sigma^2: within-group variance (it is constant in this case, but a more general model
#            with sigma_j^2)
#  tau^2: between-group variance

x11()
##### IPERPARAMETRI assegnati sulla base di INFORMAZIONI a priori, cio? a prescindere
#####   dai dati di questo specifico esperimento
par(mfrow=c(1,1))
nu0<-1  ; s20<-100   # iperparametri per la prior di sigma^2 = within-group variance
pippo=rgamma(1000,shape=nu0/2, rate=nu0/2*s20 )
hist(pippo,prob=T,main='Prior per 1/sigma^2')
mean(1/pippo);var(1/pippo)  # media e varianza di sigma^2 sono +infinito! 
eta0<-1 ; t20<-100 # iperparametri per la prior di tau^2
mu0<-50 ; g20<-25  # IMPORTANT: iperparametri di mu 
#### A priori P(\mu \in (40,60))= 0.94 circa e E(\mu)=50 che ? la media nazionale
#### Prior informativa, ricavata da info sul test somministato a tutta la nazione
dev.off()

######## MCMC analysis for school data ########
# PARAMETERS: (theta_1,...,theta_100, mu, tau^2,sigma^2)
# theta_j=mean in school j, mu= grand mean,
# tau^2 = variance between schools, sigma^2 = variance within school j (costant)
##### starting values of the MC
theta<-ybar
sigma2<-mean(sv)
mu<-mean(theta)
tau2<-var(theta)
###

### setup MCMC
set.seed(1)
S<-3000
#S<-5
THETA<-matrix( nrow=S,ncol=m)
MST<-matrix( nrow=S,ncol=3)
###

### MCMC algorithm - GIBBS SAMPLER
for(s in 1:S) 
{

  # sample new values of the thetas
  for(j in 1:m) 
  {
    vtheta<-1/(n[j]/sigma2+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
   }

  #sample new value of sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20;for(j in 1:m){ss<-ss+sum((Y[Y[,1]==j,2 ]-theta[j])^2)}
  sigma2<-1/rgamma(1,nun/2,ss/2)

  #sample a new value of mu
  vmu<- 1/(m/tau2+1/g20)
  emu<- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu)) 

  # sample a new value of tau2
  etam<-eta0+m
  ss<- eta0*t20 + sum( (theta-mu)^2 )
  tau2<-1/rgamma(1,etam/2,ss/2)

  #store results
  THETA[s,]<-theta
  MST[s,]<-c(mu,sigma2,tau2)

} 

### THETA contiene i valori simulati dei 100 theta_i, MST contiene i valori simulati di 
### mu, sigma^2, tau^2
###

mcmc1<-list(THETA=THETA,MST=MST)

stationarity.plot<-function(x,...){
S<-length(x)
scan<-1:S
ng<-min( round(S/100),10)
group<-S*ceiling( ng*scan/S) /ng
boxplot(x~group,...)               
}

### Check CONVERGENCE of the MC
# Each graph contains 10 boxplots of 500 next iterations
# (dunque 1/10 di tutti i dati simulati) of mu,sigma^2 e tau^2 
# Boxplots are similar: convergence is OK
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(MST[,1],xlab="iteration",ylab=expression(mu))
stationarity.plot(MST[,2],xlab="iteration",ylab=expression(sigma^2))
stationarity.plot(MST[,3],xlab="iteration",ylab=expression(tau^2))

library(coda)
effectiveSize(MST)
# effective sample size of MST are good!
# Good autocorrelatio functions
par(mfrow=c(1,3))
acf(MST[,1]) -> a1
acf(MST[,2]) -> a2
acf(MST[,3]) -> a3

# Stima dell'errore Monte Carlo = posterior sd / sqrt(effectivesamplesize) 
#        per mu, sigma^2, tau^2
MCERR<-  apply(MST,2,sd)/sqrt( effectiveSize(MST) )
MCERR

### POSTERIOR MEANS of  mu, sigma^2, tau^2
apply(MST,2,mean)
# Questi ultimi sono molto pi? grandi di MCERR

100*MCERR/apply(MST,2,mean) #sono tutti e 3 minori di 1%, cio? MCerr diviso la media a posteriori ? piccolo


# effective sample size of theta_j are OK!
effectiveSize(THETA) -> esTHETA
esTHETA

# Stima dell'errore Monte Carlo dei THETA = posterior sd / sqrt(effectivesamplesize) 
TMCERR<-  apply(THETA,2,sd)/sqrt( effectiveSize(THETA) )
TMCERR


##### MARGINAL POSTERIOR densities of mu, sigma^2, tau^2
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(density(MST[,1],adj=2),xlab=expression(mu),main="",lwd=2,
ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,1],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,2],adj=2),xlab=expression(sigma^2),main="", lwd=2,
ylab=expression(paste(italic("p("),sigma^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,2],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )
plot(density(MST[,3],adj=2),xlab=expression(tau^2),main="",lwd=2,
ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,3],c(.025,.5,.975)),col="gray",lty=c(3,2,3) )


##### POSTERIOR MEANS of  mu, sigma, tau
mean((MST[,1]))
mean(sqrt(MST[,2]))
mean(sqrt(MST[,3]))


##################################################
##### SHRINKAGE effect towards grand mean mu #####
##################################################
### Vediamo come le informazioni sono state condivise tra i diversi gruppi: 
###                   Bayesian hierarchical approach
### The Bayesian estimate of group-specific parameter theta_j is a convex linear combination
### of the frequentist group estimate ybar_j on one hand, and of mu (the mean of all theta_j); 
### the Bayesian estimate is pulled a bit from ybar_j towards mu by an amount depending on n_j
### This is called SHRINKAGE 
### BTW, mu is random, but E(mu)=50
### Se n_j ? GRANDE, la media empirica ? una buona stima, anche da punto di vista bayesiano;
### dunque non c'? necessit? ci CHIEDERE in PRESTITO (to BORROW) informazioni dal resto dei gruppi.
### Se n_j ? piccolo, ybar_j NON va bene come stima e "correggo" in modo tale che la stima bayesiana 
### si avvicini alla (stima) di mu, la media dei gruppi
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
#################################### theta.hat= posterior means of theta_j
theta.hat<-apply(THETA,2,mean)
## LEFT
## ybar (le stime 'classiche', che sono sample average per scuola) on the x-axis and
## theta.hat (le stime bayesiane) on the y-axis
plot(ybar,theta.hat,xlab=expression(bar(italic(y))),ylab=expression(hat(theta)))
abline(0,1)
## The slope of this line is <1, that is, high values of ybar.j correspond to
## slightly less high values of the Bayesian estimates of theta_j,
## and low values of ybar.j correspond to slightly less low values of the Bayesian 
## estimates of theta_j. This is the SHRINKAGE effect

## RIGHT
## group-specific sample sizes on the x-axis, and differences between frequentist and 
## Bayesian estimates on the y-axis.
plot(n,ybar-theta.hat,ylab=expression( bar(italic(y))-hat(theta) ),xlab="sample size")
abline(h=0)
## Groups with low sample size get shrunk the most, whereas groups with large sample size
## hardly get shrunk at all. 
## The larger the sample size for a group, the more information we have for that group
## and  the less information we need to BORROW from the rest of the population.
#dev.off()
#####

x11()
################# SCHOOLS COMPARISON ##########################
## posterior CIs for each theta_j, group-specific parameters
plot(1:100,theta.hat, cex=0.5,main="90% credible intervals",xlim=c(1,100),ylim=c(35,70),xlab='school',ylab='')
for (i in 1:100) {
probint=quantile(THETA[,i], c(0.05, 0.95))
lines(i * c(1, 1), probint)
}

# min of theta.hat is SCHOOL 5
# max of theta.hat is SCHOOL 51

# For each couple (theta_j, theta_l): posterior probability that theta_j > theta_l, that is 
# posterior prob that school j is BETTER than school l
# Matrix of these posterior probabilities, better

compare.rates <- function(x) {
  nc <- 100
  ij <- as.matrix(expand.grid(1:nc, 1:nc))
  m <- as.matrix(x[,ij[,1]] > x[,ij[,2]]) 
  matrix(colMeans(m), nc, nc, byrow = TRUE)
}

better=compare.rates(THETA)
better

# Posterior probability that school 5 is better than all the other schools
better[,5]
windows()
plot(better[,5])
# All are <=0.37
min(better[-5,5]);max(better[-5,5])
better[51,5] # probabilit? a posteriori che la scuola 5 sia meglio della scuola 51


# Posterior probability that school 51 is better than all the other schools 
better[,51]
windows()
plot(better[,51])
# All are >= 0.8398 and <=1
min(better[-51,51]);max(better[-51,51])
better[5,51] # probabilit? a posteriori che la scuola 51 sia meglio della scuola 5


better[,46]
windows()
plot(better[,46])
min(better[-46,46]);max(better[-46,46])

# According to the group sample averages, school 46 is better that school 82
ybar[c(46,82)]
# Posterior probability that school 46 is better than 82
better[82,46]

