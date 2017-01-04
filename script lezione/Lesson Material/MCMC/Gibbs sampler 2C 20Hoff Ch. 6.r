#################################################
########### Chapter 6 - Hoff's book #############
#################################################
# The data are wing lenght (in mm) of n=9 midges (moscerini)
# The model : Y_i are conditionally iid N(theta,sigma^2)
# Both theta and tau=1/sigma^2 are unknown
# Prior: pi(theta,tau)= pi(theta)* pi(tau)
#        theta ~ N(mu0, t20) (mu0 is the mean, t20 is the variance)
#        tau ~ gamma(nu0/2, (nu0*s20)/2 )

mu0<-1.9  ; t20<-0.95^2
s20<-.01 ; nu0<-1

1/s20

#data
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n<-length(y) ; mean.y<-mean(y) ; var.y<-var(y)


### Generate a Gibbs sampler MCMC of S draws
set.seed(1)
S<-1000
PHI<-matrix(nrow=S,ncol=2) #Inizializzo una matrice con 2 colonne e S righe
                           # it contains simulated values of the BIVARIATE Gibbs sampler MC 
PHI[1,]<-phi<-c( mean.y, 1/var.y) # Initial point of the chain

################################
####### Gibbs sampling #########

for(s in 2:S) {
# the bidim. vector phi contains the "current" value of the chain
# generate a new theta value from its full conditional
mun<-  ( mu0/t20 + n*mean.y*phi[2] ) / ( 1/t20 + n*phi[2] ) #updated mean 
t2n<- 1/( 1/t20 + n*phi[2] ) #updated variance
phi[1]<-rnorm(1, mun, sqrt(t2n) )

# generate a new tau=1/sigma^2 value from its full conditional
nun<- nu0+n             # it does NOT change with the iteration
s2n<- (nu0*s20 + (n-1)*var.y + n*(mean.y-phi[1])^2 ) /nun #updated parameter
phi[2]<- rgamma(1, nun/2, nun*s2n/2)

PHI[s,]<-phi         }
###

PHI

### First iterations of Gibbs sampler
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-1
plot( PHI[1:m1,],type="p",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-2
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-3
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )
###

x11()
#pdf("fig6_2.pdf",height=1.75,width=5,family="Times")

par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-5
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-15
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-100
plot( PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
       lty=1,col="gray",xlab=expression(theta),ylab=expression(tau^2))
text(  PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

dev.off()
#####

sseq<-1:1000


plot(PHI[sseq,1],PHI[sseq,2],xlab=expression(theta), ylab=expression(tau^2) ,
     xlim=range(PHI[,1]),ylim=range(PHI[,2]))


##### PLOTS of the posterior marginal densities, with credible intervals
par(mfrow=c(1,2),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
sseq<-1:1000

#image( mean.grid,prec.grid,post.grid,col=gray( (10:0)/10 ),
 #    xlab=expression(theta), ylab=expression(tilde(sigma)^2) ,xlim=range(PHI[,1]),ylim=range(PHI[,2]) )


plot(density(PHI[,1],adj=2),  xlab=expression(theta),main="",
     xlim=c(1.55,2.05),
 ylab=expression( paste(italic("p("),
     theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
abline(v=quantile(PHI[,1],prob=c(.025,.975)),lwd=2,col="gray")


plot(density(PHI[,2],adj=2), xlab=expression(tau^2),main="",
     ylab=expression( paste(italic("p("),
     tau^2,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 
abline(v=quantile(PHI[,2],prob=c(.025,.975)),lwd=2,col="gray")


dev.off()
#####

quantile(PHI[,1],c(.025,.5,.975)) # Posterior CI of theta
quantile(PHI[,2],c(.025,.5, .975)) # Posterior CI of tau
quantile(1/sqrt(PHI[,2]),c(.025,.5, .975))  # Posterior CI of sigma


###################### CONVERGENCE DIAGNOSTICS #################

#####
#pdf("fig6_7.pdf",family="Times",height=3.5,width=7)

par(mfrow=c(1,2))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(PHI[,1],xlab="iteration",ylab=expression(theta))  
plot(1/PHI[,2],xlab="iteration",ylab=expression(sigma^2))  

dev.off()

## marginal traceplots
par(mfrow=c(1,2))
plot(ts(PHI[,1]),xlab="iteration",ylab=expression(theta))
plot(ts(1/PHI[,2]),xlab="iteration",ylab=expression(sigma^2))


library(coda)

par(mfrow=c(1,2))
acf(PHI[,1]) -> tmp1
acf(1/PHI[,2]) -> tmp2


# The Effective Sample Size (ESS) of a one-dimensional component of a MCMC is
# the number of INDEPENDENT iterations to achieve the same MCerror from 
# the same marginal distribution of the MC (the stationary dist).
# For a time series x of length N, the standard error of the mean is var(x)/n where n is 
# the effective sample size. n = N only when there is no autocorrelation.
# If the ESS is small (compared to the sample size of the chain), this means that the  
# MC estimate will be rather 'poor', because the variance of the MC estimate will be 'large'
# wrt the variance of the Monte Carlo estimator. 
# The larger the ESS, the better! 

effectiveSize( PHI[,1] )
effectiveSize(1/PHI[,2] )

##########################################################################################