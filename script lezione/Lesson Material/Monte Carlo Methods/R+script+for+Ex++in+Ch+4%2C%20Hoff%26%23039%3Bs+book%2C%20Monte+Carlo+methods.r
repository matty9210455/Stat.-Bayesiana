##############################################
##          MONTE CARLO METHODS             ##
## From Chapter 4 in P. Hoff's book         ##
##############################################

# Y_1 is the vector of the numbers of children for the n_1 women WITHOUT college degrees
# Y_2 is the vector of the numbers of children for the n_2 women WITH college degrees
# The components of the two vectors are iid from Poisson with parameter theta_1 and theta_2, respectively
# PRIOR for (theta_1,theta_2: theta_1,theta_2 iid from gamma(alpha,beta)
# alpha=2, beta=1, in such a way that E(theta_i)=2, Var(theta_i)=2
# POSTERIOR of (theta_1,theta_2): gamma(a+\sum y_{1,i},beta+n_1)\times gamma(a+\sum y_{2,i},beta+n_2)
#                                =gamma(219,112)\times gamma(68,45)
# NOTATION: MC standard error = sqrt(variance of the estimate of the objective parameter) 

##Random number draws from the posterior of  theta_2, given Y_{1,2},Y_{2,2},...Y_{n_2,2}
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(2,3))
set.seed(1)
a<-68 ; b<-45
set.seed(1)
theta.support<-seq(0,3,length=100)
## 3 different groups of iid draws from the posterior of theta_2
theta.sim10<-rgamma(10,a,b)# iid sample of size M=10 from the gamma(a,b)
theta.sim100<-rgamma(100,a,b) #size M=100
theta.sim1000<-rgamma(1000,a,b) #size M=1000

xlim<-c(.75,2.25)
ylim=c(0,2.5)
lty=1

#Plot of the histogramm and kernel density estimate for the 3 groups of iid draws
hist( theta.sim10, prob=T,xlim=xlim,ylim=ylim,xlab="",main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(M),"=10",sep="")))

hist( theta.sim100, prob=T,xlim=xlim,ylim=ylim,xlab="",main="" ,ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(M),"=100",sep="")))

hist( theta.sim1000, prob=T,xlim=xlim,ylim=ylim,xlab="",main="" ,ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)
text(2.1,2.25,expression(paste(italic(M),"=1000",sep="")))


plot(density(theta.sim10),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)

plot(density(theta.sim100),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)

plot(density(theta.sim1000),xlim=xlim,ylim=ylim,xlab=expression(theta),main="",ylab="")
lines(theta.support,dgamma(theta.support,a,b),col="gray",lwd=2,lty=lty)

dev.off()



#Kernel density estimator
?density
pippo=density(theta.sim1000)
plot(density(theta.sim1000),main="",ylab="",lwd=2)
#windows()
lines(density(theta.sim1000,bw=0.01),main="",ylab="",col="blue",lwd=2)
dev.off()


##### Monte Carlo computation of the posterior mean of  theta_2
## A posteriori: theta_2 ~ gamma(a+sy=68, b+n=45)

set.seed(1)
a<-2  ; b<-1
sy<-66; n<-44

theta.sim10<-rgamma(10,a+sy,b+n)
theta.sim100<-rgamma(100,a+sy,b+n)
theta.sim1000<-rgamma(1000,a+sy,b+n)

(a+sy)/(b+n)# exact value esatto 

mean(theta.sim10) # MC estimate and MC standard error for 10 draws
sqrt(var(theta.sim10)/10)
# Probability that the absolut value of the error is larger than c 
c=0.01 
2*(1-pnorm(c/(sqrt(var(theta.sim10)/10))))

mean(theta.sim100)  #  MC estimate and MC standard error for 100 draws
sqrt(var(theta.sim100)/100)
2*(1-pnorm(c/(sqrt(var(theta.sim100)/100))))


mean(theta.sim1000)  #  MC estimate and MC standard error for 1000 draws
sqrt(var(theta.sim1000)/1000)
2*(1-pnorm(c/(sqrt(var(theta.sim1000)/1000))))

# MC computation of the posterior distribution function at 1.75
pgamma(1.75,a+sy,b+n) #valore esatto

# The MC estimate is the relative frequency of the simulated draws less than 1.75 
mean( theta.sim10<1.75)   # I'm using 10 MC draws 
mean( theta.sim100<1.75)  # 100 draws
mean( theta.sim1000<1.75) # 1000 draws

# Posterior 95% credible interval, given by two symmetric quantiles
qgamma(c(.025,.975),a+sy,b+n)
quantile( theta.sim10, c(.025,.975))
quantile( theta.sim100, c(.025,.975))
quantile( theta.sim1000, c(.025,.975))


###### Monte Carlo approximations as the sample size increases 
#pdf("fig4_2.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))

set.seed(1)
a<-2   ; b<-1
sy<-66 ; n<-44

nsim<-1000
theta.sim<-rgamma(nsim,a+sy,b+n)

#cumulative mean
#?cumsum
cmean<-cumsum(theta.sim)/(1:nsim)
cvar<- cumsum(theta.sim^2)/(1:nsim) - cmean^2
ccdf<- cumsum(theta.sim<1.75)/ (1:nsim)
cq<-NULL
for(j in 1:nsim){ cq<-c(cq,quantile(theta.sim[1:j],probs=0.975)) }

sseq<- c(1,(1:100)*(nsim/100))
cmean<-cmean[sseq] 
cq<-cq[sseq] 
ccdf<-ccdf[sseq] 

# Plots ofthe MC estimates of the mean, the df at 1.75, the CI as the sample size M increases
# REM: there is no monotonicity!
plot(sseq,cmean,type="l",xlab="# of Monte Carlo samples",ylab="cumulative mean",
     col="black")
abline(h= (a+sy)/(b+n),col="gray",lwd=2)

plot(sseq,ccdf,type="l",xlab="# of Monte Carlo samples",ylab="cumulative cdf at 1.75",col="black")
abline(h= pgamma(1.75,a+sy,b+n),col="gray",lwd=2)

plot(sseq,cq,type="l",xlab="# of Monte Carlo samples",ylab="cumulative 97.5% quantile",col="black")
abline(h= qgamma(.975,a+sy,b+n),col="gray",lwd=2)

dev.off()

########## Closer look at the MC computation of the posterior mean of theta2
a<-2   ; b<-1
sy<-66 ; n<-44

nsim<-1000
theta.sim<-rgamma(nsim,a+sy,b+n)
sseq<- c(1,(1:100)*(nsim/100))

cmean<-cumsum(theta.sim)/(1:nsim)
cmean<-cmean[sseq] 

plot(sseq,cmean,type="l",xlab="# of Monte Carlo samples",ylab="cumulative mean",
     col="black")
abline(h= (a+sy)/(b+n),col="gray",lwd=2)

dev.off()

#######################################################################
##### MC computation of the posterior probability that theta1>theta2 
set.seed(1)
a<-2 ; b<-1
sy1<-217 ;  n1<-111
sy2<-66  ;  n2<-44

a+sy1; b+n1
a+sy2; b+n2

# MC sample of size M=10000 from the joint posterior
theta1.mc<-rgamma(10000,a+sy1, b+n1)
theta2.mc<-rgamma(10000,a+sy2, b+n2)

# The probability is estimated by the number of times in the sequence {(theta_1^(j),theta_2^(j)),j=1,2,...,M}
# the first component is > than the second, divided by M
mean(theta1.mc>theta2.mc) # questa probabilit? ? uguale alla prob che theta1 sia maggiore o UGUALE a theta2
                          #  visto che (theta1,theta2) ha distribuzione assolutamente continua
  

## As an alternative, I can compute the posterior distribution of theta1/theta2;
## in yhis case, the analytic expression of its posterior is available 

#pdf("fig4_4.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))
plot(density(theta1.mc/theta2.mc,adj=2),main="",xlim=c(.75,2.25),
xlab=expression(gamma==theta[1]/theta[2]),
ylab=expression(paste(italic("p("),gamma,"|",bold(y[1]),",",bold(y[2]),")",
   sep="")) )


dev.off()

############# COMPUTATION of PREDICTIVE DISTRIBUTIONS 
# theta1.mc is a vector of size 10thousant from teh posterior distribution of theta1 
y1.mc<-rpois(10000,theta1.mc) #genero in un "colpo solo" un valore dalla Pois(theta1_j), j=1...10mila
y2.mc<-rpois(10000,theta2.mc)

mean(y1.mc>y2.mc)  #it is about 0.48

mean(y1.mc==y2.mc)  # it's about 0.21 ; 
                    #### hence the posterior predictive prob. that Y_1 is larger or equal to Y_2 is 0.69

mean(y1.mc<y2.mc)   #? circa 0.31


#pdf("fig4_5.pdf",family="Times",height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))

diff.mc<-y1.mc-y2.mc

ds<- -11:11
plot(ds,(table(c(diff.mc,ds))-1)/length(diff), type="h",lwd=3,
    xlab=expression(italic(D==tilde(Y)[1]-tilde(Y)[2])),
  ylab=expression(paste(italic("p(D"),"|",bold(y[1]),",",bold(y[2]),")",sep="")))


##### marginali a posteriori di theta1 e di theta2
x11()
par(mfrow=c(1,2))
plot(theta1.mc,theta2.mc)
plot(density(theta1.mc),xlim=c(0.5,3.5),lwd=2,main='posterior di theta1 e theta2' )
lines(density(theta2.mc),main="",xlim=c(0.5,3.5),lwd=2, col=2)
legend(2,2.5,legend=c(
    expression(paste(theta,"1",sep="")), 
    expression(paste(theta,"2",sep="")) ),
     lwd=c(2,2), 
    col=c('black','red') ,bty="n") 


#### marginali predittive esatte - bionomiali negative
par(mfrow=c(1,2))
a+sy1; (b+n1)/(1+b+n1)  #parametri della predittiva marginale di Y_1
a+sy2; (b+n2)/(1+b+n2)  #parametri della predittiva marginale di Y_2


ds=seq(0,10)
plot(ds,dnbinom(ds,a+sy1, (b+n1)/(1+b+n1)), type="h",lwd=3,xlab='y', ylab='PredMargY_1',ylim=c(0,0.35))
plot(ds,dnbinom(ds,a+sy2, (b+n2)/(1+b+n2)), type="h",lwd=3,xlab='y', ylab='PredMargY_1',ylim=c(0,0.35))


dev.off()

#################################  THE END  ##################################### 


