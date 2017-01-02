########################################################################################
## Ex in Chapter 1 in Hoff(2009). A first course in Bayesian statistical methods. 
#                                       Springer
####### Bayesian inference on a proportion #####
# Suppose we are interested in the prevalence of an infectious disease \theta 
# in a small city
# The higher the prevalence, the more public health precautious we would reccomend
# to be put into place 
# Sample: 20 individuals checked for infection
# BEFORE the sample is obtained, the number of infected individuals Y is UNKNOWN
# If \theta were known, a reasonable sampling model for Y is 
# Y | \theta ~ Bin (20, \theta)


######### LIKELIHOOD #############
# What is the "true" distribution for Y ? It depends on the "true" value of \theta!!!

setwd("C:/alessandra/Documents/Didattica/StatisticaBayesiana/mat1617")
n<-20
x<-0:n
del<-.25
plot( range(x-del), c(0,.4),xlab="number infected in the sample",
     ylab="probability",type="n")

points( x-del,dbinom(x,n,.05),type="h",col=gray(.75),lwd=3)
points( x,dbinom(x,n,.10),type="h",col=gray(.5),lwd=3)
points( x+del,dbinom(x,n,.20),type="h",col='blue',lwd=3)
legend(10,.35,legend=c(
    expression(paste(theta,"=0.05",sep="")), 
    expression(paste(theta,"=0.10",sep="")),
    expression(paste(theta,"=0.20",sep="")) ),
     lwd=c(3,3,3), 
    col=c(gray(c(.75,.5)),'blue') ,bty="n") 

# Under these assumptions, what is the probability that there is no infected invividual in the sample?
# It depends on theta!
# If theta=5%, the probability that there is no infected individual in the sample is 
dbinom(0,20,0.05)

dbinom(0,20,0.1) # Se theta=10%

dbinom(0,20,0.2) # Se theta=20%

################ PRIOR DISTRIBUTION ################
# Other studies from various parts of the country indicate that the infection 
# rate in comparable cities ranges from about 0.05 to 0.2
# with an average prevalence of 0.1
# We have to find a prior distribution on \theta consistent with this information:
# a prior pi(\theta) that assigns a substantial amount of prob to (0.05,0.2)
# and with expected value close to 0.1 
# There are infinite priors consistent with these 2 conditions !
# But for some reasons (wait a couple of classes!) we consider Beta distributions
# A priori theta \sim Beta(a,b)

a<-2 ; b<-20
a/(a+b)  #media a priori di theta
pbeta(.20,a,b) - pbeta(.05,a,b)
pbeta(.10,a,b)

(a-1)/(a-1+b-1) #moda a priori di theta (se a, b>1)

###################################################
##   CONFRONTO tra varie distribuzioni BETA      ##
###################################################
p=seq(0,1,length=500)
plot(p,dbeta(p,0.5,0.5),xlab=" ", ylab=" ",type="l", lwd=2, ylim=c(0,8),
         main="Confronto fra densità  beta")
lines(p,dbeta(p,2,3),ylim=c(0,8),col="red", lwd=2)
lines(p,dbeta(p, a,b), ylim=c(0,8),col=3,lwd=2)

X11()
########################PRIOR for this EXAMPLE ############################## 
p=seq(0,1,length=500)
plot(p,dbeta(p,a,b),xlab=" ", ylab=" ",type="l",ylim=c(0,8),main="PRIOR distribution",lwd=2)

dev.off()

################## OBSERVED DATA #####################
y<-0 ; n<-20 # Non osservo NESSUN infetto nel campione a disposizione!!!

################ POSTERIOR DISTRIBUTION ################
# a posteriori theta \sim Beta(a+y,b+n-y)
# 
# Rappresenta la nostra "incertezza" sulla proporzione di infetti nella cittadina 
# considerata, alla luce del dato y=0 (NESSUN infetto nella cittadina su 20 selezionati)

a<-2 ; b<-20


(a+y)/(a+b+n)    #media a posteriori
(a+y-1)/(a-1+b+n-1)  #moda a posteriori
pbeta(.20,a+y,b+n-y) - pbeta(.05,a+y,b+n-y)
pbeta(.10,a+y,b+n-y)


######### GRAFICO di PRIOR e POSTERIOR insieme
theta<-seq(0,1,length=500)
plot(theta, dbeta(theta,a+y,b+n-y),
     type="l",
     xlab="percentage infected in the population",
     ylab="", lwd=2, ylim=c(0,16), xlim=c(0,0.5)
     )
lines(theta, dbeta(theta,a,b),col="gray",lwd=2)
legend(.3,14,legend=c( expression(paste(italic("p"),"(",theta,")",sep="")), 
     expression(paste(italic("p"),"(",theta,"|",italic("y"),")",sep=""))  ), 
  bty="n", lwd=c(2,2),col=c("gray","black"))


windows()

#### Grafico di Verosimiglianza e Prior+Posterior
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

n<-20
x<-0:n
del<-.25

plot( range(x-del), c(0,.4),xlab="number infected in the sample",
     ylab="probability",type="n")

points( x-del,dbinom(x,n,.05),type="h",col=gray(.75),lwd=3)
points( x,dbinom(x,n,.10),type="h",col=gray(.5),lwd=3)
points( x+del,dbinom(x,n,.20),type="h",col='blue',lwd=3)
legend(10,.35,legend=c(
    expression(paste(theta,"=0.05",sep="")), 
    expression(paste(theta,"=0.10",sep="")),
    expression(paste(theta,"=0.20",sep="")) ),
     lwd=c(3,3,3), 
    col=c(gray(c(.75,.5)),'blue') ,bty="n") 

a<-2 ; b<-20
y<-0 ; n<-20

theta<-seq(0,1,length=500)
plot(theta, dbeta(theta,a+y,b+n-y),
     type="l",
     xlab="percentage infected in the population",
     ylab="", lwd=2, ylim=c(0,16)
     )
lines(theta, dbeta(theta,a,b),col="gray",lwd=2)
legend(.5,14,legend=c( expression(paste(italic("p"),"(",theta,")",sep="")), 
     expression(paste(italic("p"),"(",theta,"|",italic("y"),")",sep=""))  ), 
  bty="n", lwd=c(2,2),col=c("gray","black"))


dev.off()

################# BAYESIAN LEARNING ###################### 
a/(a+b)  #prior mean of theta
(a+y)/(b+n-y) #posterior mean of thera

(a-1)/(a-1+b-1) #prior mode of theta (se a, b>1)
(a+y-1)/(a+y-1+b+n-y-1) #posterior mode

pbeta(.20,a,b) - pbeta(.05,a,b)
pbeta(.20,a+y,b+n-y) - pbeta(.05,a+y,b+n-y)

pbeta(.10,a,b)
pbeta(.10,a+y,b+n-y)

# NOTATE che a posteriori 
# E(theta|y)=(a+y)/(a+b+n)= n /(a+b+n)*(y/n) + (a+b)/(a+b+n)*(a/(a+b))



########################### INFERENZA CLASSICA ##############################
# STIMA puntuale di theta è y/n, ma in questo caso y=0, e quindi
# la stima puntuale di theta è uguale a 0: NON ha senso!
# STIMA INTERVALLARE si riduce ad un unico punto (0): NON ha senso!

#### Adjusted Wald interval 
#    th=n/(n+4)* (y/n) + 4/(n+4)* (1/2)

a<-2 ; b<-2 
th<-  (y+a)/(n+a+b) #th=n/(n+4)* (y/n) + 4/(n+4)* (1/2)
th  # th corrisponde alla MEDIA a POSTERIORI di theta quando la prior è una Beta(2,2)!!

th+c(-1,1)*1.96*sqrt(th*(1-th)/n) # stima intervallare





