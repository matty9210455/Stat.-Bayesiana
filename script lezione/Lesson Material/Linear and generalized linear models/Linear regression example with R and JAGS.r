#############################################
## A (CONJUGATE) linear regression example ##
## Example 2.15 - S. Jackman (2009)        ##
#############################################
## In November 1993, the state of Pennsylvania conducted elections for its state legislature. 
## The result in the Senate election in the 2nd district (based in Philadelphia) was challenged 
## in court, and ultimately overturned. The Democratic candidate won 19,127 of the votes cast 
## by voting machine, while the Republican won 19,691 votes cast by voting machine, giving 
## the Republican a lead of 564 votes. However, the Democrat won 1,396 absentee ballots, while 
## the Republican won just 371 absentee ballots, more than offsetting the Republican lead based 
## on the votes recorded by machines on election day. 
## The Republican candidate sued, claiming that many of the absentee ballots were fraudulent. 
## The judge in the case solicited expert analysis from Orley Ashenfelter, an economist at 
## Princeton University. Ashenfelter examined the relationship between absentee vote margins 
## and machine vote margins in 21 previous Pennsylvania Senate elections in seven districts
## in the Philadelphia area over the preceding decade.
## OUTLIER detection problem in the contect of simple linear regression 
## in the Pennsylvania state senate election
## Response Y= abseentee vote margins (difference of Democrat percentage and Republican percentage) in the absentee ballot
## Absentee ballot = VOTO A DISTANZA
## See Jackman, from p. 87 
## COVARIATE x = machine vote margins (difference of Democrat percentage and Republican percentage)

## R package pscl where you can find the dataset
library(pscl)
data(absentee)
help(absentee) # A data frame with 22 observations on 8 variables.
summary(absentee)

attach(absentee)

## create variables for regression analysis
y <- (absdem - absrep)/(absdem + absrep)*100
x <- (machdem - machrep)/(machdem + machrep)*100
x[22]
y[22]

# The 22nd datapoint is suspect! 
## DATA PLOT ##
plot(y~x,type="n",xlim=range(x),ylim=range(y),
     xlab="Democratic Margin, Machine Ballots (Percentage Points)",
     ylab="Democratic Margin, Absentee Ballots (Percentage Points)")

ols <- lm(y ~ x, subset=c(rep(TRUE,21),FALSE)) ## OLS analysis of the dataset, but DROP data point 22

abline(ols,lwd=2) ## overlay ols (retta di regressione stimata, senza il dato 22)
points(x[-22],y[-22],pch=1) ## data
points(x[22],y[22],pch=16) ## disputed data point
text(x[22],y[22],"Disputed\nElection",cex=.75,adj=1.25)
axis(1)
axis(2)


summary(ols)


############## DATA ANALYSIS WITHOUT the 22nd observation
y=y[1:21]
x=x[1:21]
X=as.matrix(cbind(rep(1,21),x))
n=length(y)
k=dim(X)[2] #lo abbiamo chiamato p alla lavagna!

######## PARAMETERS
## beta0, beta1, sigma^2

#### ORDINARY LEAST SQUARE ESTIMATION OF beta

inv=function(X)
{
# RETURN THE INVERSE OF THE SYMMETRIC MATRIX X
EV=eigen(X)
EV$vector%*%diag(1/EV$values)%*%t(EV$vector)
}

inv(t(X)%*%X); solve(t(X)%*%X)#use either these commands to compute the inverse of a square matrix

betahat=inv(t(X)%*%X)%*%t(X)%*%y
betahat

#### FREQUENTIST UNBIASED ESTIMATION OF sigma2
S2=t(y-X%*%betahat)%*%(y-X%*%betahat) ## SOMMA dei RESIDUI al quadrato
sigma2hat=S2/(n-k)
sigma2hat
sqrt(sigma2hat)

#### ESTIMATION OF THE VARIANCE OF betahat
diag(as.real(sigma2hat)*inv(t(X)%*%X))
#### ESTIMATION OF THE sd OF betahat
sqrt(diag(as.real(sigma2hat)*inv(t(X)%*%X)))


################### BAYESIAN ANALYSIS #############################
## CONJUGATE PRIOR: beta|sigma^2 ~ N(b0, sigma^2*B0)
##        sigma^2 ~ inv-gamma (nu0/2, nu0*(sigma0^2)/2)
##                                tau0=1/sigma0^2
## LIKELIHOOD: Gaussian data  Y ~ N_n()

## b0 = (0,1), i.e. I believe that  beta0=0 and beta1=1 on average a priori,
## equivalent to E(Y_i|x_i)=x_i
## M=B0^{-1}

nu0=6.2
tau0=1/47.07

b0=c(0,1)
B0= diag(c(2.25,0.0022),2,2) #Si veda Jackman (2009) p. 108 e seguenti per 
                             # la scelta degli iperparametri in B0
# A PRIORI, we assume that beta0 and beta1 are UNCORRELATED, i.e. independent. 
(B0/tau0)*(nu0/(nu0-2)) # marginal priopr variance of vector beta
b0[2]-3*sqrt((B0[2,2]/tau0)*(nu0/(nu0-2)));b0[2]+3*sqrt((B0[2,2]/tau0)*(nu0/(nu0-2)))
# equivale ad avere assunto che a priori P(beta1<0)=0.01, cio? che P(beta1>0)=0.99
# per beta0, stiamo assumendo che in (-25,25) ci sia il 95% della massa a priori 

nu0=6.2
tau0=1/47.07

a= nu0/2
b= nu0/(2*tau0)

nu0;tau0
a;b

# PRIOR MEAN OF sigma2
b/(a-1)

# PRIOR VARIANCE OF sigma2
b^2/(a-1)^2/(a-2) #sigma2 ? molto dispersa!

# PRIOR MEAN and VARIANCE of tau=1/sigma2
a/b; a/(b^2)

# M=B0^{-1}
#

######### POSTERIOR DISTRIBUTION #########
M=inv(B0)
Bn=inv(M+ t(X)%*%X)
bn= as.vector(Bn%*%( M%*%b0 +t(X)%*% y ))
## la posterior ? cos? caratterizata:
## beta | sigma2, X, dati ~ N(bn, sigma^2 * Bn) 
## sigma^2| dati ~  inv-gamma (nun/2, nun*(sigman^2)/2)

nun=(nu0+n)
T=inv(B0+ inv(t(X)%*%X))

sigma2n=(1/nun*(nu0/tau0+S2+ t(b0-betahat)%*%T%*%(b0-betahat) ))


Bn   
bn           #media a posteriori di beta
Bn*as.numeric(sigma2n) #matrice di scala a posteriori di beta
             # la covarianza A POSTERIORI tra beta0 e beta1 ? diversa da 0
Bn*as.numeric(sigma2n)*(nun/(nun-2)) # posterior covariance matrix of beta
# POSTERIOR covariance of beta0 and beta1 is NOT 0, it is about -0.35 
diag(as.numeric(sigma2hat)*inv(t(X)%*%X)) #frequentist estimate of the variance of the beta estimator 

(nun*sigma2n/2)/((nun/2)-1) # posterior mean of sigma2

# Confrontiamo matrice di covarianza (marginale) di beta a posteriori e a priori
Bn*as.numeric(sigma2n)*(nun/(nun-2)) #matrice di covarianza a posteriori di beta
(B0/tau0)*(nu0/(nu0-2)) # matrice di covarianza della marginale a priori di beta

#################### Plot della t-student bivariata #####################
### Plot the marginal posterior and prior for the vector beta
library(mvtnorm)
u1<-seq(-20, 20, by=0.1)
u2<-seq(0, 2, by=0.01)
#lower=c(-1,-1)
#fgrid <- function(xx, yy, f)
#{
#z <- matrix(nrow=length(xx), ncol=length(yy))
#for(m in 1:length(xx)){
#for(n in 1:length(yy)){
#z[m,n] <- f(c(xx[m], yy[n])) } }
#z
#}

z <- matrix(nrow=length(u1), ncol=length(u2))
for(mm in 1:length(u1)){
for(nn in 1:length(u2)){
z[mm,nn] <- dmvt( c(u1[mm], u2[nn]) , bn,Bn*(sigma2n),df=nun)
}
}

x11()

par(mfrow=c(1,2),mgp=c(1.70,.70,0))
contour(u1, u2, z, levels=c(-3.5,-5,-6.5,-8,-9.5),xlim = range(u1),
        ylim = range(u2), main="Posterior density of beta",
xlab=expression(beta0), ylab=expression(beta1))


z2 <- matrix(nrow=length(u1), ncol=length(u2))
for(mm in 1:length(u1)){
for(nn in 1:length(u2)){
z2[mm,nn] <- dmvt( c(u1[mm], u2[nn]) , b0,B0/tau0,df=nu0)
}
}

contour(u1, u2, z2, levels=c(-3.5,-5,-6.5,-8,-9.5),
        xlim = range(u1),ylim = range(u2), main="Prior density of beta",
         xlab=expression(beta0), ylab=expression(beta1))


############ Predictions 
### Plot the predictive distribution for a "new" observation corresponding to the suspect observation.
### We compare the whole predictive distribution, given x=observed covariate, to the observed value of y.
### The predictive distribution of Y_new, given x=x_{22}, is t-student (unidim in this case)
### with parameter (xnew*bn, sigma2pred, nuo+n) 
xnew=c(1,-1.45) 
mu=as.numeric(xnew %*%bn) # -6.16% expected value of the predicted Ynew
mu           # the recorded value is 58%
sigma2pred=(sigma2n*(1+ t(xnew)%*%Bn%*% xnew ))

u1=seq(-35,73,by=0.1)
z <- array(length(u1))
for(mm in 1:length(u1)){
z[mm] <- dmvt(u1[mm], mu, sigma2pred,df=nun,log=FALSE)
}

x11()
plot(u1,z,xlim=c(-35,73),type='l', xlab='Y_new', ylab=' ',
     main='Predictive distribution of Y_new at xnew=-1.45')
abline(v=58.0079,col='red')

### The observed value is in the right tail of the predictive distribution. It is clear that 
### there is a huge difference between the model (we build from the other 21 datapoints) and the 22-nd 
### datapoint. This last point is an outlier wrt the model considered so far.

###############################################################
## PRIOR g di ZELLNER: beta|sigma^2 ~ N(b0, sigma^2*B0)
##        sigma^2 ~ inv-gamma (nu0/2, nu0*(sigma0^2)/2)
##                                tau0=1/sigma0^2
## con B0=c(X^T X)^(-1), nu0=0
## cio?  beta|sigma^2 ~ N(b0, sigma^2 c(X^T X)^(-1)) 
##        \pi(sigma^2) \prop 1/ sigma2
###############################################################
b0=c(0,1) # as before
# If c=1, same weight assigned to the prior and the sample
# If c=100, the prior weight corresponds to 1% of the sample weight
c=1   #c=100
B0=c* inv(t(X)%*%X)
 

nu0=0
tau0=1/47.07 # as before

a= nu0/2
b= nu0/(2*tau0)

nu0;tau0
a;b

######### POSTERIOR DISTRIBUTION ##########
M=inv(B0)
Bn=inv(M+ t(X)%*%X)
bn= as.vector(Bn%*%( M%*%b0 +t(X)%*% y ))
## the posterior:
## beta | sigma2, X, dati ~ N(bn, sigma^2 * Bn) 
## sigma^2| dati ~ inv-gamma (nun/2, nun*(sigman^2)/2)

nun=(nu0+n)
T=inv(B0+ inv(t(X)%*%X))

sigma2n=(1/nun*(nu0/tau0+S2+ t(b0-betahat)%*%T%*%(b0-betahat) ))


Bn   
bn         #media a posteriori di beta
Bn*as.numeric(sigma2n) #matrice di scala a posteriori di beta
Bn*as.numeric(sigma2n)*(nun/(nun-2)) #matrice di covarianza a posteriori di beta
## Per la var a posteriori abbiamo ottenuto valori leggermente pi? alti del caso coniugato vero e proprio
diag(as.numeric(sigma2hat)*inv(t(X)%*%X)) #stima varianza frequentista di beta

(nun*sigma2n/2)/((nun/2)-1) #media a posteriori di sigma2

############ Plot della t-student bivariata - posterior marginal of the vector beta
library(mvtnorm)
u1<-seq(-25, 13, by=0.1)
u2<-seq(0, 2, by=0.01)

z <- matrix(nrow=length(u1), ncol=length(u2))
for(mm in 1:length(u1)){
for(nn in 1:length(u2)){
z[mm,nn] <- dmvt( c(u1[mm], u2[nn]) , bn,Bn*as.numeric(sigma2n),df=nun)
}
}

x11()
par(mfrow=c(1,1),mgp=c(1.70,.70,0))
contour(u1, u2, z, levels=c(-3.5,-5,-6.5,-8,-9.5),xlim = range(u1),
        ylim = range(u2), main="Posterior density of beta - prior di Zellner",
xlab=expression(beta0), ylab=expression(beta1))



############ Predictions corrispondenti alla prior g di Zellner
xnew=c(1,-1.45)
mu=as.numeric(xnew %*%bn) # -6.16% expected value of the predicted Ynew
mu          # the recorded value is 58%
sigma2pred=(sigma2n*(1+ t(xnew)%*%Bn%*% xnew ))

u1=seq(-35,73,by=0.1)
z <- array(length(u1))
for(mm in 1:length(u1)){
z[mm] <- dmvt(u1[mm], mu, sigma2pred,df=nun,log=FALSE)
}

x11()
plot(u1,z,xlim=c(-35,73),type='l', xlab='Y_new', ylab=' ', ylim=c(0,0.03),
     main='Predictive distribution of Y_new at xnew=-1.45 - Zellner prior')
abline(v=58.0079,col='red')

# La varianza della predittiva ? leggermente pi? grande nel caso della prior di Zellner
# per c=1

#### TRY c=100!!! What happens?

###########################################################################################
################################################
###                 JAGS use                 ###
###         Just Another Gibbs Sampler       ###
###           Alessandra Guglielmi           ###
################################################
##    http://mcmc-jags.sourceforge.net/   
##    WEB SITE: info and instructions on JAGS downloading 
################################################ 
## JAGS is Just Another Gibbs Sampler. It is a program for the analysis of Bayesian models
## using Markov Chain Monte Carlo (MCMC) which is not wholly unlike 
## OpenBUGS (http://www.openbugs.info). 
## JAGS was written with three aims in mind: to have an engine for the BUGS language that runs on Unix;
## to be extensible, allowing users to write their own functions, distributions, and samplers; 
## and to be a platform for experimentation with ideas in Bayesian modelling.
## JAGS is designed to work closely with the R language and environment for statistical
## computation and graphics (http://www.r-project.org). 
## You will find it useful to install the coda package for R to analyze the output. 
## You can also use the rjags package to work directly with JAGS from within R. 


rm(list=ls())


########## How to call JAGS from R? ##########

########## 
#Load the library
setwd("C:/alessandra/Documents/Didattica/StatisticaBayesiana/mat1617")
library(rjags)   # to interface R with JAGS
library(pscl)
data(absentee)
#help(absentee)
summary(absentee)

attach(absentee)

## create variables for regression analysis
y <- (absdem - absrep)/(absdem + absrep)*100
x <- (machdem - machrep)/(machdem + machrep)*100
x.sosp=x[22]
y.sosp=y[22]
y=y[1:21]
x=x[1:21]
#### 21 (bivariate) data points 

## JAGS takes a user's description of a Bayesian model for data, and returns an MCMC sample 
## of the posterior distribution

## Define the data (in this case, datapoints (y_i), the covariate vector (x_i), the sample size,
##                   the new covariate value xstar) in a list 
data = list(y=y[1:21],
                x=x[1:21],
                n=21,
                xstar=x.sosp)

## A list of initial value for the MCMC algorithm 
# that JAGS will implement
inits = function() {list(beta=c(0,0),
                   sigma=5,
                   ystar=0) }
## JAGS can automatically inizialize the chain, but the efficiency of the MCMC can be improved
##  if we intelligently provide reasonable starting values, i.e. values in the midst of the posterior,
## e.g. MLE, or values close to the MLE  

## The Bayesian model is in the text file "regression.bug". This file, in addition to the list "data" 
## is taken as an input by JAGS for generating the MCMC in 3 steps.
## The first step (jags.model) gets all the info into JAGS and lets JAGS figure out appropriate sampler 
## for the model.
## The second step (update) runs the chain for a burn-in period.
## The third step (coda.samples)runs and records the MCMC sample we will subsequently examine (using R).

modelRegress=jags.model("regression.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
 
### Command jags.model compiles and initializes the model described in the text file "regression.bug"
### (it is a text file, but we use *.bug to understand what file it is, but we could also name it 
### *.txt); as an input, we have "data" and JAGS then run the sampler for n.chain iterations. 
### By default n.adapt=1000.
### When a JAGS model is compiled, it may require an initial sampling phase during which the samplers
### adapt their behaviour to maximize their efficiency (e.g. a Metropolis-Hastings random walk
### algorithm may change its step size). The sequence of samples generated during this adaptive phase
### is not a Markov chain, and therefore may not be used for posterior inference on the model. 
### File "regression.bug" specify both the likelihodd and the prior 
### beta has 2 components which are iid Gaussian, with mean =0 and variance=10thousand
### sigma is Uniform on the interval (0,100)
### In this case, the PRIOR IS NOT CONJUGATE.
### From rjags manual:
# jags.model is used to create an object representing a Bayesian graphical model, 
# specified with a BUGS-language description of the prior distribution, and a set of data.
### USAGE:
# jags.model(file, data=sys.frame(sys.parent()), inits, n.chains = 1, n.adapt=1000, quiet=FALSE)
# Se inits ? omessa, i valori iniziali saranno specificati automaticamente. ATTENZIONE:
# talvolta la catena si blocca subito se inits NON ? BEN specificato (anche nel caso in cui venga calcolato automaticamente)
# Durante le prime n.adapt iterazioni, la catena ? adattativa, per es. potrebbe NON essere markoviana,
# o la proposal distribution potrebbe cambiare

#####save(modelRegress,file='modelregress.Rdata') 

update(modelRegress,n.iter=19000) #this is the burn-in period
## this function DOES NOT record the sampled values of parameters

## THIRD STEP: we tell JAGS to generate MCMC samples that we will use to represent the posterior distribution 
variable.names=c("beta", "sigma", "ystar") #parameters - see the file .bug - that will have 
                                           # their values recorded
n.iter=50000 
thin=10  
## the final sample size, i.e. the number of iterations to build the ergodic average, will be 5K

outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)
# lancio la catena per altre n.iter iterazioni, specificando l'eventuale thinning
# salvo l'intera catena dei parametri monitorati; 
# the OUTPUT is mcmc.list object - coda-formatted object; 
# it needs to be converted into a matrix, in order to be "readable".
# 
##################################################################################
# SUMMING UP, to produce the output (the MCMC chain) using rjags 4 steps are required: 
#
# 1. Define the model using the BUGS language in a separate file.
# 2. Read in the model file using the jags.model function. This creates an object of 
#    class ?jags?.
# 3. Update the model using the update method for ?jags? objects. This constitutes a
     ?burn-in? period.
# 4. Extract samples from the model object using the coda.samples function. This creates 
#    an object of class ?mcmc.list? which can be used to summarize the posterior 
#    distribution. The coda package also provides convergence diagnostics to check that 
#    the output is valid for analysis (see Plummer et al 2006).

#### save(outputRegress,file='Jackman_regr_output.Rdata')

#### ANALISI DELL'OUTPUT
library(coda)       # per leggere catena mcmc
library(plotrix)    # per fare plot CIs

##### CARICA le traiettorie della catena
#### load('Jackman_regr_output.Rdata')


data.out=as.matrix(outputRegress) # trasform the mcmc.list into a matrix,
data.out=data.frame(data.out)     # or, better, into a dataframe (easiest to handle in R)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain # this is the final sample size


#############################################
##  Summary statistics for the posterior
#############################################
summary(data.out)
head(data.out)

## Use either packege CODA, or standard tools in R
par(mfrow=c(1,3))
acf(data.out[,'beta.1.'],lwd=3,col="red3",main="autocorrelation of beta1")
acf(data.out[,'beta.2.'],lwd=3,col="red3",main="autocorrelation of beta2")
acf(data.out[,'sigma'],lwd=3,col="red3",main="autocorrelation of sigma_res")

#par(mfrow=c(1,1))

###Some nice graph to presents the posterior samples

#sub-chain containing the beta sample
beta.post <- data.out[,1:2]
#posterior mean of the beta parameters
beta.bayes  <- apply(beta.post,2,"mean")
beta.bayes
#sub-chain whit the sigma_res samle
sig.post= data.out[,'sigma']

#### Stime ai minimi quadrati
pippo= lm(y ~ x)
summary(pippo)

pippo$coefficients
beta.bayes

## Representation of the posterior chain of  beta0
chain <- beta.post[,1]
#Divide the plot device in three sub-graph regions
#two square on the upper and a rectangle on the bottom
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot of beta0")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation of beta0")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior of beta0",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2)
## Posterior credible interval of beta0
quantile(chain,prob=c(0.025,0.5,0.975))

## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))


x11()
##Representation of the posterior chain of beta1
chain <- beta.post[,2]
#Divide the plot device 
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot of beta1")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation of beta1")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior of beta1",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2)
#### Posterior credible interval of  beta1
quantile(chain,prob=c(0.025,0.5,0.975))

## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))



x11()
##Representation of the posterior chain of sigma_res
chain <-sig.post
#Divide the plot device 
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot of sigma_res")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation of sigma_res")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior of sigma_res",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2)
## Posterior Credible bound of sigma_res
quantile(chain,prob=c(0.05,0.5,0.95))

## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother"),lwd=c(2,2,2), col=c("red","red","blue"),lty=c(1,2,1))



#close all the graphic windows
graphics.off()

###############################################
### Predictive distibution of Y_new at x=xstar=suspected value  
### "ystar" column in the output contains simulated values from N(beta0+beta1*xstar,1/tau)
### for each value (beta0,beta1,tau) simulated a posteriori, i.e.
### "ystar" is a MCMC sample from the predictive distribution of Ynew corresponding at xstar
ystar.pred=data.out[,'ystar'] 


## Representation of the chain y.star
chain <- ystar.pred
#Divide the plot device in three sub-graph regions
#two square on the upper and a rectangle on the bottom
layout(matrix(c(1,2,3,3),2,2,byrow=T))
#trace-plot of the posterior chain
plot(chain,type="l",main="Trace plot of y.star")
# autocorrelation plot
acf(chain,lwd=3,col="red3",main="autocorrelation of ystar")
#Histogram
hist(chain,nclass="fd",freq=F,main="Posterior of ystar",col="gray") 
## Overlap the kernel-density 
lines(density(chain),col="blue",lwd=2,xlim=c(0,60))
## Posterior credible interval of ystar
quantile(chain,prob=c(0.025,0.5,0.975))

## Display the posterior credible interval on the graph
abline(v=quantile(chain,prob=c(0.025)),col="red",lty=2,lwd=2)
abline(v=quantile(chain,prob=c(0.5)),col="red",lty=1,lwd=2)
abline(v=quantile(chain,prob=c(0.975)),col="red",lty=2,lwd=2)
## Add the legend to the plot
#Aggiungiamo una retta relativa all'osservazione di 
#y[22]=y.sops
abline(v=y.sosp,col="magenta",lwd=2)

legend("topright",legend=c("posterior median", "95% Credible bounds","kernel density smoother","observed value at xstar"),lwd=c(2,2,2,2), col=c("red","red","blue","magenta"),lty=c(1,2,1,1))

### The observed value at xstar is OUTSIDE the 95% CI of the predictive distribution, 
### i.e. it is rather extreme w.r.t. the predictive, i.e. it is in the tails of the predictive distribution...
### IT IS pretty suspect!!
 

################################################
####  Credible band for the regression line  ###
################################################
## For any x in a grid of point, we compute 90% CI of Y_new(x)| x
### Valori sui quali andremo a valutare la retta di regressione bayesiana, 
### cio? la retta y = E[ beta[1]|dati]+E[ beta[2]|dati]*x
 
x.gr=seq(-50,100,length=200)

#Pred is a matrix: each column contains MCMC values from the predictive distribution of Y_new(x)| x
# each row is a draw from the predictive of the regression line (as a function of x)
Pred <- matrix(ncol=200,nrow=n.chain)
for(i in 1:n.chain){
	#print(i)
	Pred[i,]=beta.post[i,1]+beta.post[i,2]*x.gr+rnorm(1,mean=0,sd=sig.post[i]) #stesso errore per ogni valore di x.gr
}

x11()
layout(matrix(c(1,2,3,3),1,1,byrow=T))
### Plot the mean (over columns) of  all these values 
predit=apply(Pred,2,"mean")
plot(x.gr,predit,type="l",col="red",lwd=2)#this is a line: E(beta0|data)+E(beta_1|data)*x
points(x,y,col="green",pch="*",cex=2) #datapoints in green!
gra=gray(1:100/100)
gra=rep(gra,10)
### Questa che ho disegnato con i comandi qui sopra ? una retta, perch? ? la retta che associa 
### ad ogni valore di x.gr il valore E(Y^*(x.gr)|dati)= E[E(Y^*(x.gr)|par,dati)|dati]
###                                                   = E[ beta[1]+beta[2]*x.gr |dati]
###                                                   = E[ beta[1]|dati]+E[ beta[2]|dati]*x.gr

#### Disegno le prime 600 rette di regressione simulate 
for(i in 1:600){lines(x.gr,Pred[i,],col=gra[i])}  ### These are lines, since we used the same simulated value
                                                  ### rnorm(1,mean=0,sd=sig.post[i]) for all x.gr
predit.qua=apply(Pred,2,"quantile",prob=c(0.05,0.95)) #quantiles of each predictive distribution of  Y_new(x)| x       
### For any value of x.gr in the grid, we compute 3 quantile of the (simulated) distribution of  
###       y=beta0+beta1*x.gr+err 
### and then we plot it as a function of x.gr
lines(x.gr,predit.qua[1,],col="red",lwd=3,lty=2) ### THIS IS NOT necessarily a line! 
lines(x.gr,predit.qua[2,],col="red",lwd=3,lty=2) ### THIS IS NOT necessarily a line! 
lines(x.gr,predit,type="l",col="red",lwd=3)  ### this IS a line, i.e. it is the Bayesian regression line
                                             ### computed before
### Invece le curve dei quantili sono curve (NON SONO NECESSARIAMENTE RETTE) che rappresentano le bande di credibilit? 
### della retta di regressione  y = E[ beta[1]|dati]+E[ beta[2]|dati]*x
### Se "affetto" questo plot con una retta x= x0, ottengo la distribuzione predittiva
### corrispondente a x0 

points(x.sosp,y.sosp,pch=16,col="magenta") ## disputed data point
text(x.sosp,y.sosp,"Disputed\nElection",cex=.75,adj=1.25,col="magenta")

#######################################################################################


