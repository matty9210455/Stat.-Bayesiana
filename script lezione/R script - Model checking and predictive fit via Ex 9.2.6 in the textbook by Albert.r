#### LINEAR REGRESSION with iid N(0,sigma^2) errors
################################
# Section 9.2.6 An Example
################################
library(LearnBayes) 

data(birdextinct)
attach(birdextinct)

birdextinct
# The DATASET concerns EXTINCTION of BIRDS; 
# Measurements on breedings ("in cova") pairs of landbird species were collected from 16 islands about Britain 
# over several decades. For each pairs (62 in total), the dataset records:
# SPECIES = name of bird species
# TIME = average time of extinction on the islands
# NESTING = average number of nesting pairs
# SIZE = size of the species, 1 or 0 if large or small
# STATUS = status of the species, 1 or 0 if resident or migrant.
## LINEAR MODEL for TIME, using NESTING, SIZE, STATUS as covariates.

plot(birdextinct)

# Conviene passare al logaritmo di TIME
logtime=log(time)
plot(nesting,logtime)
# POSITIVE CORRELATION between these variables
 
# OUTLIER: probably datapoints with log(time)>3
out = (logtime > 3)
#aggiungo le specie relativa agli OUTLIER sull'ultimo grafico che ho generato
text(nesting[out], logtime[out], label=species[out], pos = 2)	

x11()
# Command JITTER adds a small amount of noise to a numeric vector
par(mfrow=c(1,2))
plot(jitter(size),logtime,xaxp=c(0,1,1))
plot(jitter(status),logtime,xaxp=c(0,1,1))

##### ML estimates#####
# "lm" returns the model matrix and the response, as well
fit=lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
summary(fit)
# NESTING is a strong effect 
# estimate of the effect of SIZE is negative: animali piccoli (SIZE=1) hanno tempi di estinzione
# pi? brevi

fit$x
fit$y
plot(fit$residual)

##### CAMPIONAMENTO della POSTERIOR con REFERENCE PRIOR (impropria)
##### POTETE GUARDARLO da SOLI - INSIEME VEDIAMO l'analisi con la prior g di Zellner
# La funzione di LearnBayes si chiama BLINREG = Bayesian Linear Regression 
help(blinreg)
# Mi restituisce un campione Monte Carlo - dunque iid - dalla distribuzione finale di 
# theta=(beta,sigma) con prior(sigma^2)=1/sigma^2
# INPUT: vettore risposta y, matrice disegno X, numerosit? del campione in uscita
theta.sample=blinreg(fit$y,fit$x,2)
theta.sample

#Prima campiono tau da una gamma((n-p)/2,(s^2)/2) con s^2=sum(fit$residual^2),
#poi campiono beta da una gaussiana p-variata di media MLE(beta) e
#matrice di varianza e covarianza ((X^t X) tau)^(-1) 
theta.sample=blinreg(fit$y,fit$x,5000)

par(mfrow=c(2,2))
hist(theta.sample$beta[,1],main="INTERCEPT",xlab=expression(beta[0]),prob=T)
hist(theta.sample$beta[,2],main="NESTING",xlab=expression(beta[1]),prob=T)
hist(theta.sample$beta[,3],main="SIZE",xlab=expression(beta[2]),prob=T)
hist(theta.sample$beta[,4],main="STATUS",xlab=expression(beta[3]),prob=T)
#check rapido che covareiate sono diverse da 0
X11()
par(mfrow=c(1,2))
#La distr. FINALE di tau=1/sigma^2 ? gamma((n-p)/2, ((n-p)/2)*S^2)
#Ecco l'istogramma:
hist(1/((theta.sample$sigma)^2),main="PRECISIONE TAU",xlab=expression(tau),prob=T)
hist(theta.sample$sigma,main="ERROR SD",xlab=expression(sigma),prob=T)


#Le medie e le mediane a posteriori dei beta (che coincidono perch?
#le distrib finali sono simmetriche) sono uguali alle stime MLE - a meno
#di errori di simulazione - perch? abbiamo usato la prior di Jeffreys
apply(theta.sample$beta,2,quantile,c(.05,.5,.95))
apply(theta.sample$beta,2,mean)
fit$coefficients

quantile(theta.sample$sigma,c(.05,.5,.95))
#stima bayesiana di sigma^2 - media a posteriori - ? maggiore della stima ai minimi quadrati
mean(theta.sample$sigma)
#stima frequentista ai mimini quadrati di sigma
sqrt(sum(fit$residuals^2)/(62-4))


###### Estimating mean extinction times
# Vogliamo stimare E(y_new|beta, x1)=x_1*beta, dove x_1 ? una "nuova" matrice di covariate;
# si noti che questa ? semplicemente  una combinazione lineare di beta.
cov1=c(1,4,0,0) #A
cov2=c(1,4,1,0) #B
cov3=c(1,4,0,1) #C
cov4=c(1,4,1,1) #D
X1=rbind(cov1,cov2,cov3,cov4)
# La funzione blinregexpected prende in INPUT la nuova matrice disegno 
# e il campione generato dalla posterior e ne fa la combinazione lineare
mean.draws=blinregexpected(X1,theta.sample)

c.labels=c("A","B","C","D")
par(mfrow=c(2,2))
for (j in 1:4)
hist(mean.draws[,j], main=paste("Covariate set",c.labels[j]),xlab="log TIME",prob=T)
 
######## Calcolo delle PREDITIVE - Predicting extinction time
# Prima abbiamo simulato la distribuzione finale di E(y_new|X1) - che ? una funzione dei 4 prametri di regressione.
# Ora invece vogliamo la distribuzione predittiva di y_new, cio?
# la legge condizionale di y_new, dato il vettore delle osservazioni y.
cov1=c(1,4,0,0)
cov2=c(1,4,1,0)
cov3=c(1,4,0,1)
cov4=c(1,4,1,1)
X1=rbind(cov1,cov2,cov3,cov4)
#La funzione blinregpred genera un valore dalle verosimiglianze valutate
#in ognuno dei valori di theta generati dalla posterior, cio? generata da una 
#gaussiana di dim m=4 con media X_1 beta^(j) e matrice di var =(X^tX *tau^(j))^(-1)
pred.draws=blinregpred(X1,theta.sample)

c.labels=c("A","B","C","D")
par(mfrow=c(2,2))
for (j in 1:4)
hist(pred.draws[,j],main=paste("Covariate set",c.labels[j]),xlab="log TIME",prob=T)


###########################################################################
############################ Zellner's g-prior ############################
help(blinreg)
# The output is a Monte Carlo sample (iid) from the posterior of 
# theta=(beta,sigma)
# Gives a simulated sample from the joint posterior distribution of the regression vector 
# and the error standard deviation for a linear regression model with a (noninformative or) Zellner g prior.
# INPUT: response vector y, design matrix X, num. of iterations, prior list with components c0 and beta0 of Zellner's g prior

# Zellner g prior N_p(b0, sigma^2 (c0 (Xt X)^{-1})) * (1/sigma^2)
prior=list(b0=c(0,0,0,0), c0=10)  #c0\in [10,100]
theta2.sample=blinreg(fit$y,fit$x,5000,prior=prior)

# La distribuzione finale si fattorizza come al solito;
# la distribuzione di beta, dato tau e y, ? 
# N_p(c0/(c0+1)(b0/c0 + hatbeta), sigma^2 * c0/(c0+1)* (Xt X)^{-1})) ),
# mentre la distribuzione di tau, dato y, ? gamma(n/2, s^2/2+ ...)

x11()
par(mfrow=c(2,2))
hist(theta2.sample$beta[,1],main="INTERCEPT",xlab=expression(beta[0]),prob=T)
hist(theta2.sample$beta[,2],main="NESTING",xlab=expression(beta[1]),prob=T)
hist(theta2.sample$beta[,3],main="SIZE",xlab=expression(beta[2]),prob=T)
hist(theta2.sample$beta[,4],main="STATUS",xlab=expression(beta[3]),prob=T)

X11()
par(mfrow=c(1,2))
hist(1/((theta2.sample$sigma)^2),main="PRECISIONE TAU",xlab=expression(tau),prob=T)
hist(theta2.sample$sigma,main="ERROR SD",xlab=expression(sigma),prob=T)


#Le medie e le mediane a posteriori dei beta 
apply(theta2.sample$beta,2,quantile,c(.05,.5,.95))
apply(theta2.sample$beta,2,mean)
fit$coefficients # stime frequentiste 

quantile(theta2.sample$sigma,c(.05,.5,.95))
#stima bayesiana di sigma^2 - media a posteriori - ? maggiore della stima ai minimi quadrati
mean(theta2.sample$sigma)
#stima frequentista ai mimini quadrati di sigma
sqrt(sum(fit$residuals^2)/(62-4))


###### Estimating mean extinction times
# Vogliamo stimare E(y_new|beta, x1)=x_1*beta, dove x_1 ? una "nuova" matrice di covariate;
# si noti che questa ? semplicemente  una combinazione lineare di beta.
cov1=c(1,4,0,0) #A
cov2=c(1,4,1,0) #B
cov3=c(1,4,0,1) #C
cov4=c(1,4,1,1) #D
X1=rbind(cov1,cov2,cov3,cov4)
# La funzione blinregexpected prende in INPUT la nuova matrice disegno 
# e il campione generato dalla posterior e ne fa la combinazione lineare
mean.draws=blinregexpected(X1,theta2.sample)

c.labels=c("A","B","C","D")
par(mfrow=c(2,2))
for (j in 1:4)
hist(mean.draws[,j], main=paste("Covariate set",c.labels[j]),xlab="log TIME",prob=T)

 
######## Calcolo delle PREDITIVE  con la PRIOR g di ZELLNER - Predicting extinction time
# Prima abbiamo simulato la distribuzione finale di E(y_new|X1) - che ? una funzione dei 4 prametri di regressione.
# Ora invece vogliamo la distribuzione predittiva di y_new, cio?
# la legge condizionale di y_new, dato il vettore delle osservazioni y.
cov1=c(1,4,0,0)
cov2=c(1,4,1,0)
cov3=c(1,4,0,1)
cov4=c(1,4,1,1)
X1=rbind(cov1,cov2,cov3,cov4)
#La funzione blinregpred genera un valore dalle verosimiglianze valutate
#in ognuno dei valori di theta generati dalla posterior, cio? generata da una 
#gaussiana di dim m=4 con media X_1 beta^(j) e matrice di var =(X^tX *tau^(j))^(-1)
pred.draws=blinregpred(X1,theta2.sample)

c.labels=c("A","B","C","D")
par(mfrow=c(2,2))
for (j in 1:4)
hist(pred.draws[,j],main=paste("Covariate set",c.labels[j]),xlab="log TIME",prob=T)

#####################################################################################


############## MODEL CHECKING via POSTERIOR PREDICTIVE DISTRIBUTIONS ################
# From the book Gelman et al. (2013) - $Y_i^{new}$ is "the replicated data 
# that could have been observed, or, to think predictively, as the data that
# we would see tomorrow if the experiment that produced $y_i$ today were 
# replicated with the same model and the same value of $\theta$ that produced
# the observed data"
# If the model fits, then replicated data generated under the model should look similar to observed data.
# To put another way, the observed data should look PLAUSIBLE under the posterior predictive distribution. 
# This is a self-consistency test! 
# Any systematic differences between the simulations and the data indicate 
#      potential failings of the model. 
#### As an ALTERNATIVE, apply the CROSS-VALIDATION approach - see BELOW!


help(blinregpred)
# Simulo la distribuzione predittiva di Y_1,...,Y_62, ottenendo un campione
# iid da questa, lungo quanto era il campione dalla finale, e ne calcolo i 
# quantili di ordine 0.05 e 0.95
pred.draws=blinregpred(fit$x,theta2.sample)
pred.sum=apply(pred.draws,2,quantile,c(.05,.95))
# Plot the 90% CI and the datapoints
par(mfrow=c(1,1))
ind=1:length(logtime)
matplot(rbind(ind,ind),pred.sum,type="l",lty=1,col=1,xlab="INDEX",ylab="log TIME")
points(ind,logtime,pch=19)

# Consider species datapoints outsice these interavals as OUTLIERS
out=(logtime>pred.sum[2,]| logtime<pred.sum[1,])
text(ind[out], logtime[out], label=species[out], pos = 4)

################### BAYESIAN PREDICTIVE p-values ####################
# min( P(Y_i^{new}>y_i|x_i,data), P(Y_i^{new}<y_i|x_i, data)) 
# Compute how many simulated values from the predictive are > y_i (over the total num of simulations):
# this is an estimate of P(Y_i^{new}>y_i|data, x_i)
ind_pvalue=t((t(pred.draws)>logtime))
pred.suptail=apply(ind_pvalue,2,sum)/5000
prob.tail<-rep(0, 62) 
for (i in 1 : 62)
  prob.tail[i]=min(pred.suptail[i],1-pred.suptail[i])
## this is the BAYESIAN PREDICTIVE p-value of item i

# if this number is close to 1/2, the corresponding predictive "explain" well the datapoint
# if this number is close to 0, the datapoint is an outlier for the model
x11()
par(mfrow=c(1,2))
plot(ind, prob.tail)
abline(h=0.05)

plot(ind,logtime,pch=19)
# OUTLIERS: those species such that the BAYESIAN PREDICTIVE p-value is smaller than 0.05 (or than 0.1)
out2=(prob.tail<0.05)
text(ind[out2], logtime[out2], label=species[out2], pos = 4)




############ Model checking via bayes residuals ##############
# In Chaloner & Brant (1988) si definisce OUTLIER (in un modello di regressione)
# come quella osservazione con "errore osservato sorprendente", oppure
# pi? in generale, ogni osservazione che NON ? stata prodotta dal meccanismo che ha 
# generato la maggior parte dei dati"
# A priori la prob. che una osservazione sia un outlier ? 
# P(|eps_i|/sigma>k)=2Phi(-k); bisogner? calcolare A POSTERIORI la prob di questo 
# evento, e confrontarla con quella a priori
# NON ? un metodo molto CONVINCENTE...

#### Richiamo i comandi 
#library(LearnBayes)
#data(birdextinct)
#attach(birdextinct)
#logtime=log(time)
#fit=lm(logtime~nesting+size+status,data=birdextinct,x=TRUE,y=TRUE)
#theta.sample=blinreg(fit$y,fit$x,5000)
####

#help(bayesresiduals)
##Calcoliamo la probabilit? predittiva che una osservazione sia classificata
##come OUTLIER con k=2 - formula di Chaloner & Brant (1988)
#
#prob.out=bayesresiduals(fit,theta.sample,2)
#par(mfrow=c(1,1))
#plot(nesting,prob.out,cex=1)
##Ritengo OUTLIER quelle specie per cui queste probabilit? sono maggiori di 
##0.35, per esempio
##out = (prob.out > 0.35)
#
#text(nesting[out], prob.out[out], label=species[out], pos = 4)

##############################################################
############# Predictive Bayesian residuals ##################
## Compute mean and sd of all n predictive distributions first. 
## Then compute the difference between the predictive mean and observed datapoint, over 
## the predictive sd. This is the "predictive Bayesian residual".
## Consider as an OUTLIER all data such that the corresponding residual is LARGE 
## e.g. its absolute value is larger than 3 
# Simulo la distribuzione predittiva di Y_1,...,Y_62; la prior ? quella di Zellner, quindi 
# ? una caso particolare della coniugata, e quindi riesco a simulare un campione iid dalla posterior,
# usando la funzione blinreg(fit$y,fit$x,5000,prior=prior) - FATTO DI SOPRA!

pred.mean=apply(pred.draws,2,mean)
pred.sd=apply(pred.draws,2,sd)
 
bres= (pred.mean- logtime)/pred.sd   # residui bayesiani
out2 = (abs(bres) > 2) #as a reference value we take 2
# or 1.8

par(mfrow=c(1,1))
plot(ind,bres,cex=1)
text(nesting[out2], bres[out2], label=species[out2], pos = 4)

# Predictive goodness-of-fit: SUM of the predictive Bayesian residuals
# to compare different models 
# The "best" model is the one with the smallest value for this index

sum(bres^2)

###################################################################
#######################  LPML and  CPO_i  #########################
#################### a CROSS-VALIDATION approach ##################
###################################################################
# Conditional density of each Y_i, given parameters, is
# Gaussian with mean X*beta and variance sigma^2
mean.draws.data=blinregexpected(fit$x,theta2.sample)
dim(mean.draws.data)

cpo=seq(1,62)
for (i in 1:62){
cpo[i]=1/mean(1/dnorm(fit$y[i],mean.draws.data[,i], theta2.sample$sigma))
}

1/mean(1/dnorm(fit$y[2],mean.draws.data[,2], theta2.sample$sigma)) # this is  CPO[2] of item i=2

LPML=sum(log(cpo))
LPML


### Repeat the analysis WITHOUT covariate "status"
##### Stima ai MINIMI QUADRATI (stime frequentiste) #####
fit2=lm(logtime~nesting+size,data=birdextinct,x=TRUE,y=TRUE)
summary(fit2)

fit2$x
fit2$y


##### SAMPLING from the psoterior with Zellner g prior
# Mi restituisce un campione Monte Carlo - dunque iid - dalla distribuzione finale di 
# theta=(beta,sigma) con prior(sigma^2)=1/sigma^2
# INPUT: vettore risposta y, matrice disegno X, numerosit? del campione in uscita
prior3=list(b0=c(0,0,0), c0=10)  #c0\in [10,100]

theta3.sample=blinreg(fit2$y,fit2$x,5000, prior=prior3)

##### LPML
mean.draws3.data=blinregexpected(fit2$x,theta3.sample)
dim(mean.draws3.data)

cpo2=seq(1,62)
for (i in 1:62){
cpo2[i]=1/mean(1/dnorm(fit$y[i],mean.draws3.data[,i], theta3.sample$sigma))
}


LPML2=sum(log(cpo2))
LPML2
### The best model is the one with LARGEST LPML: which model is better, then?

LPML; LPML2
