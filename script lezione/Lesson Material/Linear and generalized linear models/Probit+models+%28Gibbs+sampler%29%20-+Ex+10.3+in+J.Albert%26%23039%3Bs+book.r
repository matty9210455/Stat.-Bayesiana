# Sampling from a TRUNCATED distribution, i.e
# conditioning to the r.v. to belong to an interval (strinctly contained in the support).
# Gaussian r.v. X with, parameters mu=1, sigma=1
# TRUNCATED in A=(a1,a2)=(0,+\infinity), 
# that is X, but conditioning X to assume values in A

rm(list=ls())
set.seed(3)
mu=1
sigma=1

Fa1 <- pnorm(0,mu,sigma)
Fa2 <- 1
u1 <- runif(2000)
u1 <- u1*(Fa2-Fa1)+Fa1
#u1=runif(2000,min=Fa1,max=Fa2)
#par(mfrow=c(1,2))
x=qnorm(p=u1,mu,sigma)# la funzione qnorm calcola phi^{-1}, l'inversa della f.r. di una gaussiana standard

hist(x,freq=F,col="red")
curve( dnorm(x,mu,sigma)/(1-pnorm(0,mu,sigma)),from=0,to=6,col="blue",add=T,lwd=3)

##########################################################################
## Example 10.3 dal libro di J. Albert "Bayesian computations with R".
## GLMs for binary data: a regression model with PROBIT link 
########################################################################## 
#
# Data from DONNER party, a group of 'pioneers' who crossed Sierra Nevada, 
# in a wagon train (in carovana)in 1846-47.  Most pioneers starved to death.
# Binary response representing death/survival. 
# We analyze this dataset via a GLM with PROBIT link. 
# n=45
#
# Analisi bayesiana col package LearnBayes
library(LearnBayes)

data(donner)
attach(donner)
names(donner)
age
male
survival
 
# Response = survival: 0-> morto, 1-> vivo
# p_i = survival probability 
# COVARIATE: age (min age: 15 years)
#            male: 0-> donna, 1-> uomo

#matrice disegno con la prima colonna di '1' per includere l'intercetta
X=cbind(1,age,male)

# La funzione glm calcola la MLE per i parametri di regressione;
# è necessario specificare la distribuzione nella famiglia esponenziale e la link function
# Si deve togliere 'intercetta', perchè è già stata aggiunta alla matrice disegno
fit=glm(survival~X-1,family=binomial(link=probit))
# Equivalently: glm(survival ~ male + age,family=binomial(link = probit)) 
summary(fit)

################# CALCOLO della POSTERIOR con PRIOR IMPROPRIA ##############
# La funzione bayes.probit implementa un Gibbs sampler per il calcolo della posterior 
# del parametro (beta, Z), beta= tre parametri di regressione, Z= vettore delle n latenti
# PRIOR: qui abbiamo usato una prior improria per beta, ma è anche possibile assegnare
# una prior gaussiana con un certo vettore di medie e una certa matrice di precisione - vedere HELP
# PIU' SOTTO: PRIOR PROPRIA "semiconiugata"

help(bayes.probit)
print(bayes.probit)
## Il punto cruciale è la simulazione delle full-conditionals di Z_i| beta, y_i
## Nella funzione, questo è il passo z = rtruncated(N, LO, HI, pnorm, qnorm, X %*% beta, 1)
## IN GENERALE, per simulare una v.a. Z da una N(mu,1) troncata sull'intervallo (a,b), 
## si scrive:
## Z=\Phi^{-1}(\phi(a-mu)+ U*(\Phi(b-mu)-\Phi(a-mu)) ) + mu, con U~ U(0,1)
## a=0, b=+\infty
u1 <- runif(2000)
mu=1
Fa1 <- pnorm(0,mu) #\Phi(a-mu)
Fa2 <- 1  #\Phi(b-mu)
u1 <- runif(2000)
u1 <- u1*(Fa2-Fa1)+Fa1
x=qnorm(p=u1)+mu
hist(x,freq=F,col="red")

## L'altro passo è la full conditional dei beta, che equivale alla posterior di un modello lineare gaussiano 
## quando i dati sono le variabili latenti Z_i

# La funzione vuole in input il vettore delle risposte binarie y, la matrice disegno X, e 
# il numero di iterazioni del Gibbs sampler m
m=10000
fitBayes=bayes.probit(survival,X,m) # PRIOR IMPROPRIA
?bayes.probit
# Gives a simulated sample from the joint posterior distribution of the regression vector for 
# a binary response regression model with a probit link and a informative normal(beta, P) prior. 
# bayes.probit(y,X,m,prior=list(beta=0,P=0))
## Se NON specifico la lista "prior", assumo la prior IMPROPRIA per i parametri di regressione beta,
## cioè proporzionale ad una costante... 

# In uscita dà beta, la matrice dei valori simulati dei parametri di regressione
# e log.marg, una stima Monte Carlo della "log marginal likelihood of the model" 

# Medie e varianze a posteriori dei coefficienti di regressione stimate con metodo MCMC
apply(fitBayes$beta,2,mean)
apply(fitBayes$beta,2,sd)
# Sono simili alle stime frequentiste (vedi full-conditionals dei beta)
fit$coefficients

par(mfrow=c(1,3))
plot(density(fitBayes$beta[,1]), xlab='beta0 - intercept',main=' ') 
plot(density(fitBayes$beta[,2]), xlab='beta1 - age',main=' ') 
plot(density(fitBayes$beta[,3]), xlab='beta2 - male',main=' ') 
# Sia beta1 che beta2 sono significativi, e in entrambi i casi la posterior concentra 
#   massa sui valori negativi.
# La probabilità di sopravvivere decresce al crescere dell'età ed è minore per gli uomini
#   invece che per le donne.
# La probabilità di sopravvivenza, fissate AGE e MALE, è \Phi(beta0+beta1*age+beta2*male),
# e quindi è funzione dei parametri di regressione; vediamo la distribuzione a posteriori di 
# questa quantità, per un uomo (male=1), al variare di age.
# Coincide con la banda di previsione della sopravvivenza predittiva
#     di un "nuovo" paziente maschio.
a=seq(15,65)
X1=cbind(1,a,1)
#C'è una funzione del package che fa questo:
p.male=bprobit.probs(X1,fitBayes$beta)
# Disegnamo il grafico della mediana a posterior di queste quantità (al variare dell'età),
# insieme ai quantili di ordine 0.05 e 0.95
par(mfrow=c(1,1))
plot(a,apply(p.male,2,quantile,.5),type="l",ylim=c(0,1), xlab="age",ylab="Probability of Survival")
lines(a,apply(p.male,2,quantile,.05),lty=2)
lines(a,apply(p.male,2,quantile,.95),lty=2)


###################################################
# Section 10.3.2  Proper priors and model selection
###################################################

library(LearnBayes)
data(donner)
y=donner$survival
X=cbind(1,donner$age,donner$male)

# PRIOR of beta: N_p with mean=beta00 and PRECISION matrix P0
# Fisso la media pari al vettore nullo e la matrice di varianza c*(X^T X)^(-1)
#  (è simile alla prior g di Zellner per il caso lineare classico)
beta00=c(0,0,0); c0=100 #la prior ha un peso dell'1% rispetto al campione
## Provate con c0=1 e c0=10
P0=t(X)%*%X/c0
#prior= list with components beta, the prior mean, and P, the prior PRECISION matrix
#         P è quella che abbiamo indicato con B0^(-1)
inv=function(X)
{
# RETURN THE INVERSE OF THE SYMMETRIC MATRIX X
EV=eigen(X)
EV$vector%*%diag(1/EV$values)%*%t(EV$vector)
}
B0=inv(P0)
B0    #B0= c0*(t(X)%*%X)^{-1}


# bayes.probit function implements a Gibbs sampler to compute the posterior  of  
# (beta, Z) parameter, beta= 3 regression parameters, Z= vector of  n latent variables
# Gaussian PRIOR: full conditionals have closed-form analytic expressions

help(bayes.probit)
print(bayes.probit)
## Il punto cruciale è la simulazione delle full-conditionals di Z_i| beta, y_i
## Nella funzione, questo è il passo z = rtruncated(N, LO, HI, pnorm, qnorm, X %*% beta, 1)
## IN GENERALE, per simulare una v.a. Z da una N(mu,1) troncata sull'intervallo (a,b), 
## si scrive:
## Z=\Phi^{-1}(\phi(a-mu)+ U*(\Phi(b-mu)-\Phi(a-mu)) ) + mu, con U~ U(0,1)
## a=0, b=+\infty
u1 <- runif(2000)
mu=1
Fa1 <- pnorm(0,mu) #\Phi(a-mu)
Fa2 <- 1  #\Phi(b-mu)
u1 <- runif(2000)
u1 <- u1*(Fa2-Fa1)+Fa1
x=qnorm(p=u1)+mu
hist(x,freq=F,col="red")

## L'altro passo è la full conditional dei beta, che equivale alla posterior di un modello lineare gaussiano 
## quando i dati sono le variabili latenti Z_i

# La funzione vuole in input il vettore delle risposte binarie y, la matrice disegno X, e 
# il numero di iterazioni del Gibbs sampler m e i parametri della PRIOR gaussiana per beta 
# Gives a simulated sample from the joint posterior distribution of the regression vector for 
# a binary response regression model with a probit link and a informative normal(beta, P) prior. 
# bayes.probit(y,X,m,prior=list(beta=0,P=0))
## Se NON specifico la lista "prior", assumo la prior IMPROPRIA per i parametri di regressione beta,
## cioè proporzionale ad una costante... 

# In uscita dà beta, la matrice dei valori simulati dei parametri di regressione
# e log.marg, una stima Monte Carlo della "log marginal likelihood of the model" 



m=10000
fitBayes=bayes.probit(survival,X,m,list(beta=beta00,P=P0)) # PRIOR PROPRIA

# Posterior mean and variances of regression coefficients via the MCMC
apply(fitBayes$beta,2,mean)
apply(fitBayes$beta,2,sd)
# Le medie a posteriori sono abbastanza simili  alle stime frequentiste, perchè costruite a partire 
# dal Gibbs Sampler in cui campiono i beta da una Gaussiana, con media data dalla media ponderata 
# della stima di max ver (che pesa 100 volte di più della media prior) e della media a priori 
fit$coefficients

par(mfrow=c(1,3))
plot(density(fitBayes$beta[,1]), xlab='beta0 - intercept',main=' ') 
plot(density(fitBayes$beta[,2]), xlab='beta1 - age',main=' ') 
plot(density(fitBayes$beta[,3]), xlab='beta2 - male',main=' ') 
mean(fitBayes$beta[,2]<0)
mean(fitBayes$beta[,3]<0)
# Both beta1 and beta2 are significative, and in both cases the marginal posterior is concentrated  
#   on negative values.
# Survival probability decreases ad age increases and it is smaller for men than for women.
# Survival probability, given AGE and MALE values, is \Phi(beta0+beta1*age+beta2*male),
# and therefore it is a function of the regression parameters; let us compute its posterior d., 
# for a man (male=1), as age varies.
# Coincide con la banda di previsione della sopravvivenza predittiva
#     di un "nuovo" paziente maschio.
a=seq(15,65)
X1=cbind(1,a,1)
#C'è una funzione del package che fa questo:
p.male=bprobit.probs(X1,fitBayes$beta)
# Disegnamo il grafico della mediana a posterior di queste quantità (al variare dell'età),
# insieme ai quantili di ordine 0.05 e 0.95
par(mfrow=c(1,1))
plot(a,apply(p.male,2,quantile,.5),type="l",ylim=c(0,1), xlab="age",ylab="Probability of Survival")
lines(a,apply(p.male,2,quantile,.05),lty=2)
lines(a,apply(p.male,2,quantile,.95),lty=2)



################################################################
# Model choice - we compute the log of the marginal density of the data for 4 different models
# and make 2-by-2 comparisons.
# Choose the model corresponding to the highest marginal (NON è un criterio OTTIMO):

# Modello 1: tutte e tre le covariate (compresa l'intercetta)
bayes.probit(y,X,1000,list(beta=beta00,P=P0))$log.marg

# Modello 2: intercetta, male
bayes.probit(y,X[,-2],1000,
   list(beta=beta00[-2],P=P0[-2,-2]))$log.marg

# Modello 3: intercetta, age
bayes.probit(y,X[,-3],1000,
   list(beta=beta00[-3],P=P0[-3,-3]))$log.marg

# Modello 4: solo intercetta
bayes.probit(y,X[,-c(2,3)],1000,
   list(beta=beta00[-c(2,3)],P=P0[-c(2,3),-c(2,3)]))$log.marg
#################################################################
## Model 1 is the best according to this criterion 

###############################################################################################
#### Package MCMCpack  http://mcmcpack.berkeley.edu/
#### Paper on JStatSoft http://www.jstatsoft.org/article/view/v042i09
#### Better download the package from https://cran.r-project.org/web/packages/MCMCpack/index.html 
#### In alternativa: arm, MCMCglm
#
# MODELS implemented in the package: linear regression (with Gaussian errors), 
# a hierarchical longitudinal model with Gaussian errors, 
# a probit model, a logistic regression model, a one-dimensional item response theory model, 
# a K-dimensional item response theory model, a normal theory factor analysis model, 
# a mixed response factor analysis model, an ordinal factor analysis model, 
# a Poisson regression, a tobit regression, a multinomial logit model, 
# a dynamic ecological inference model, a hierarchial ecological inference model,
# and an ordered probit model. 
# The package also contains densities and random number generators for 
# commonly used distributions that are not part of the
# standard R distribution, a general purpose Metropolis sampling algorithm, 
# and some utilities for visualization and data manipulation.

#####################################################################################
###### https://cran.r-project.org/web/views/Bayesian.html  is a site for R packages 
######         for Bayesian inference
#####################################################################################

library(MCMCpack)
?MCMCprobit

library(LearnBayes)
data(donner)
y=donner$survival
X=cbind(1,donner$age,donner$male)

## mlfit <- glm(y ~ X, family = binomial(link = probit), x = TRUE, y = TRUE)

# Input: response vector, covariate (design) matrix and number of simulations to run
# Output: MCMC sample forom the posterior distribution of beta 
# The fitted prior is Gaussian, beta ~ N_p(b0,B0^(-1))
## Here we use Zellner g prior: beta ha distribuzione Gaussiana con media beta00 e matrice di varianza c*(X^T X)^(-1)
 
beta00=c(0,0,0); c0=100 #la prior ha un peso dell'1% rispetto al campione
## Provate con c0=1 e c0=10
P0=t(X)%*%X/c0

bfit <- MCMCprobit(y ~ donner$age +as.factor(donner$male), b0=beta00, B0=P0,marginal.likelihood="Chib95")
 # CAREFUL: here B0 denotes the PRECISION matrix for beta

bfitout<-as.matrix(bfit)

apply(bfitout,2,mean)
apply(bfitout,2,sd)


par(mfrow=c(1,3))
plot(density(bfitout[,1]), xlab='beta0 - intercept',main=' ') 
plot(density(bfitout[,2]), xlab='beta1 - age',main=' ') 
plot(density(bfitout[,3]), xlab='beta2 - male',main=' ') 

