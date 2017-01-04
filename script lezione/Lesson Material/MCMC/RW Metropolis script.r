######################################################
# EXAMPLE 1 - "simple" Metropolis-Hastings algorithm #
###################################################### 
# Random Walk Metropolis-Hastings chain 
#   with invariant distribution = N(0,1);
# proposal density: q(x,y)=f(y-x), where f(u)=1/(2alpha)I_(-alpha,alpha)),
#   i.e.          Y= x+U(-alpha,alpha).
# The state space is E=\Re, an open connected set, and 
#   f(u)>0 in U(0).
# We'll run the algorithm under different values for alpha, e.g. alpha=1,10,100
# The chain starts at 0
# Tipically, for large values of alpha, the acceptance ratio will be small, 
#   i.e. I'll reject the candidate point often, 
# since I'm sampling the candidate point from a region where 
#   the exact posterior density is very low  

norm<-function (n, alpha) 
{
        vec <- vector("numeric", n+1)
        x = 0
#        x <- -2*alpha  # RUN the algorithm by yourself with this value as the initial point!
        accept=0
        vec[1] <- x
        for (i in 2:n) {
                innov <- runif(1, -alpha, alpha)
                y <- x + innov
                aprob <- min(1, dnorm(y)/dnorm(x))
                u <- runif(1)
# If condition "(u < aprob)" is NOT met, we'll skip command "x <- y", 
#  so that the MC does not move from x 
                if (u < aprob) 
                     {
                        x <- y
                        accept1=accept+1
                        accept=accept1
                     }
                vec[i] <- x
        }

        vec[n+1]= accept/n
        
        vec
}

# Number of iterations of the MC
k = 10000
##################### alpha=1
normvec<-norm(k,1)
normvec1 = normvec[1:k]
rate1 = normvec[(k+1)]
par(mfrow=c(1,2))
plot(ts(normvec1))
hist(normvec1,30, prob=T)
val=seq(-3,3,0.1)
points(val,dnorm(val),type='l',col='red')

x11()
#plot(density(normvec1),bw='nrd')
#points(val,dnorm(val),type='l',col='red')

burnin=1000 # usually we get rid of the first burnin iterations, 
            # where the chain hasn't reached stationarity yet 
# The ergodic Th says we can use all iteration from any initial point x. However, 
# we are averaging some function of some finite number of samples, so that out average will 
# be a better approximation if we start at a typical point in the density we are sampling from
# and if this density is closer to the target distribution
b1=burnin+1
normveclast1=normvec[b1:k]

par(mfrow=c(1,1))
plot(density(normveclast1))
points(val,dnorm(val),type='l',col='red')

rate1

x11()
################## alpha=10
normvec<-norm(k,10)
normvec10 = normvec[1:k]
rate10 = normvec[(k+1)]
par(mfrow=c(1,2))
plot(ts(normvec10))
hist(normvec10,30,prob=T)
points(val,dnorm(val),type='l',col='red')


#windows()
#normveclast10=normvec10[b1:k] # traceplot discarding the first burnin iterations 
#plot(density(normveclast10))
#points(val,dnorm(val),type='l',col='red')

rate10

################## alpha=100
x11()
normvec<-norm(k,100)
normvec100 = normvec[1:k]
rate100 = normvec[(k+1)]
par(mfrow=c(1,2))
plot(ts(normvec100))
hist(normvec100,30,prob=T)
points(val,dnorm(val),type='l',col='red')

#windows()
#normveclast100=normvec100[b1:k]  # traceplot discarding the first burnin iterations 
#plot(density(normveclast100))
#points(val,dnorm(val),type='l',col='red')

rate100

## COMPARISON among kernel density plots of the simulated MCs for different values of alpha
windows()
par(mfrow=c(1,3))
plot(density(normvec1),main="alpha=1")
points(val,dnorm(val),type='l',col='red')
plot(density(normvec10),main="alpha=10")
points(val,dnorm(val),type='l',col='red')
plot(density(normvec100),main="alpha=100")
points(val,dnorm(val),type='l',col='red')


## COMPARISON among traceplots of the simulated MCs for different values of alpha
windows()
par(mfrow=c(1,3))
plot(ts(normvec1[1:1000]))
plot(ts(normvec10[1:1000]))
plot(ts(normvec100[1:1000]))

## COMPARISON among acceptance rates of the simulated MCs for different values of alpha
rate1;rate10;rate100


# end


##################################################################
# EXAMPLE 2 - Ex in 6.7 in Albert
##################################################################
# Section 6.7 Learning about a Normal Population from Grouped Data
##################################################################
#Example of a MH algorithm
library(LearnBayes)
# Data are heights (in inches) of men from a local college 
# The interest is in making inference about the mean and sd of the heights; 
# X_i \sim n(\mu,\sigma^2) 
# 1 inch= 2.54 cm 
# However data are grouped: out of n=211 men, n_1=14 of them have height less than 66, 
# n_2=30 have height between 66 and 68, ect.  

# d is a list defining the intervals and the frequencies
d=list(int.lo=c(-Inf,seq(66,74,by=2)),int.hi=c(seq(66,74,by=2), Inf),f=c(14,30,49,70,33,15))
d


help(groupeddatapost)
# it computes the log posterior density of (\mu,log \sigma) for normal sampling 
#       where the data is observed in grouped form
# Laa prior per (mu,log(sigma)) ? IMPROPRIA, proporzionale ad una cost.
# In pratica, ? un modello multinomiale, in cui il parametro (= probabilit?) vettoriale
# ? la massa assegnata da una N(mu,sigma^2) ad una partizione di 6 intervalli

#INPUT: theta=valore di mu e log(sigma) in cui calcolo la log-posterior
#       data= dataframe dei dati in forma raggruppata
#OUTPUT: log-posterior density 

# EX: compute the log-posterior density at  mu=70, log(sigma)=1
groupeddatapost(c(70,1),d)


help(laplace)
#################################
# It is an iterative method (Newton's method) to compute, for a general posterior density,
#   an estimate of the posterior mode and the associated variance-covariance matrix, 
# INPUT: the log-posterior density, the initial point of the iterative method, 
#         parameters of the log-posterior density (the vector d in this case) 

start=c(70,1) #punto iniziale - costruito da dati "fittizi"

fit=laplace(groupeddatapost,start,d)

fit
# fit$var is the matrix  (\tilde I_n)^{-1}, where \tilde I_n is the  
#   "generalized observed Fisher information matrix"  (-d^2/ d d log(pi(mu,log(sigma)|x)) #d qui indica la DERIVATA!
# Approximately the marginal posterior density of mu is N(70.17,0.035) and 
#     the marginal posterior density of  lambda=log(sigma) is  N(0.974, 0.003)

diag(fit$var)
modal.sds=sqrt(diag(fit$var))

####################### METROPOLIS-HASTINGS  ALGORITHM ################################
# Build the MC according to a RANDOM WALK Metropolis algorithm:
# the bivariate proposal density q is N(val_prec,c^2 V)
# where the scale factor c=2, V  is the covariance matrix (\tilde I_n)^{-1};
# this is the covariance of the Gaussian approximation (centered at the posterior mode)
# The PROPOSAL distribution has to be close to the true posterior, i.e.
#   it should be an approximation of the posterior!!! 
# The scale factor c determines the behaviour of the MC 
######################proposal=list(var=fit$var,scale=2)
proposal=list(var=fit$var,scale=2) #2
proposal


# Proposal distributions: N(0,\Sigma, Uniform on the ball of radius R, t-Student 
# (or more general mixtures of normal misture);
# usually, variance matrix of the Gaussian proposal is
# c^2 * inverse of the "generalized observed Fisher information" matrix
# for  c >0; here  c =1/2, 1, 2 work "well"

# La funzione del package LearnBayes ? rwmetrop
# Simulates iterates of a random walk Metropolis chain for an arbitrary real-valued posterior density defined by the user
#help(rwmetrop)
#print(rwmetrop)
#INPUT: log-posterior = groupeddatapost
#       proposal = N(0,c^2*V)
#       start= initial point of the MC
#       number of iterations=10000
#       d=parametri della log-posterior = iperparametri + dati
# the acceptance probability alpha(x,y) is the ratio of the posterior densities (random walk MH)

fit2=rwmetrop(groupeddatapost,proposal,start,10000,d)

#OUTPUT : 
# par = a matrix of simulated values where each row corresponds to a simulated value of the vector parameter
# accept= the acceptance rate of the algorithm

fit2
fit2$accept

# The acceptance rate is good? 
# The proposal should be more "diffuse" than posterior density , in order to explore 
# the whole support of the posterior density. 
# If c is too small, the acceptance rate will be high; basically, a small scale c 
# means that I'm sampling the candidate point y from a density with small variance,
# so that y will be close to x, and therefore alpha(x,y)=pi(y)/pi(x) will assume values close to 1.
# This means that the MC explores the state space weakly, i.e. the MC moves slowly through 
# the state space, and the algorithm can turn out to be inefficient. 
# On the other hand, if c is too large, alpha(x,y) will often be small, and may candidate points will be rejected,
# so that the overall acceptance rate r will be too small, and also in this case the algorithm is inefficient.
# Il tasso di accettazione BUONO per un RW-MH deve essere compreso in [0.25, 0.45] circa, o [0.23,0.5). 
# Il valore 0.23... ? il limite quando la dimensione dello spazio converge a infinito.

# For a random walk MH, remember that \alpha(x,y)=\pi(y)/\pi(x) indicates how probable the new proposed sample y
# is with respect to the current sample x.  If we attempt to move to a point that is more probable 
# than the existing point (i.e. a point in a higher-density region of \pi), 
# we will always accept the move. However, if we attempt to move to a less probable point, 
# we will sometimes reject the move, and the more the relative drop in probability, the more likely 
# we are to reject the new point. 

# Per quanto riguarda independence chain Metropolis, tassi di accettazione che
# portano ad algoritmi efficienti sono decisamente pi? alti; teoricamente,
# pi? ? alto il rapporto di accettazione e migliore ? il tasso di convergenza 
# della catena.

#Calcolo la media ergodica componente per componente
#help(apply)
#Medie a posteriori (calcolate col metodo MCMC) di mu e log(sigma)
post.means=apply(fit2$par,2,mean)
#Standard deviation delle 2 distribuzioni a posteriori marginali (MCMC)
post.sds=apply(fit2$par,2,sd)
post.means
post.sds

#Confronto fra medie e varianze dell'approx gaussiana con quelle a posteriori
cbind(c(fit$mode),modal.sds)
cbind(post.means,post.sds)
#In questo caso l'approssimazione col metodo di Laplace ? buona


#Disegnamo un countour plot della posterior esatta (nota a meno della costante di normalizzazione)
#  e confrontiamo con i punti della MC dal 5001 in poi
mycontour(groupeddatapost,c(69,71,.6,1.3),d)
points(fit2$par[5001:10000,1],fit2$par[5001:10000,2])

#Provate a cambiare scale alla proposal, per es. scale=0.1 
# e scale=10

windows()
## traceplots delle due marginali
par(mfrow=c(1,2))
plot(ts(fit2$par[5001:10000,1]),ylab="mu")
plot(ts(fit2$par[5001:10000,2]),ylab="sigma")


##################################################
# Section 6.8 Example of Output Analysis
##################################################
# MONITORING the CONVERGENCE of the MC to the STATIONARY distribution
# If we simulate realizations from a MC, then under 'broad' conditions, EVENTUALLY
# the simulated values will be marginally distributed acconding to the invariant distr.
# and the ergodic means will converge to correspondent integrals of the invariant distribution.
# However, initial points of the MC, if included in the computation of the ergodic mean, 
# may influence the value of the ergodic mean itself, yielding a poor approximation. 
# For this reason we usually discard the first non-stationary portion of the chain, 
# called BURN-IN.
# Therefore we need to be able to detect when the marginal behaviour of the MC is 
# sufficiently close to stationarity, and harvest all subsequent realizations as a DEPENDENT sample 
# from the stationary distribution.  
# We diagnose convergence retrospectly, by guessing for how long to run the simulation and
# then trying to determinate if some latter portion of the chain can be considered stationary.
# The sample size required depends on the EFFICIENCY of the MC; moreover this sample size 
# increases with an increasing level of serial dependence of the chain (autocorrelation of the MC).

library(LearnBayes)
d=list(int.lo=c(-Inf,seq(66,74,by=2)),
        int.hi=c(seq(66,74,by=2), Inf),
        f=c(14,30,49,70,33,15))

library(coda)
library(lattice)
start=c(70,1)
fit=laplace(groupeddatapost,start,d)

# Let us choose an UNLUCKY initial point for the chain, and a small scale factor for the proposal c=0.2; 
start=c(65,1)
proposal=list(var=fit$var,scale=0.2)
bayesfit=rwmetrop(groupeddatapost,proposal,start,10000,d)

bayesfit$accept
apply(bayesfit$par,2,mean)
apply(bayesfit$par,2,sd)


#Assegno i nomi alle colonne di bayesfit$par
dimnames(bayesfit$par)[[2]]=c("mu","log sigma")
xyplot(mcmc(bayesfit$par),col="black")
x11()
# BURN-IN: about 600 iterations
xyplot(mcmc(bayesfit$par[1:1000,]),col="black")
# If we forget the first 2000 iterations
xyplot(mcmc(bayesfit$par[-c(1:2000),]),col="black")
# Sembra che gi? dopo le prime 2000 iterazioni si sia giunti a convergenza
# a 'good ' traceplot is a fat hairy caterpillar !
# a 'snake-like' caterpillar does NOT mean in general the the MC hasn't reached convergence,
#    but we need to increase the sample size because there is high correlation 

# For a given sample size, the ACCURACY of our inferences is dependent on the EFFICIENCY 
# of our posterior samples. We'll see that the efficiency decreases with an increasing level 
# of AUTOCORRELATION.  

par(mfrow=c(2,1))
# Plots of AUTOCORRELATION between X_n and X_{n+k} by CODA
autocorr.plot(mcmc(bayesfit$par[-c(1:2000),]),auto.layout=FALSE)
#Dall'HELP della funzione "autocorr": High autocorrelations within chains 
#indicate slow mixing (the support of the target distribution has not been fully explored yet)
#and, usually, slow convergence. 
#It may be useful to 
#thin out a chain with high autocorrelations before calculating summary 
#statistics: a thinned chain may contain most of the information, but take up 
#less space in memory. Re-running the MCMC sampler with a different 
#parameterization may help to reduce autocorrelation

summary(mcmc(bayesfit$par[-c(1:2000),]))
# We estimate the (square root of the) variance of the MCMC estimator, that is sigma^2_f / T
# col metodo delle batch means:
batchSE(mcmc(bayesfit$par[-c(1:2000),]), batchSize=50)
# Si ottiene un valore diverso da Naive SE e da Time-series SE !!

effectiveSize(mcmc(bayesfit$par[-c(1:2000),])) # su 8000 iterazioni

# An estimate of the square root of sigma^2_f / T is called Monte Carlo standard error



# In order to check the convergence of the MC :
# - analysis of the traceplots (also to understand which burn-in we should consider). 
# - MCerror/posterior sd has to be small (less than 0.1 o 1 o 5%);
#      l'MCerror viene stimato col metodo delle batch means o quello usato in "summary".
# - convergence diagnostics
# In order to check the mixing of the chain:
# - autocorrelation plot: the autocorrelation should be 'small' as the lag increases.  
#   In fact, if there is a strong dependence between X_n and X_{n+k}, this means
#   that the value of X_{n+k} strongly depends on the previous value X_n, and therefore
#   the chain does not mix 'well', since the chain would be bound to stay in some region 
#    of the support set, and not 'covering' the whole support of the target distribution. 
# Slow mixing usually denotes a 'slow' convergence:
# one solution is to reparametrize the parameters in order to reduce autocorrelation. 
# the other solution (more popular!) is to perform a process known as THINNING whereby only 
# every kth value from the MC is actually stored for inference, so that autocorrelation among 
# successive value of the chain is diminished. 
# Note that this only represents an efficiency gain in terms of storing and post-processing 
# the sample: for the same computational cost of simulation, the full sample size will always
# contain more information. 
# Il numero di iterazioni da salvare (al meno del burnin e del thinning)
# per calcolare poi le medie ergodiche di stabilisce in base al MCerror 
# che si vuole ottenere.
# L'Effective Sample Size (ESS) di una componente unidimensionale  di una MCMC ? 
# il numero di iterazioni indipendenti che dovrei effettuare dalla stessa distribuzione a posteriori
# per ottenere lo stesso MCerror; se l'ESS ? piccolo, questo indica che la 
# stima sar? 'povera', perch? la standard deviation del parametro sar? grande 
# per tenere costante l'MCerror - pi? ? grande l'ESS e meglio ?...


#Ripeto l'algoritmo di M-H con un "migliore" punto iniziale e una migliore "scale"
start=c(70,1)
proposal=list(var=fit$var,scale=2.0)
bayesfit2=rwmetrop(groupeddatapost,proposal,start,10000,d)

bayesfit2$accept
apply(bayesfit2$par,2,mean)
apply(bayesfit2$par,2,sd)

dimnames(bayesfit2$par)[[2]]=c("mu","log sigma")
xyplot(mcmc(bayesfit2$par),col="black")
x11()
xyplot(mcmc(bayesfit2$par[1:30,]),col="black")
sim.parameters=mcmc(bayesfit2$par[-c(1:2000),])
xyplot(mcmc(bayesfit2$par[-c(1:2000),]),col="black")

par(mfrow=c(2,1))
autocorr.plot(sim.parameters,auto.layout=FALSE)
summary(sim.parameters)
batchSE(sim.parameters, batchSize=50)
#L'ultima colonna di "summary" e il comando batchSE stimano entrambi 
#la radice quadrata di (sigma_f)^2/n che e' la varianza asintotica dell'
#errore Monte Carlo, ma batchSE ? una stima MENO precisa della stima "time series"
 
#Kernel density delle posterior marginali 
densplot(mcmc(bayesfit2$par[-c(1:2000),1]))
densplot(mcmc(bayesfit2$par[-c(1:2000),2]),ylim=c(0,8))


effectiveSize(mcmc(bayesfit2$par[-c(1:2000),]))
