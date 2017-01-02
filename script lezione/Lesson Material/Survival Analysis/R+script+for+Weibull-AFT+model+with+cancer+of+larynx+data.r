### Example 13.1.2 di BIDA: AFT,  Weibull distribution
### data: cancer of the larynx  

rm(list=ls())


########## 
#Load the library
setwd("C:/alessandra/Documents/Didattica/StatisticaBayesiana/mat1617")
library(rjags)   # to interface R with JAGS


in_data = read.table("data_larynxcancer.txt",header=TRUE)
in_data
attach(in_data)

### analisi frequentista -> ILARIA
library(survival)
# i tempi vengono trasformati in un oggetto "Surv"
# che include anche le censure a destra
tempi90=t[1:90] #90 datapoint: 50 observed, e altri 40 censurati a destra
cens90=cens[1:90]
tempi90[cens90 >0]=cens90[cens90 >0]
tempi90 # ora il vettore tempi90 contiene sia i dati esatti che quelli censurati

s_cens=1-cens90
s_cens[s_cens<0]=0
s_cens # ora ho definito un vettore degli indicatori di "censura": 
       # l'i-esimo valore =1 se il tempo è stato osservato esattamente 
       # e =0 se il tempo è censurato a destra

tempi<-Surv(tempi90,s_cens) # i tempi vengono trasformati in un oggetto "Surv"
                            # che include anche le censure a destra

#Stimatore di Kaplan e Maier, sull'intero
#insieme dei dati, no covariate:
fit0=survfit(tempi~1,data=in_data[1:90,])
plot(fit0)

#Stimatore di Kaplan-Meier considerando Yr come covariata
fit1=survfit(tempi~as.factor(Yr),data=in_data[1:90,])
summary(fit1)
#Se log=T, la funzione di sopravvivenza sull'asse delle y e' plottata
#in scala logaritmica
plot(fit1,log=F)
# NON è per niente chiaro l'effetto di Yr sulla sopravvivenza
barplot(table(Yr))
boxplot(t~as.factor(Yr))

# Infatti: 
plot(fit1,fun="cloglog",xlab='log-tempi di sopravvivenza',ylab='log H(t)',col=c(2,3,4,5,6,7,8,9,10))
##Per i dettagli sul'opszione fun vedere help plot.survfit
#cloglog" creates a complimentary log-log survival plot 
#(f(y) = log(-log(y)) along with log scale for the x-axis). 



#Stimatore di Kaplan-Meier considerando stage come covariata
fit2=survfit(tempi~as.factor(stage),data=in_data[1:90,])
summary(fit2)
#Se log=T, la funzione di sopravvivenza sull'asse delle y e' plottata
#in scala logaritmica
plot(fit2,log=F)
## da questo plot l'effetto di Stage è piu chiaro 
############################################################################
## Analisi del modello bayesiano AFT con jags

## data on 90 male patients with cancer of larynx: 50 survival times are observed, 40 are right censored
## survival times(in months) 
## COVARIATE: stage of the disease at diagnosis (1,2,3,4) (stage)
##            year of diagnosis (Yr)
##            age at  diagnosis (age)
## Extra 16 rows: covariates values for 16 NEW patients; we will compute their predictive distribution
## In the bug file we will standardize Age and Yr
## For people with the same age and year of diagnosis, we are interested in the median lifetimes for individuals
## in Stage 2,3 or 4 relative to Stage 1 
## We are also interested in the predictive 5-months survival for 16 NEW patients


# A list of initial value for the MCMC algorithm 
# that JAGS will implement
inits = function() {list(beta=c(2,-0.14,-0.5,-1.5,-0.2,0), tau=1)}
# Al contrario di WinBUGS, in JAGS l'inizializzazione viene fatta automaticamente
# ATTENZIONE: qualche volta l'inizializzazione "automatica" NON funziona!

 
modelAFT=jags.model("AFT_larynx_model.txt",data=in_data,n.chains=1) 
#modelAFT=jags.model("AFT_larynx_model.txt",data=in_data,inits=inits,n.adapt=50000,n.chains=1) 
# Faccio un update del modello per il numero di iterazioni specificato SENZA salvare nulla 
update(modelAFT,49000)


variable.names=c("beta", "sigma", "rm", "prob", "mu")
# Parameters to be monitored: beta (6 components), sigma, mu (linear predictor for each patient in the dataset)
#                       rm(vector of  exp(beta[i]) - relative median of two patients)
#                       prob (posterior prob that  beta[i]  >0)
n.iter=50000 
thin=10  

outputAFT=coda.samples(model=modelAFT,variable.names=variable.names,n.iter=n.iter,thin=thin)
# salvo l'intera catena dei parametri monitorati (si tratta di una lista mcmc)


#save(outputAFT,file='BIDA_AFT_output.Rdata')

#### ANALISI DELL'OUTPUT
library(coda)       # per leggere catena mcmc
library(plotrix)    # per fare plot CIs

##### CARICA le traiettorie della catena
#### load('BIDA_AFT_output.Rdata')

data.out=as.matrix(outputAFT) # trasformo il dataframe in matrice (piu' facile da maneggiare)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain

#############################################
##  Summary statistics for the posterior
#############################################
summary(data.out)

## I could use CODA, but her ewe use standard R tools 
par(mfrow=c(2,3))
acf(data.out[,'beta.1.'],lwd=3,col="red3",main="autocorrelation of beta1")
acf(data.out[,'beta.2.'],lwd=3,col="red3",main="autocorrelation of beta2")
acf(data.out[,'beta.3.'],lwd=3,col="red3",main="autocorrelation of beta3")
acf(data.out[,'beta.4.'],lwd=3,col="red3",main="autocorrelation of beta4")
acf(data.out[,'beta.5.'],lwd=3,col="red3",main="autocorrelation of beta5")
acf(data.out[,'beta.6.'],lwd=3,col="red3",main="autocorrelation of beta6")

par(mfrow=c(2,3))
plot(ts(data.out[,'beta.1.']),xlab="t",ylab="beta1")
plot(ts(data.out[,'beta.2.']),xlab="t",ylab="beta2")
plot(ts(data.out[,'beta.3.']),xlab="t",ylab="beta3")
plot(ts(data.out[,'beta.4.']),xlab="t",ylab="beta4")
plot(ts(data.out[,'beta.5.']),xlab="t",ylab="beta5")
plot(ts(data.out[,'beta.6.']),xlab="t",ylab="beta6")

# Kernel density plot delle marginal posterior dei beta
par(mfrow=c(2,3))
plot(density(data.out[,'beta.1.'],adj=2),  xlab=expression(beta[1]),main="")
abline(v=quantile(data.out[,'beta.1.'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.2.'],adj=2),  xlab=expression(beta[2]),main="")
abline(v=quantile(data.out[,'beta.2.'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.3.'],adj=2),  xlab=expression(beta[3]),main="")
abline(v=quantile(data.out[,'beta.3.'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.4.'],adj=2),  xlab=expression(beta[4]),main="")
abline(v=quantile(data.out[,'beta.4.'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.5.'],adj=2),  xlab=expression(beta[5]),main="")
abline(v=quantile(data.out[,'beta.5.'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(data.out[,'beta.6.'],adj=2),  xlab=expression(beta[6]),main="")
abline(v=quantile(data.out[,'beta.6.'],prob=c(.025,.975)),lwd=2,col="red")

# IC a posteriori dei parametri beta_1, ...,beta_6
apply(data.out[,1:6],2,"quantile",prob=c(0.025,0.5,0.975))


windows()
par(mfrow=c(1,1))
plot(ts(data.out[,'sigma']),ylab="sigma")

par(mfrow=c(1,2))
plot(density(data.out[,'sigma'],adj=2),  xlab=expression(sigma),main="")
abline(v=quantile(data.out[,'sigma'],prob=c(.025,.975)),lwd=2,col="red")
plot(density(1/data.out[,'sigma'],adj=2),  xlab=expression(alpha),main="Parametro di ageing della Weibull")
abline(v=quantile(1/data.out[,'sigma'],prob=c(.025,.975)),lwd=2,col="red")
# marginal posterior of alpha is concentrated on alpha values >1, pointing out that an IFR Weibull distribution 
#   "explains" well the data 

# rm[2]=exp(beta[2]): relative median  of a patient with Stage 2 wrt to a patient (same age and Yr) with Stage 1, 
# i.e. ratio of medians of two individuals with same age and Yr, but one with Stage 2 (numerator) and the other (denom) with Stage 1  
# rm[3]=exp(beta[3]) è la mediana relativa (rapporto) di 2 pazienti 
# con stessa age e Yr, ma il primo ha Stage 3 e il secondo ha Stage 1
# rm[5]=exp(beta[5]): relative median of two individuals who are 1 (sample) standard deviation apart in age (=10.8years),
#             but with the same stage and year of diagnosis 
rmchain=cbind(rm_2=data.out[,'rm.2.'],rm_3=data.out[,'rm.3.'],
         rm_4=data.out[,'rm.4.'],rm_5=data.out[,'rm.5.'],rm_6=data.out[,'rm.6.'])
apply(rmchain,2,"quantile",prob=c(0.025,0.5,0.975))
# se andiamo ad analizzare la posterior di rm[2], questa è a cavallo del valore 1:
# non c'è evidenza che Stage 2 produce un funzionale mediana che è più piccolo di quello di Stage 1;
# c'è invece più certezza che un paziente in Stage 4 abbia un funzionale mediana che è più piccolo 
# di quello di Stage 1, a parità delle altre covariate. 
# Analogamente, un paziente con Yr pari a +1 sd empirica della covariate Yr ha funzionale mediana 
# che è più piccolo del paziente di riferimento


# Posterior probability that  beta_i>0
probchain=cbind(data.out[,'prob.2.'],data.out[,'prob.3.'],
               data.out[,'prob.4.'],data.out[,'prob.5.'],data.out[,'prob.6.'] )
apply(probchain,2,mean) # sono significative quelle covariate che corrispondono a valori
                            # molto vicino a 0 o molto vicini ad 1
apply(probchain,2,sd)
# C'è incertezza sulla significatività di beta_1 (riferimento di Stage 1)
# e su beta_2; gli altri beta_i sono più significativi, soprattutto beta_4 e beta_6 che
# non contenevano (beta_4 quasi) lo 0 nell'IC di livello 95%.
# La posterior di beta_6 è concentrata su valori negativi, in accordo con i dati, per cui
# la mediana dei tempi di sopravvivenza era più alta per Yr = 72 e 73 (e per 70 che ha solo 2 dati)


### Predictions for "new" patients  
### Covariate contenute nelle ultime 16 righe (dalla riga 91 alla riga 106) di data_in
### Stage = 1,2,3,4; Age = 50,70; Yr = 71,77
in_data[91:106,]

names(data.out)
mu=data.out[, 7:112] #sono i predittori lineari di "tutti" i 106 pazienti
     #in realtà, mi interessano solo quelle degli ultimi 16 pazienti per fare previsione
S=matrix(5000,16)
med=matrix(5000,16)

## (legge predittiva della) probabilità di sopravvivenza a 
##     5 mesi di "nuovi" pazienti
S = exp(-log(2)*exp(( log(5) - mu[,91:106])/data.out[,'sigma']))
## (legge predittiva dei) funzionali mediana di "nuovi" pazienti
med <- exp(mu[91:106])  

apply(S,2,"quantile",prob=c(0.025,0.5,0.975)) # il valore 50% è una stima della 
                         # PREDICTIVE 5-months survival probabilities for the new patients 
                                           

apply(med,2,"quantile",prob=c(0.025,0.5,0.975)) # IC a posteriori del funzionale mediana per i nuovi pazienti


#######################################################################################
######################################################################################

