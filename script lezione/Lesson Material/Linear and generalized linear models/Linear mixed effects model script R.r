################### EXAMPLE on Linear Mixed Models, HOFF's book - Sect 11.1, 11.3
#### FUNZIONI che mi servono dopo
### log-density of the multivariate normal distribution
ldmvnorm<-function(X,mu,Sigma,iSigma=solve(Sigma),dSigma=det(Sigma)) 
{
  Y<-t( t(X)-mu)
  sum(diag(-.5*t(Y)%*%Y%*%iSigma))  -
  .5*(  prod(dim(X))*log(2*pi) +     dim(X)[1]*log(dSigma) )
                                 
}
###


### sample from the multivariate normal distribution
rmvnorm<-function(n,mu,Sigma)
{
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 )
  {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  res
}
###


### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
     Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
     S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}
###

############################################################################################
########## DATA NESTED in GROUPS (hierarchy), LMM (linear mixed-effects model) ############
# LEGGO I DATI: 
# mathscore dataset of 10-th grade childrem from 100 different large urban public high schools
# the data concerns 1993 students
setwd("C:/alessandra/Documents/Didattica/StatisticaBayesiana/mat1516")
Y=read.table("school.mathscore.txt") #### ATTENZIONE!!! su Beep il nome dei file è diverso
                                     #### Questo è il dataset SENZA covariata 

### medie dei risultati PER GRUPPO
m = length(unique(Y[,1]))
n<-sv<-ybar<-rep(NA,m) 
for(j in 1:m) 
{ 
  ybar[j]<- mean(Y[Y[,1]==j,2])   #mean(Y[[j]])
  sv[j]<-var(Y[Y[,1]==j,2])
  n[j]<-sum(Y[,1]==j)
}

#########################################################################################################
## LOAD the dataset: covariate SES is an indicator of the student's family socioeconomic status 
odat=read.table("Dati_Ex_Sect11.3libroHoff_conCOVARIATA.txt") # in realtà questo dataset contiene school.mathscore.txt,
                                   # "school.ch11.txt" è più completo
                                   #### THIS IS the dataset INCLUDING covariate values - SEE the BEEP page

group<-odat$sch_id # school id 
unique(group) # 100 different schools

######################  COVARIATES   ###############################
#  X is a list of m=100 matrices: each matrix has a number of rows = number of students in the school (included in the sample) 
#  each matrix has a number of columns = 2: the first column contains 1's, the second contains the value 
#  of covariates ses (CENTRED wrt the group mean) of each student in the school

X<-list() ; 
for(j in 1:m) 
{
  xj<-odat$stu_ses[Y[,1]==j] #Y[,1] è la prima colonna, contiene l'indice della scuola 
  xj<-(xj-mean(xj))
  X[[j]]<-cbind( rep(1,n[j]), xj  )
}
####### la covariata è già centrata, e l'istruzione seguente, già commentata, NON SERVE
####ses_cen=odat$stu_ses-mean(odat$stu_ses)       
ses_cen=odat$stu_ses
summary(ses_cen)

### Scatteplot of the data (mathscores and SESs)
plot(odat$stu_ses,Y[,2],xlab='ses',ylab='mathscore')

### Least Squares Estimates within each group (j=1,...,100)
### perchè potrebbe essere che la relazione tra ses e voto dipenda dalla scuola
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[Y[,1]==j,2] ~ -1+X[[j]] ) # ci tolgo l'intercetta, perchè è già contenuta nelle matrici X
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
                } 

par(mar=c(2.75,2).75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
# LEFT panel: 100 different LS regression lines AND the average of these lines 
# (l'intercetta è la media di tutte le 100 intercette ai min quadrati, e analogamente la slope)
plot( range(ses_cen),range(Y[,2]),type="n",xlab="SES", 
   ylab="math score")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[2],lwd=2)

# MIDDLE panel
plot(n,BETA.LS[,1],xlab="sample size",ylab="intercept") # intercepts of the different regression lines vs 
                                                        # group (school) sizes
abline(h= BETA.MLS[1],col="black",lwd=2)

# RIGHT panel
plot(n,BETA.LS[,2],xlab="sample size",ylab="slope") #slopes of the different regression lines vs 
                                                        # group (school) sizes
abline(h= BETA.MLS[2],col="black",lwd=2)
## BLACK lines represent AVERAGE values (average line, average intercept, average slope)
## Schools with the highest sample sizes haave regression coefficients that are generally close to the average,
## whereas schools with extreme coefficients are generally those with low sample sizes.
## Bayesian solution: stabilize the estimates for small sample size schools by SHARING INFORMATION 
## ACROSS GROUPS, using a hierrarchical model.  
## Many regression lines show a positive slope, pointing out that, for many groups, as SES increases, 
## mathscore increases as well. 
## However, there are more than 15 schools with negative slope!
sum(BETA.LS[,2]<0)

dev.off()

####################################################################
## HIERARCHICAL REGRESSION MODEL - LMM LINEAR MIXED effects MODEL ##
####################################################################
p<-dim(X[[1]])[2]
# mu0, the prior expectation of theta is fixed equal to the average of the corresponding (frequentist)
# regression parameters  - similarly for the matrix Lambda0
#
theta<-mu0<-apply(BETA.LS,2,mean)
nu0<-1 ; s2<-s20<-mean(S2.LS)
eta0<-p+2 ; #così la prior per la matrice Sigma è diffusa
L0=matrix(nrow=2,ncol=2)
L0[1,1]=cov(BETA.LS)[1,1]
L0[1,2]=cov(BETA.LS)[1,2]
L0[2,1]=cov(BETA.LS)[2,1]
L0[2,2]=cov(BETA.LS)[2,2]

####  Inizializzazione della MC (Gibbs Sampler)
Sigma<-S0<-L0 #<-as.matrix(cov(BETA.LS)); 
BETA<-BETA.LS
THETA.b<-S2.b<-NULL
iL0<-solve(L0) ; iSigma<-solve(Sigma)
Sigma.ps<-matrix(0,p,p)
SIGMA.PS<-NULL
BETA.ps<-BETA*0
BETA.pp<-NULL
set.seed(1)
mu0[2]+c(-1.96,1.96)*sqrt(L0[2,2]) #prior IC per theta_2 (la slope media)


##### PARAMETERS: beta_1,...,beta_m, theta, SIGMA, sigma^2 (m=100)
##### dim(beta_i)=dim(theta)=2, SIGMA 2-by-2 matrix 
# sigma ^2 = variance of each response variable (we assume it constant to simplify)
#    la sua prior è una inv-gamma(nu0/2,sigma_0^2 nu0/2) 
#    con nu0=1 e sigma_0^2= media delle varianze campionarie nei gruppi 
# beta_1,...,beta,m |theta,SIGMA  iid N_2(\theta,SIGMA)
# theta ~ N_2(mu0,L0), con mu_0 pari al vettore delle medie nei gruppi delle stime ai minimi quadrati 
#                       delle intercette e delle slope, mentre L0 è la matrice delle covarianze empiriche 
#                       (nei gruppi) delle stime ai minimi quadrati delle intercette e delle slope
# SIGMA  ~ inv-Wishart(eta0,S_0^{-1}) con eta_0=4, così E(SIGMA)=1/(eta_0-p-1)* S_0= S_0 = matrice delle 
#                                                  covarianze empiriche delle stime ai minimi quadrati

##### Gibbs Sampler cycle
for(s in 1:10000) {  #10000 è un po' troppo, ci mette 5 minuti
  ##update beta_j 
  for(j in 1:m) 
  {  
    Vj<-solve( iSigma + t(X[[j]])%*%X[[j]]/s2 )
    Ej<-Vj%*%( iSigma%*%theta + t(X[[j]])%*%Y[Y[,1]==j,2]/s2 )
    BETA[j,]<-rmvnorm(1,Ej,Vj) 
  } 
  ##

  ##update theta
  Lm<-  solve( iL0 +  m*iSigma )
  mum<- Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum))
  theta<-t(rmvnorm(1,mum,Lm))
  ##

  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1, eta0+m, solve( S0+t(BETA-mtheta)%*%(BETA-mtheta) )  ) 
  ##

  ##update s2
  RSS<-0
  for(j in 1:m) { RSS<-RSS+sum( (Y[Y[,1]==j,2]-X[[j]]%*%BETA[j,] )^2 ) }
  s2<-1/rgamma(1,(nu0+sum(n))/2, (nu0*s20+RSS)/2 )
  ##
  ##store results
  if(s%%10==0) 
  { 
    cat(s,s2,"\n")
    S2.b<-c(S2.b,s2);THETA.b<-rbind(THETA.b,t(theta))
    Sigma.ps<-Sigma.ps+solve(iSigma) ; BETA.ps<-BETA.ps+BETA
    SIGMA.PS<-rbind(SIGMA.PS,c(solve(iSigma)))
    BETA.pp<-rbind(BETA.pp,rmvnorm(1,theta,solve(iSigma)) )
  }
  ##
}
## FINE CICLO
## Ho fatto un thinning di 10, quindi il mio final sample size è 1000


#save.image("data.f11-3")
#load("data.f11-3")

##### Convergence diagnostics
library(coda)
effectiveSize(S2.b)        # sigma^2
effectiveSize(THETA.b[,1]) #theta_1
effectiveSize(THETA.b[,2])  #theta_2

apply(SIGMA.PS,2,effectiveSize)  #\Sigma

par(mfrow=c(2,2))
tmp<-NULL;for(j in 1:dim(SIGMA.PS)[2]) { tmp<-c(tmp,acf(SIGMA.PS[,j])$acf[2]) }
dev.off()

acf(S2.b)
acf(THETA.b[,1])
acf(THETA.b[,2])

##################################################################################
# Plot of the marginal posterior of theta2, the population-slope, i.e. the fixed effect of covariate ses
#  and the posterior predictive distribution of a NEW school (not included in the data), but from the SAME 
#  population of schools. We sample from the "model$ of the beta parameters, when theta, Sigma are the 
#  current values in the chain. 
# This NEW parameter has been denoted by BETA.pp[,2] in the Gibbs sampler above

par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))

plot(density(THETA.b[,2],adj=2),xlim=range(BETA.pp[,2]), 
      main="",xlab="slope parameter",ylab="posterior density",lwd=2)
lines(density(BETA.pp[,2],adj=2),col="gray",lwd=2)
legend( -3 ,1.0 ,legend=c( expression(theta[2]),expression(tilde(beta)[2])), 
        lwd=c(2,2),col=c("black","gray"),bty="n") 

quantile(THETA.b[,2],prob=c(.025,.5,.975))
# 95% posterior CI of the 'average' slope theta_2 is very different from the prior CI that is (-3.86, 8.6)


# Posterior probability that theta2 >0 is very large, but it does not indicate that 
# any given within-school slope cannot be negative 

# Let us compute the probability that the slope tilde_beta2 (for a NEW school NOT INCLUDED in the dataset) 
#  is negative; this value is computed simulating tilde_beta2 from slope population distribution, when 
#  parameters are the current values of (\theta, \Sigma)
mean(BETA.pp[,2]<0) # è un valore piccolo, ma NON è nullo!

BETA.PM<-BETA.ps/1000
plot(range(ses_cen),range(Y[,2]),type="n",xlab="SES",
   ylab="math score")
for(j in 1:m) {    abline(BETA.PM[j,1],BETA.PM[j,2],col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,2]),lwd=2 )
### Posterior expectations of the 100 school -specific regression lines 
### i.e. the lines y=E(beta_{1,j}|dati) + E(beta_{2,j}|dati) * x_ij                         

# Guardate l'analoga figura con le 100 rette di regressione distinte (ai min quadrati)
# RIGHT PLOT: there is a shrinkage of the 'extreme' frequentist regression lines towards 
#  the across-group average line (black line); in fact, the Bayesian estimates we got are 
#  convex linear combinations of the (unique) prior mean line and the within-group frequentist estimates 
# Since we SHARED information across groups, hardly any of the slopes are negative!

#dev.off()

windows()
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot( range(ses_cen),range(Y[,2]),type="n",xlab="SES", 
   ylab="math score",main="Frequentist LS within-group regression lines")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[2],lwd=2)
plot(range(ses_cen),range(Y[,2]),type="n",xlab="SES",
   ylab="math score",main="Bayesian estimates - hierarchical model")
for(j in 1:m) {    abline(BETA.PM[j,1],BETA.PM[j,2],col="gray")  }
abline( mean(THETA.b[,1]),mean(THETA.b[,2]),lwd=2 )

#####################################################################################
