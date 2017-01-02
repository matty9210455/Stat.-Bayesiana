
#Per essere sicuri di cancellare tutto!
rm(list=ls())

########## Here I use JAGS!!!!
#Load the library
library(rjags)   # to interface R with JAGS 




## Rheumatoid arthritis is an autoimmune disease characterized by 
## chronic synovial inflammation and destruction of cartilage and bone 
## in the joints. The Rotterdam Early Arthritis Cohort (REACH) study was
## initiated in 2004 to investigate the development of rheumatoid arthritis
## in patients with early manifestations of joint impairment.
## Information regarding basic patient characteristics, serological 
## measurements, and patterns of disease involvement at baseline has 
## been gathered in 681 recruited patients. It is of interest
## to know which of the following 12 factors are potentially 
## associated with the development of rheumatoid arthritis considered as
## a binary (yes/no) outcome:
## ACCP (cyclic citrullinated peptide antibody), 
## age,
## ESR (erythrocyte sedimentation rate), 
## DC (duration of complaints in days),
## stiffness (duration of morning stiffness in minutes), 
## RF (rheumatoid factor),
## gender,
## Sym (symmetrical pattern of joint inflammation;yes/no), 
## SJC (swollen joint count), 
## TJC (tender joint count), 
## BCPH (bilateral compression pain in hands; yes/no), 
## BCPF (bilateral compression pain in feet; yes/no).

##The standard approach to analyze these data would be to use 
## logistic/probit regression combined with some off-the-shelf variable
## selection method. The F -to-out backward selection with p D 0:05 yields a
## model with the following variables: ACCP, ESR, DC, Sym, SJC, and BCPH.
## The model with the most favorable value of the AIC selected after an 
## exhaustive model evaluation contains two extra variables: RF and stiffness. 
## Which of these models provide the best approximation to the true 
## underlying relationships is, if at all possible, difficult to assess. 



reach <- read.table("REACH_data.txt",header=T)
#the design matrix
x <- as.matrix(reach[,1:12])
## the response vector
Y <- as.vector(reach[,13])




#number of observations
N <- dim(x)[1]
#number of covariates
p <- dim(x)[2]



###Parameters of the  quasi-spike slab prior
#We fix the shape of the inverse gamma to be 2 (finte mean)
# while, as suggestd Iswaran Rao (2003) jasa, we fix v1=1

a=2
v1=1
#We consider setting 2 of hyperparameters
var_sp <- 0.00027
var_sl  <- 2.7

##Then as a conseguence, 
## Remark if X~t-student with 
## scale parameter sigma and nu>2 df
##  var(X)=sigma*(nu/(nu-2)) 

##Then
b <- var_sl*(a-1)/v1
b
v0 <- var_sp*(a-1)/b 
v0


##### delta
s1 <- sqrt(b*v0/a)
s2 <- sqrt(b*v1/a)



rr <- (s2/s1)^(2/(2*a+1))
dd <- sqrt(2*a*(1-rr)/(rr/s2^2-1/s1^2 ) )

dd
###########################################




library("MCMCpack")

#  In the next lines, I will use some  plot-tricks  of R, 
# to draw the quasi-spike and slab distribution mistrures  we will usre
# as prior  in this code 


curve(0.5*dinvgamma(x,shape=a,scale=v1*b)+0.5*dinvgamma(x,shape=a,scale=v0*b),from=0,to=10,type="n",ylab="prior density")

curve(0.5*dinvgamma(x,shape=a,scale=v0*b),from=0.000001,to=0.2,add=T,col="red",lty=2,lwd=2)
text(2,0.05,"quasi-spike component",col="red")

curve(0.5*dinvgamma(x,shape=a,scale=v1*b),from=0,to=10,add=T,col="blue",lty=2,lwd=2)
text(7,0.05,"slab component",col="blue")





#############################################
### To plot the quasi-spike and slab prior we need the 
### distribution of a t-student with scale parameter sigma
### and center mu

tdist <- function(x,nu,sig=1,mu=0){
	out <- lgamma((nu+1)/2)-lgamma(nu/2)-1/2*log(pi*nu)-log(sig)-((nu+1)/2)*log( 1+1/nu*( (x-mu)/sig)^2)
	return(exp(out))
}

##Then we compute the scale parameter for the quasi-spike and the slab
## prior
s1 <- sqrt(b*v0/a)
s2 <- sqrt(b*v1/a)

##And we try to visualize
curve(0.5*tdist(x,nu=2*a,sig=s1)+0.5*tdist(x,nu=2*a,sig=s2),from=-5,to=5,type="l",ylab="prior density")

curve(0.5*tdist(x,nu=2*a,sig=s1),from=-5,to=5,add=T,col="red",lty=2,lwd=2)
text(1.8,7,"quasi-spike component",col="red")
curve(0.5*tdist(x,nu=2*a,sig=s2),from=-5,to=5,add=T,col="blue",lty=2,lwd=2)
text(-2,0.3,"slab component",col="blue")
curve(0.5*tdist(x,nu=2*a,sig=s1)+0.5*tdist(x,nu=2*a,sig=s2),from=-5,to=5,add=T)





###We will use jags to implement the model


data_win <-list(N=N, p=p, Y = Y,x=as.matrix(x),v0=v0,v1=v1,a=a,b=b)


## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
	list( beta0=1.0, beta=rep(1,p), tau=1.0,
	     g=rep(0,p),.RNG.seed=321,.RNG.name = 'base::Wichmann-Hill') 
}




#########################

# The function  jags.model() compile the jags model, moreover it performs
#and adaptive burn-in for 'nchain' parallel chains.
#As input we have to provide: the model, the data, 
#the initial condition of the chains, the number of chains, and the number of
#adaptations (burn inn).

model=jags.model("NMIG_probit.bug",data=data_win,n.adapt=1000,inits=inits,n.chains=1) 






####Posterior parameter WinBUGS  has to save 
param <- c("beta0","beta","sigma2","g","mdl")


## Some other information for the 
## MCMC algorihm
#number of iterations


###The probit model is computationally more intensive
## I run a chain previusly with the usual choice of
nit <- 50000
#thinning
thin <-10

#to be fast here
#nit=100
#thin=1




##The command coda.samle() calls jags from R passing the data and initial
#value just defined

output_nmig=coda.samples(model=model,variable.names=param,n.iter=nit,thin=thin)

#save(output_nmig,file="nmig_set1.dat")

load("nmig_set1.dat")
# the output is an mcmc object of the library coda

str(output_nmig)

### A qalitative analysis of the convergence of the posteror chais

plot(output_nmig)
### To work wit the posterior chain it is better 
### to cast output to be an array object of R

output_nmig <- as.matrix(output_nmig)

###Some variable selection thecniques:

# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)

##We save the posterior chain of the inclusion variable in post.g
colnames(output_nmig)
post.g <- output_nmig[,14:25]

# then we compute the sample mean , column by column
apply(post.g,2,"mean")
post_mean_g <- apply(post.g,2,"mean") 

## to produce a bar plot of the  posterior inclusion probabilities
#postscript("reach_nmig_set2.ps")
barplot(post_mean_g,names.arg=colnames(x),main="Posterior inclusion probabilities set 2")
#Finally we can add a horizontal line 
abline(h=0.5,col="red",lty=2,lwd=3)
#dev.off()
# The variable included in the median probability model are:
mp_model<-post_mean_g>0.5
#mp_model is a logical vector with TRUE if the corresponding
# parameters is included in to the model FALSE otherwise
mp_model

## Frequency of the selected variables
post_mean_g[mp_model]


##trick to see whech variable choose with the...
paste(colnames(x),1:p,sep="")[mp_model]


## close all the graphical device
dev.off()


######
######  Highest posterior density model (HPD)
# The following function has as input an integer and a  base p.
# As a result it returns a vector containing the representation in
# in the base p of the integer n.

##Note that n%%base indicates n mod base and %/% indicates integer division. 

as.binary <- function(n,base=2)
{ 
	if(n==0){return(0)}
	out <- NULL
	while(n > 0) {
		out <- c(n%%base , out)
		n <- n %/% base
       	}
	names(out) <- NULL
   return(out)
} 


##In the following function I will use the as.binary function
## and two time a paste() trick to write down automatically
## in a nice form wich variable are included in a visited model

##It have as input the number of the model n and the number of covariate n.cov

wich.variable <- function(n,n.cov){
	bin <- rev(as.binary(n))
	n.bin <- length(bin)
	logic.bin <- rep(FALSE,n.cov+1)
	for(i in 1:n.bin){
		logic.bin[i]=(bin[i]==1)
	}
	out <- paste(paste("x",0:p,sep="")[logic.bin],collapse="_")

	return(out)

}






### We start to analize how many models have been visited 
## by the posterior chain:
length(unique( output_nmig[,"mdl"]))

## Now we compute the posterior frequency of the visited  models
visited_models<-table(output_nmig[,"mdl"])
visited_models




#We can visualize the table of the visited models
x11()
barplot(visited_models)


## Let's sort the table to see which are the "top ten"  
top <- sort(visited_models,decreasing=T)[1:10]
top


## In the following lines we will use the as.binary with a "paste()" trick 
## to visualize wich variable is included in the most visited models
numeric_top <- as.numeric(names(top))

for(j in 1:10){
names(top)[j]=wich.variable(numeric_top[j],p)
}


top


 
#############################
#########  Th Hard Shrinckage (HS; hard thresholding/selection shrinkage)
#variables are included whenever the absolute value of the estimated 
# coefficient (e.g., posterior mean) exceeds some threshold value.

# Remark we will  base the interval decision criterion for HS on a 
# one standard deviation interval around the posterior mean.

#posterior of the beta parmeter


##We save the posterior chain of the inclusion variable in post.g
colnames(output_nmig)
post.beta <- output_nmig[,1:12]

# then we compute the sample mean , column by column
mean.beta.post <- apply(post.beta,2,"mean")
sd.beta.post <- apply(post.beta,2,"sd")

require(plotrix)

#postscript("HS_decision_nmig_set2.ps")
plotCI(x=1:p, y=mean.beta.post, uiw=sd.beta.post,lwd=1.5, main="Decision intervals for HS, nmig-set1")
abline(h=0,col="blue")
#dev.off()
##




####################



