########## Here I use JAGS!!!!
#Load the library
library(rjags)   # to interface R with JAGS 


reach <- read.table("REACH_data.txt",header=T)
#the design matrix
x <- as.matrix(reach[,1:12])
## the response vector
Y <- as.vector(reach[,13])



#number of observations
N <- dim(x)[1]
#number of covariates
p <- dim(x)[2]


###Parameters of the  spike slab prior Set 1
c_ss<-10
intersect<-0.1
tau_ss<-intersect/sqrt(2*log(c_ss)*c_ss^2/(c_ss^2-1))



## With this choice of hyperameter, c_ss and intersection 
## the variance of the quasi-spike prior is
tau_ss^2
##While the variance of the slab is 
(tau_ss*c_ss)^2

#  In the next lines, I will use some  plot-tricks  of R, 
# to draw the quasi-spike and slab distribution mistrures  we will usre
# as prior  in this code 

curve(0.5*dnorm(x,mean=0,sd=tau_ss)+0.5*dnorm(x,mean=0,sd=tau_ss*c_ss),from=-5,to=5,type="n",ylab="prior density")


curve(0.5*dnorm(x,mean=0,sd=tau_ss),from=-5,to=5,add=T,col="red",lty=2,lwd=2)
text(-1.7,11,"quasi-spike component",col="red")
curve(0.5*dnorm(x,mean=0,sd=tau_ss*c_ss),from=-5,to=5,add=T,col="blue",lty=2,lwd=2)
text(2,0.5,"slab component",col="blue")
curve(0.5*dnorm(x,mean=0,sd=tau_ss)+0.5*dnorm(x,mean=0,sd=tau_ss*c_ss),from=-5,to=5,add=T,col="black",lwd=2)
text(1.5,4,"spike and slab prior",col="black")

#########################################





###We will use jags to implement the  ssvs model 


data_win <-list(N=N, p=p, Y = Y,x=as.matrix(x),tau_ss=tau_ss,c_ss=c_ss)

## A list of initial value for the MCMC algorithm 
# that WinBUGS will implement
inits = function() {
	list( beta0=0.0, beta=rep(0,p),
	     g=rep(0,p),.RNG.seed=321,.RNG.name = 'base::Wichmann-Hill') 
}




#########################

# The function  jags.model() compile the jags model, moreover it performs
#and adaptive burn-in for 'nchain' parallel chains.
#As input we have to provide: the model, the data, 
#the initial condition of the chains, the number of chains, and the number of
#adaptations (burn inn).

model=jags.model("SSVS_probit.bug",data=data_win,n.adapt=1000,inits=inits,n.chains=1) 








####Posterior parameters WinBUGS  has to save 
param <- c("beta0","beta","sigma2","g","mdl")


## Some other information sfor the 
## MCMC algorihm
#number of iterations
nit <- 50000
#thinning
thin <-10





##The command coda.samle() calls jags from R passing the data and initial
#value just defined

output=coda.samples(model=model,variable.names=param,n.iter=nit,thin=thin)




#save(output,file='ssvs_set3.Rdata')   # we save the chain

load('ssvs_set3.Rdata')

# the output is an mcmc object of the library coda
str(output)

### A qalitative analysis of the convergence of the posteror chais

###We give a look at the trace plot and at the density summary
### of the posterior chains for each of the parameters
plot(output,ask=T)


#We can close now the graphical device
dev.off()


### To work with the posterior chains it is better 
### to cast the output to be an array object of R

output <- as.matrix(output)

###Some variable selection thecniques:

# The median probability moodel (MPM)
# pick variables with estimated posterior inclusion probabilities 
# higher than 0.5
# Notice that the estimated posterior inclusion probabilities are the
# posterior means of the gamma variables (in the code we called g)

##We save the posterior chain of the inclusion variable in post.g
post.g <- output[,14:25]

# then we compute the sample mean, column by column
apply(post.g,2,"mean")
post_mean_g <- apply(post.g,2,"mean") 

## to produce a bar plot of the  posterior inclusion probabilities
#postscript("inclusion_set1.ps")
barplot(post_mean_g,names.arg=paste("x",1:p,sep=""),main="Posterior inclusion probabilities set 1")
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
paste("x",1:p,sep="")[mp_model]


## close all the graphical device
dev.off()


######
######  Highest posterior density model (HPD)

# pick a model with the highest estimated posterior probability 

# Recall that we have enumerated the model 
# using a binary coding
# mdl=2^g0+2^g1+...+2^gp
# the visited model are saved in the chain

plot(output[,"mdl"])

# for example at iteration 235 the chain explored the 
# model 
output[235,"mdl"]
#  What were the regressors whose corresponding  gamma's were
#  different from zero

# we will need a binary representation of the number   

# The following function has as input an integer and a  base p.
# The result is a vector containing the representation in
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


### Then we apply the function as.binary to output[235,"mdl"]

as.binary(output[235,"mdl"])

## It is better to read the output of as.binary in rverse order
rev(as.binary(output[235,"mdl"]))


##In the following function I will use the as.binary function
## and  a paste() trick two time to write down automatically, and in
## in a nice form, which variable are included in a visited models.

##The input of the funcrion is the model-number n and the number of covariate n.cov

which.variable <- function(n,n.cov){
	bin <- rev(as.binary(n))
	n.bin <- length(bin)
	logic.bin <- rep(FALSE,n.cov+1)
	for(i in 1:n.bin){
		logic.bin[i]=(bin[i]==1)
	}
	out <- paste(paste("x",0:p,sep="")[logic.bin],collapse="_")
	return(out)

}


which.variable(output[235,"mdl"] ,12)




### We start to analize how many models have been visited 
## by the posterior chain:
length(unique( output[,"mdl"]))

## Now we compute the posterior frequency of the visited  models
visited_models<-table(output[,"mdl"])
visited_models




#We can visualize the table of the visited models
x11()
barplot(visited_models)


## Let's sort the table to see which are the "top ten"  
top <- sort(visited_models,decreasing=T)[1:10]
top


## In the following lines we will use the which.variable function 
## to visualize which  are the variables included in the most visited models
numeric_top <- as.numeric(names(top))

for(j in 1:10){
names(top)[j]=which.variable(numeric_top[j],p)
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
colnames(output)
post.beta <- output[,1:12]

# then we compute the sample mean , column by column
mean.beta.post <- apply(post.beta,2,"mean")
sd.beta.post <- apply(post.beta,2,"sd")

require(plotrix)

#postscript("HS_decision_set1.ps")
plotCI(x=1:p, y=mean.beta.post, uiw=sd.beta.post,lwd=1.5, main="Decision intervals for HS, set1")
abline(h=0,col="blue")
#dev.off()
##




####################


###Boxplot of the posterior parameter
boxplot(as.data.frame(post.beta))
abline(h=0,col="blue",lwd=3,lty=2)



