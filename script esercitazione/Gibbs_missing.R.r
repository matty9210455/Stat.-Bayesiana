#########################################
#### Gibbs sampler for missing data  ###
#########################################
rm(list=ls())
# Data
y <- matrix(c(-1,-1,1,1,2,-2,2,-2,NA,NA,NA,NA,-1,1,-1,1,NA,NA,NA,NA,2,-2,2,-1), 
            12, 2)
dim(y)
y
n <- nrow(y)
n
colnames(y) <- c("y1", "y2")
plot(y, pch = 19)

# vector of the augmented observations
y.aug <- y

# initial value for Sigma
Sig <- diag(2)
# mu is the mean, set equal to 0.

# Indexes related to the missing data
miss <- which(is.na(y), arr.ind = T)
miss

# arr.ind: should array indices be returned when x is an array?
righe <- miss[,1] # rows
colonne <- miss[,2] # columns
righe   
colonne 
# Number of missing data
n.mis <- length(righe)
n.mis

## Prior distribution for matrix Sigma
n0 <- 4   # nu
Psi <- 2* diag(2)
Psi

# Number of iterations
iteff <- 5000
burnin <- 10000
thin <- 10
niter <- iteff*thin + burnin
niter

# Parameter: correlation between the components of the vector
rho <- vector(length = iteff)
missing.new <- NULL

set.seed(313)
# Gibbs sampler
for(iter in 1:niter)
{
  ## Step 1: simulate the missing data
	for(i in 1:n.mis)
  {
		if(colonne[i] == 1)
    {
		  ## We simulate the Observations
		  ## related to Missing information on the First Component:
		  ## We need to simulate [Y_1 | Y_2]
		  ## They will be Gaussian with mean and variance,
      ## as we saw at the blackboard
			mu1 <- Sig[1,2]/Sig[2,2]*y[righe[i],2]
			s2_1 <- Sig[1,1] - Sig[1,2]^2/Sig[2,2] # varianza
			y.aug[righe[i],1] <- rnorm(1, mean = mu1, sd = sqrt(s2_1))
		} 
		else
    {
      ## We simulate the Observations
      ## related to Missing information on the Second Component:
      ## We need to simulate [Y_2 | Y_1]
      ## They will be Gaussian with mean and variance,
      ## as we saw at the blackboard
      mu2 <- Sig[2,1]/Sig[1,1]*y[righe[i],1]
			s2_2 <- Sig[2,2] - Sig[2,1]^2/Sig[1,1]
			y.aug[righe[i],2] <- rnorm(1, mean = mu2, sd = sqrt(s2_2))
		}
	}

  ## Step 2: simulate matrix Sigma, conditionally on 
  ## the vector of data augmented with the missing values
  # Compute S:
	S <- t(y.aug)%*%y.aug
	
	# print(S)
	# Parameters of the function to draw from the Wishart distribution:
	# n 	 integer sample size.
	# df 	 numeric parameter, “degrees of freedom”.
	# Sigma 	 positive definite (p * p) “scale” matrix, the matrix parameter 
	#        of the distribution.
 
	## The function rWishart is in the library of R
	Sig <- rWishart(1, df = n + n0, Sigma = solve(S+Psi))[, ,1]
	## Invert the matrix to obtain a sample from an invWishart
	Sig <- solve(Sig)

	if(iter>burnin & (iter-burnin)%%thin==0)
	{ ## Compute the correlation parameter at each iteration
	  rho[(iter-burnin)/thin] <- (Sig[1,2]/sqrt(Sig[1,1]*Sig[2,2]))
	  missing.new <- c(missing.new, y.aug[5,2])
	}
	if(iter%%10000 == 0) cat("*** Iteration number ", iter,"/", niter,"\n")
}
# End Gibbs

x11()
par(mfrow = c(1,2))
hist(rho, nclass="fd", freq = F)
lines(density(rho), col = "darkblue", lwd = 2)
plot(rho, type = 'l')
# The posterior of rho has two peaks: it is as if in our data there were two "groups".
# The data for which the correlation between the first component and 
# the second is positive,
# And those for which the correlation is negative.


# Let's have a look at the missing observations that we sampled!
head(missing.new)
x11()
par(mfrow = c(1,2))
hist(missing.new, prob = T, main = "Missing value: y_(5,2)")
lines(density(missing.new), col = "darkblue", lwd = 2)
plot(missing.new, type = 'l')

# Check convergence with CODA...

###############################################################################################
## Try with the Metropolis Hastings!
library(mvtnorm)
library(bayesm)

## Log posterior up to a normalizing constant
lposterior <- function(y, Sig, n0, Psi)
{
	out = 0
	for(i in 1:4)
  {
	  ## The first 4 observations have no missing data
	  ## So they Fully contribute to the likelihood
	  ## (bivariate normal)
		out = out + dmvnorm(y[i,], mean=c(0,0), sig=Sig, log=T)
	}
	for( i in 5:8)
  {
	  # These observations have a missing component
	  # (The second) --> the contribution to the likelihood
	  # Is partial and given by a univariate normal (MARGINAL)
		out = out + dnorm(y[i,1], mean = 0, sd = sqrt(Sig[1,1]), log = T)
	}
	for(i in 9:12)
  {
	  # These observations have a missing component
	  # (The first) --> the contribution to the likelihood
	  # Is partial and given by a univariate normal (MARGINAL)
	  out = out + dnorm(y[i,2], mean = 0, sd = sqrt(Sig[2,2]), log = T)
	}
	# Prior information: 
	# use the package "bayesm" that contains the function 
	# lndIWishart() to compute the log of the inv-Wishart density
	out = out + lndIWishart(nu = n0, V = Psi, IW = Sig)
	if(iter%%1000 == 0) cat("*** Iteration number ", iter,"/", niter,"\n")
	return(out)
}

#### Metropolis Hastings ####
## Starting point
Sig <- diag(1,2)
# Posterior chain of \rho, correlation
rho_MH <- vector()

## Parameters of the proposal

## The proposal is given by the prior (!it is not an optimal choice!)
Psiinv <- solve(Psi)
np = 2

niter = 5000 

system.time( # --> time
for(iter in 1:niter)
{
	# propose a new value for the matrix
	Delta <-  rWishart(1, df = np, Sigma = Psiinv)[,,1]
	Delta <- solve(Delta)

  # log of the acceptance rate 
  log_acp <- lposterior(y, Delta, n0, Psi) + lndIWishart(nu = np, V = Psi, IW = Sig) # numeratore
	log_acp <- log_acp - lposterior(y, Sig, n0, Psi) - lndIWishart(nu = np, V = Psi, IW = Delta)  # denominatore
	lu <- log(runif(1))
	
	## Accept or reject?
	if(lu < log_acp)
	{
		Sig = Delta
	}

	## compute the correlation at each iteration
	rho_MH[iter] <- (Sig[1,2]/sqrt(Sig[1,1]*Sig[2,2]))
}
)

x11()
par(mfrow = c(1,2))
hist(rho_MH, nclass="fd", freq = F)
lines(density(rho_MH), col = "red", lwd = 2)
lines(density(rho), col = "darkgreen", lwd = 2)

plot(rho_MH, type = 'l')
# save(rho_MH,file="rho_H.Rdat")

graphics.off()
## END