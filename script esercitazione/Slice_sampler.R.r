###########################
### Slice sampler #########
###########################
rm(list=ls())
# Consider the following density: 
# d/Gamma(1/d) *exp(-x^(d))

# Define a function to compute the density
dens <- function(x,d)
{
  return(d/gamma(1/d)*exp(-x^(d)))
}

# Define a function to compute the density in logarithmic scale
ldens <- function(x,d)
{
  return(log(d) - lgamma(1/d) - x^d)
}

# Draw the density for d = 1/2
x11()
curve(dens(x, 1/2), from = 0, to = 30, add = F, col = "red", n = 500, lwd = 2)

### Slice sampler to sample from f()

## As Input parameters we have:
## - The number of iterations "niter"
## - The initial value of the Markov Chain "x0"
## - The parameter "d" of the density
slice.dens <- function(niter, x0 = 1, d = 1/2)
{
  # vector containig the auxiliary r.v. Unif(0,f(x))
  u.out <- vector(length = niter)
  # vector containing the r.v. of interest 
  x.out <- vector(length = niter) 
  for(j in 1:niter)
  {
    # sample U ~ Unif(0,f(x))   
    u.out[j] <- runif(1, min = 0, max = dens(x0, d))
    z = log(u.out[j]) # log scale
    # In this example, the set S
    # has the form (0, app) where app is
    app <- (log(d) - lgamma(1/d) - z)^(1/d)  

    ## A new value fo the chain is sampled
    ## from a Unif (0, app)    
    x0 <- runif(n = 1, min = 0, max = app)
    x.out[j] <- x0
  }	
  return(list(x = x.out, u = u.out))
}


niter <- 6000
d <- 1/2
## d <- 1/20 
set.seed(1905)
chain <- slice.dens(niter, x0 = 1, d = d)

x11()
hist((chain$x), nclass="fd", freq = F, col = "gray", main = "Hist of the chain", xlim = c(0,60))
curve(dens(x,d), n = 500, from = 0, add = T, col = "red", lwd = 1.6)

# Traceplot
plot(chain$x, type = "l") 
# it has a very heavy right tail: this is due to the distribution
# From which we are sampling

# Acf plot
par(mfrow=c(2,1))
acf(chain$x)
acf(chain$u)
# Note: the chain are quite correlated... 
# a solution would be the Thinning

plot(cumsum(chain$x)/(1:length(chain$x)), type = 'l')
# Note: you can see that at first the chain did not converged. 
# We throw away the initial iterations
# that we do not want to use for estimation! BURNIN
burnin = 1000
thin = 5
niter
chain.mat = cbind(chain$x[ seq(burnin + 1, niter, by = thin)],chain$u[ seq(burnin + 1, niter, by = thin)] )
dim(chain.mat)

require(coda)
chain.mcmc = mcmc(chain.mat, start = burnin +1, end = niter, thin = thin)

plot(chain.mcmc)
cumuplot(chain.mcmc)

par(mfrow = c(1,2))
acf(chain.mat[,1])  
acf(chain.mat[,2])  

graphics.off()
#### Graphic visualization of the Slice sampler
# draw the density
x11()
curve(dens(x,d), n = 500, from = 0, to = 5, add = F, col = "blue", lwd = 2, ylim = c(0,0.5))
# The asterisk represents the value x0 
points(1,0, pch="*", cex = 2, col="blue")
# segment (0,f(x0))
segments(x0 = 1, x1 = 1, y0 = 0, y1 = dens(1,d), lty=2)
# first value for  U ~ Unif(0,f(x0))
points(1, chain$u[1], pch = 4, cex = 1)
segments(x0 = 1, x1 = 1, y0 = 0, y1 = chain$u[1], col = "red", lwd = 1.5)

# Set S 
segments(x0 = (-log(chain$u[1]) + log(d)- log(gamma(1/d)))^(1/d), x1 = 0, y0 = chain$u[1], y1 = chain$u[1])
##  x1
points(chain$x[1],chain$u[1],lty=2)
segments(x0 = 1, x1 = chain$x[1], y0 = chain$u[1], y1 = chain$u[1], col = "red", lwd = 2)
points(chain$x[1],0,pch="*",cex=2,col="blue") 
segments(x0 = chain$x[1], x1 = chain$x[1], y0 = 0, y1 = dens(chain$x[1],d), lty = 2)


points(chain$x[1] ,chain$u[2],pch=4,cex=1)
segments(x0=chain$x[1],x1=chain$x[1],y0=chain$u[1],y1=chain$u[2],col="red",lwd=2)
segments(x0=(-log(chain$u[2])+log(d)-log(gamma(1/d)))^(1/d),x1=0,y0=chain$u[2],y1=chain$u[2],lty=2)
points(chain$x[2],chain$u[2])
segments(x0=chain$x[1],x1=chain$x[2],y0=chain$u[2],y1=chain$u[2],col="red",lwd=2)
points(chain$x[2],0,pch="*",cex=2,col="blue") 

segments(x0=chain$x[2],x1=chain$x[2],y0=0,y1=dens(chain$x[2],d),lty=2)
points(chain$x[2] ,chain$u[3],pch=4,cex=1)
segments(x0=chain$x[2],x1=chain$x[2],y0=chain$u[2],y1=chain$u[3],col="red",lwd=2)
segments(x0=(-log(chain$u[3])+log(d)-log(gamma(1/d)))^(1/d),x1=0,y0=chain$u[3],y1=chain$u[3],lty=2)
points(chain$x[3],chain$u[3])
segments(x0=chain$x[2],x1=chain$x[3],y0=chain$u[3],y1=chain$u[3],col="red",lwd=2)
points(chain$x[3],0,pch="*",cex=2, col="blue")

for(i in 1:30)
{
  points(chain$x[2+i] ,chain$u[3+i],pch=4,cex=1)
  segments(x0=chain$x[2+i],x1=chain$x[2+i],y0=chain$u[2+i],y1=chain$u[3+i],col="red",lwd=2)
  points(chain$x[3+i],chain$u[3+i])
  segments(x0=chain$x[2+i],x1=chain$x[3+i],y0=chain$u[3+i],y1=chain$u[3+i],col="red",lwd=2)
  points(chain$x[3+i],0,pch="*",cex=2,col="blue")
  Sys.sleep(1.)
}

# The chain must explore thoroughly all the space below the graph!
points(chain$x, chain$u, cex = 0.3, col = "yellow")


################ End ################################