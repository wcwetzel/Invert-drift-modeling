#### simulate drift data
## 10 sep 2011

# drift is a linear function of N(t) and rational function of initial food f.
# dN/dt = f(N) = (b/f) * N
# N(t) = N0 * exp( (b/f) * t)

# parameters:
N0 = 113 # initial number of herbivores in a channel
f = 1 # initial food abundance
b = -0.1 # coefficient

# plot deterministic drift model through time
curve(N0 * exp(b * x / f), 0 , 55) 

## now stochastic simulation
library("GillespieSSA") # package for simulating stochastic 
	# continuous time processes, ODEs. with discrete state variables
	# function ssa is the workhorse

# dN/dt = F(N,f) = emigration rate
# Pr(N -> N-1) = h(N, f) = individual transition probability
# N * h(N, f) = F(N, f)
# dN/dt = N * h(N, f)
# h(N, f) = b * N / f simple, linear starting point
# dN/dt = b * N^2 / f

# define parameters and model for GillespieSSA
# (this takes parms from above)
parms = c(b = abs(b), f = f) # parameters in dN/dt, b is |b| b/c nu is -1
x0 = c(N=N0) # initial N
nu = matrix(-1) # stage change vector, change in N caused by 1 event
a = c("b*N/f")# propensity function, Pr of 1 event in infintesimal time 
	# interval [t, t + dt)
tf = 100 # time at which to end simulation
method = "D" # direct method: slow, stable
simName = "linear drift"

# now finally run the sim
out = ssa(x0, a, nu, parms, tf, method, simName, verbose=TRUE, 
	consoleInterval = 1)

plot(out$data[,2] ~ out$data[,1]) # plot sim results
curve(N0 * exp(b * x / f), 0 , 100, add=TRUE, col='red') # add
	# deterministic curve


# Now simulate data
# Run sim 30 times and save N at t = 10

f = sort(rep(1:10, 30))
data.sim = data.frame(rep=NA, N=NA, f=NA)
for(i in 1:300){
	parms = c(b = abs(b), f = f[i])
	out = ssa(x0, a, nu, parms, tf, method, simName, verbose=FALSE)
	index = which(abs(out$data[,1]-10)==min(abs(out$data[,1]-10)))
	data.sim[i,] = c(i, out$data[,2][index], out$args$parms['f'])
}

hist(data.sim[,'N'])
mean(data.sim[,'N'])
var(data.sim[,'N'])

plot(data.sim[,'N'] ~ data.sim[,'f'], ylab='Herbivores remaining', 
	xlab='Food index')

plot(I(113 - data.sim[,'N']) ~ data.sim[,'f'], ylab='Herbivores remaining', 
	xlab='Food index')
