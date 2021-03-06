#### simulate MODEL 1
## 10 sep 2011

# drift is a linear function of N(t) and rational function of initial food f.
# MODEL 1
# dN/dt = f(N) = -m* N 
# N(t) = N0 * exp( -m * t)

# parameters:
N0 = 113 # initial number of herbivores in a channel
f = 1 # initial food abundance
m = 0.1 # coefficient

## now stochastic simulation
library("GillespieSSA") # package for simulating stochastic 
	# continuous time processes, ODEs. with discrete state variables
	# function ssa is the workhorse


# define parameters and model for GillespieSSA
# (this takes parms from above)
parms = c(m = m, f = f) # parameters in dN/dt
x0 = c(N=N0) # initial N
nu = matrix(-1) # stage change vector, change in N caused by 1 event
pf = c("m*N")# propensity function, Pr of 1 event in infintesimal time 
	# interval [t, t + dt)
tf = 100 # time at which to end simulation
method = "D" # direct method: slow, stable
simName = "model1"

# Now simulate data for a range of initial food densities
# Run sim 30 times and save N at t = 10

flevels = seq(1, 9.5, by=0.5) # experimental levels of food
nperf = 30 # number of replicates per food
f = sort(rep(flevels, nperf))
data.sim1 = data.frame(rep=NA, N=NA, f=NA)
for(i in 1:length(f)){
	parms = c(m = m, f = f[i])
	out = ssa(x0, pf, nu, parms, tf, method, simName, verbose=FALSE)
	index = which(abs(out$data[,1]-10)==min(abs(out$data[,1]-10)))
	data.sim1[i,] = c(i, out$data[,2][index], out$args$parms['f'])
}


plot(data.sim1[,'N'] ~ data.sim1[,'f'], ylab='Herbivores remaining', 
	xlab='Food index')

# plot(I(113 - data.sim[,'N']) ~ data.sim[,'f'], ylab='Herbivores drifting', 
	# xlab='Food index')
