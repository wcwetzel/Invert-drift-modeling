#### simulate MODEL 3

# MODEL 3


# parameters:
N0 = 113 # initial number of herbivores in a channel
f = 1 # initial food abundance
c = 0.01
m = 0.1 # coefficient

## now stochastic simulation
library("GillespieSSA") # package for simulating stochastic 
	# continuous time processes, ODEs. with discrete state variables
	# function ssa is the workhorse


# define parameters and model for GillespieSSA
# (this takes parms from above)
parms = c(m = m, c=c, f = f) # parameters in dN/dt
x0 = c(N=N0) # initial N
nu = matrix(-1) # stage change vector, change in N caused by 1 event
pf = c("m*N / (c + f)")# propensity function, Pr of 1 event in infintesimal time 
	# interval [t, t + dt)
tf = 100 # time at which to end simulation
method = "D" # direct method: slow, stable
simName = "drift"

# Now simulate data for a range of initial food densities
# Run sim 30 times and save N at t = 10

f = sort(rep(1:10, 30))
data.sim3 = data.frame(rep=NA, N=NA, f=NA)
for(i in 1:300){
	parms = c(m = m, c = c, f = f[i])
	out = ssa(x0, pf, nu, parms, tf, method, simName, verbose=FALSE)
	index = which(abs(out$data[,1]-10)==min(abs(out$data[,1]-10)))
	data.sim3[i,] = c(i, out$data[,2][index], out$args$parms['f'])
}


 plot(data.sim3[,'N'] ~ data.sim3[,'f'], ylab='Herbivores remaining', 
	 xlab='Food index')

# plot(I(113 - data.sim[,'N']) ~ data.sim[,'f'], ylab='Herbivores drifting', 
	# xlab='Food index')
