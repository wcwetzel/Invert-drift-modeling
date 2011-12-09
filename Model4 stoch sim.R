#### simulate MODEL 4

# MODEL 4


# parameters:
N0 = 113 # initial number of herbivores in a channel
f = 1 # initial food abundance
m = 0.01 # coefficient of drift
a = 0.001 # coef for density effect on drift

## now stochastic simulation
library("GillespieSSA") # package for simulating stochastic 
	# continuous time processes, ODEs. with discrete state variables
	# function ssa is the workhorse


# define parameters and model for GillespieSSA
# (this takes parms from above)
parms = c(m = m, a=a, f=f) # parameters in dN/dt
x0 = c(N=N0) # initial N
nu = matrix(-1) # stage change vector, change in N caused by 1 event
pf = c("m*N + a*N^2")# propensity function, Pr of 1 event in infintesimal time 
	# interval [t, t + dt)
tf = 100 # time at which to end simulation
method = "D" # direct method: slow, stable
simName = "model4"

# Now simulate data for a range of initial food densities
# Run sim 10 time steps and save N at t = 10

#nperfN = 1 # number of replicates per food/N combo
fN = expand.grid(f = seq(1, 5, by=1), N = seq(113, 213, by=20))
fN = rbind(fN, fN, fN, fN)

data.sim4 = data.frame(rep=NA, N=NA, f=NA, N0=NA)
for(i in 1:nrow(fN)){
	x0 = c(N = fN[i,2])
	parms = c(m = m, a = a, f = fN[i,1])
	out = ssa(x0, pf, nu, parms, tf, method, simName, verbose=FALSE)
	index = which(abs(out$data[,1]-10)==min(abs(out$data[,1]-10)))[1]
	data.sim4[i,] = c(i, out$data[,2][index], out$args$parms['f'], x0)
	}


 plot(data.sim4[,'N'] ~ data.sim4[,'f'], ylab='Herbivores remaining', 
	 xlab='Food index')

# plot(I(113 - data.sim[,'N']) ~ data.sim[,'f'], ylab='Herbivores drifting', 
	# xlab='Food index')
