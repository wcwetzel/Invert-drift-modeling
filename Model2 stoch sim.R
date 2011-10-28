#### simulate MODEL 2

# MODEL 2


# parameters:
N0 = 113 # initial number of herbivores in a channel
f = 1 # initial food abundance
m = 0.1 # coefficient
b = 0.01


## now stochastic simulation
library("GillespieSSA") # package for simulating stochastic 
	# continuous time processes, ODEs. with discrete state variables
	# function ssa is the workhorse


# define parameters and model for GillespieSSA
# (this takes parms from above)
parms = c(m = m, f = f, b=b) # parameters in dN/dt
x0 = c(N=N0) # initial N
nu = matrix(-1) # stage change vector, change in N caused by 1 event
pf = c("m*N - b * f * N")# propensity function, Pr of 1 event in infintesimal time 
	# interval [t, t + dt)
tf = 100 # time at which to end simulation
method = "D" # direct method: slow, stable
simName = "model2"

# Now simulate data for a range of initial food densities
# Run sim 30 times and save N at t = 10

flevels = seq(1, 9.5, by=0.5) # experimental levels of food
nperf = 30 # number of replicates per food
f = sort(rep(flevels, nperf))
data.sim2 = data.frame(rep=NA, N=NA, f=NA)
for(i in 1:length(f)){
	parms = c(m = m, f = f[i])
	out = ssa(x0, pf, nu, parms, tf, method, simName, verbose=FALSE)
	index = which(abs(out$data[,1]-10)==min(abs(out$data[,1]-10)))
	data.sim2[i,] = c(i, out$data[,2][index], out$args$parms['f'])
}


 plot(data.sim2[,'N'] ~ data.sim2[,'f'], ylab='Herbivores remaining', 
	 xlab='Food index')

#plot(I(113 - data.sim[,'N']) ~ data.sim[,'f'], ylab='Herbivores drifting', 
	#xlab='Food index')
