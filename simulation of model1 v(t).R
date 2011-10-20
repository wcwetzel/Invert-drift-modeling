# This is to check that my solutions to model1's ODEs are correct
# 10 Oct 2011

# E(N(t)) is expected number of individuals left in the stream
# V(t) is the variance in the number left

# ODEs for dE(N)/dt and dV/dt:
# dE(N)/dt = f(E(N))
# dV/dt = -2f'(N)V + f(N)
# from Carl Boettiger

# model1 and my analytical solutions:
# dN/dt = -mN
# N(t) = N0*exp(-mt)
# dV/dt = 2mV -mN
# V(t) = exp(2mt) * (N0*exp(-3mt)/3 + c)

# initial number and parameters:
N0 = 100
V0 = N0/3 # not sure what V0 should be!
m = 0.2 # emigration rate
t = seq(0,10, by=0.1) # timesteps

# analytical solutions for N(t) and V(t):
Nt1 = N0 * exp(-m * t)
Vt1 = exp(2*m*t) * (N0*exp(-3*m*t)/3 + (V0 - N0/3))

# Numerically simulate ODEs to test analytical solutions:
require(deSolve)

model1 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N					# dN/dt
		dv.dt = 2 * m * V - m * N		# dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}

parms = c(N0 = N0, V0 = V0, m = m)
start = c(N = N0, V = V0)
# run the model:
out.model1 = as.data.frame(lsoda(start, times=t, model1, parms))


# graphically compare numerical and analytical results
# solid line = analytical results
# dashes = numerical results
par(mfrow=c(1,2))
plot(Vt1 ~ t, type='l', 
	ylim=c(min(c(Vt1, out.model1$V)), max(c(Vt1, out.model1$V))))
points(out.model1$time, out.model1$V, type='l', lty=3, lwd=4)
plot(Nt1 ~ t, type='l')
points(out.model1$time, out.model1$N, type='l', lty=3, lwd=4)




