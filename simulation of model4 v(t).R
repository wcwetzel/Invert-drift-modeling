# This is to check that my solutions to model4's ODEs are correct
# 10 Oct 2011

# E(N(t)) is expected number of individuals left in the stream
# V(t) is the variance in the number left

# ODEs for dE(N)/dt and dV/dt:
# dE(N)/dt = f(E(N))
# dV/dt = -2f'(N)V + f(N)
# from Carl Boettiger

# model1 and my analytical solutions:
# dN/dt = -mN - aN^2
# N(t) = m / ((m/N0 + a) * exp(t*m) - a)
# dV/dt = -2(-m - 2aN)V - mN - aN^2
# V(t) = N(t)(m + aN(t)) / (2(m + 2aN(t))) + 
#		(V0 - N0(m + aN0) / (2(m + 2aN0))) exp(2t(m + 2aN(t)))

# initial number and parameters:
N0 = 100
V0 = 100 # not sure what V0 should be!
m = 0.2 # emigration rate
a = 0.001
t = seq(0,10, by=0.1) # timesteps

# analytical solutions for N(t) and V(t):
Nt4 = m / ((m/N0 + a) * exp(t*m) - a)
Vt4 = Nt4 * (m + a*Nt4) / (2*(m + 2*a*Nt4)) + 
		(V0 - N0*(m + a*N0) / (2*(m + 2*a*N0))) * exp(2*t*(m + 2*a*Nt4))

# also solve model1 for comparison
Nt1 = N0 * exp(-m * t)
Vt1 = Nt1 / 2 + (V0 - N0/2) * exp(2 * m * t)


# Numerically simulate ODEs to test analytical solutions:
require(deSolve)

model1 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N - a * N^2						# dN/dt
		dv.dt = (2*m + 4*a*N) * V - m*N - a*N^2		# dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}

parms = c(N0 = N0, V0 = V0, m = m, a=a)
start = c(N = N0, V = V0)
# run the model:
out.model4 = as.data.frame(lsoda(start, times=t, model1, parms))


# graphically compare numerical and analytical results
# solid line = analytical results
# dashes = numerical results
# green dashes = model1 analytical results
par(mfrow=c(1,2))
# plot variance
plot(Vt4 ~ t, type='l', 
	ylim=c(min(c(Vt4, out.model4$V)), max(c(Vt4, out.model4$V))))
points(out.model4$time, out.model4$V, type='l', lty=3, lwd=4)
points(Vt1 ~ t, type='l', lty=2, lwd=1, col='green')

# plot E(N)
plot(Nt4 ~ t, type='l')
points(out.model4$time, out.model4$N, type='l', lty=3, lwd=4)
points(Nt1 ~ t, type='l', lty=2, lwd=1, col='green')
# green should overlay black when a = 0


