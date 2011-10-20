# This is to check that my solutions to the ODEs for models 1,2,3 are correct
# 19 Oct 2011

# E(N(t)) is expected number of individuals left in the stream
# V(t) is the variance in the number left

# ODEs for dE(N)/dt and dV/dt:
# dE(N)/dt = f(E(N))
# dV/dt = -2f'(N)V + f(N)
# from Carl Boettiger

# model1 and my analytical solutions:
# dN/dt = -mN
# N(t) = N0*exp(-mt)
# dV/dt = 2mV(t) -mN(t)
# V(t) = exp(2mt) * (V0 + N0 * (exp(-3 * m * t) - 1) / 3)

# model2 and my analytical solutions:
# dN/dt = -mN + bFN
# N(t) = N0*exp(t(bF-m))
# dV/dt = -2(bF - m)V(t) - mN(t) + bfN(t)
# V(t) = exp(2t(m-bF)) * (V0 + N0(exp(3t(bf-m))-1)/3)

# model3 and my analytical solutions:
# dN/dt = -mN / (c + F)
# N(t) = N0*exp(-mt / (c + F))
# dV/dt = 2mV(t) / (c + F) - mN(t) / (c + F)
# V(t) = exp(2mt / (c + F)) * (V0 + N0 * (exp(-3*m*t/(c +F)) - 1)/3)






# initial number and parameters:
N0 = 100
V0 = N0/3 # not sure what V0 should be!
m = 0.2 # emigration rate
t = seq(0,10, by=0.1) # timesteps
f = 1
b = 0.1
c = 0.1

# analytical solutions for N(t) and V(t) for all 3 models:
Nt1 = N0 * exp(-m * t)
Vt1 = exp(2*m*t) * (N0*exp(-3*m*t)/3 + (V0 - N0/3))
Nt2 = N0 * exp(t*(b*f - m))
Vt2 = exp(2*t*(m-b*f)) * (V0 + N0 * (exp(3 * t * (b*f - m)) - 1)/3)
Nt3 = N0 * exp(-m * t / (c + f))
Vt3 = exp(2 * m * t / (c + f)) * (V0 + N0 * (exp(-3*m*t/(c + f)) - 1)/3)

# plot analytical results
par(mfrow=c(1,2))
plot(Nt1 ~ t, type='l', lty=2) # model1 = black
points(Nt2 ~ t, type='l', lty=3, col=2) # model2 = red
points(Nt3 ~ t, type='l', lty=4, col=4) # model3 = blue

plot(Vt1 ~ t, type='l', lty=2) # model1 = black
points(Vt2 ~ t, type='l', lty=3, col=2) # model2 = red
points(Vt3 ~ t, type='l', lty=4, col=4) # model3 = blue


# Numerically simulate ODEs to test analytical solutions:
require(deSolve)

parms = c(N0 = N0, V0 = V0, m = m, f=f, c=c, b=b)
start = c(N = N0, V = V0)

model1 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N					# dN/dt
		dv.dt = 2 * m * V - m * N		# dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}
out.model1 = as.data.frame(lsoda(start, times=t, model1, parms))

model2 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N	+ b*f*N				# dN/dt
		dv.dt = -2*(b*f - m)*V - m*N + b*f*N		# dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}
out.model2 = as.data.frame(lsoda(start, times=t, model2, parms))

model3 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N / (c + f)					# dN/dt
		dv.dt = 2*m*V / (c + f) - m*N / (c + f)	# dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}
out.model3 = as.data.frame(lsoda(start, times=t, model3, parms))


# graphically compare numerical and analytical results
par(mfrow=c(3,2))
plot(Nt1 ~ t, type='l', lty=2, col='blue')
points(out.model1$N ~ t, type='l', lty=4, col='red')

plot(Vt1 ~ t, type='l', lty=2, col='blue')
points(out.model1$V ~ t, type='l', lty=4, col='red')

plot(Nt2 ~ t, type='l', lty=2, col='blue')
points(out.model2$N ~ t, type='l', lty=4, col='red')

plot(Vt2 ~ t, type='l', lty=2, col='blue')
points(out.model2$V ~ t, type='l', lty=4, col='red')

plot(Nt3 ~ t, type='l', lty=2, col='blue')
points(out.model3$N ~ t, type='l', lty=4, col='red')

plot(Vt3 ~ t, type='l', lty=2, col='blue')
points(out.model3$V ~ t, type='l', lty=4, col='red')
