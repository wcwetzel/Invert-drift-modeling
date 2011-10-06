# This is to check that my solution to model4's ODE is correct
# 4 Oct 2011

# model4: dN/dt = -mN - aN^2
# model1: dN/dt = -mN
# note that model4 reduces to model1 when a=0

N0 = 113
m = 0.5
a = 0.001
t = seq(0,20, by=0.1)

Nt4 = m /( (m/N0 + a) * exp(t*m) - a)


plot(Nt4 ~ t, type='l', col='blue')

Nt1 = N0 * exp(-m * t)

points(Nt1 ~ t, type='l', col='red', lty=2)

require(deSolve)

model4 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N - a * N^2
		res = dn.dt
		list(res)
	})
}

parms = c(N0 = N0, m = m, a = a)

times = t

start = c(N = N0)

out.model4 = as.data.frame(lsoda(start, times, model4, parms))

points(out.model4$time, out.model4$N, type='l', col='blue', lty=3, lwd=4)

# IT WORKS!!
# red line is model 1
# blue line is model4 using my analytical solution of ODE
# thick blue dashes is model4 solved numerically
# model4 should equal model1 when a=0