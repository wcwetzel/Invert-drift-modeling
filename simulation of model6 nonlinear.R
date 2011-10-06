# This is to check that my solution to model4's ODE is correct
# 4 Oct 2011

# model6: dN/dt = (-mN - aN^2) / (c + F)
# model1: dN/dt = -mN
# note that model6 reduces to model1 when c=0 and F=1 or c=1 and F=0

N0 = 113
f = 1
m = 0.1
a = 0
c = 1
t = seq(0,20, by=0.1)

Nt6 = m /( (m/N0 + a) * exp(t * m / (c + f)) - a)


plot(Nt6 ~ t, type='l', col='blue', ylim=c(0,120))

Nt1 = N0 * exp(-m * t)

points(Nt1 ~ t, type='l', col='red', lty=2)

require(deSolve)

model6 = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = (-m*N - a * N^2) / (c + f)
		res = dn.dt
		list(res)
	})
}

parms = c(N0 = N0, m = m, a = a, c = c, f = f)

times = t

start = c(N = N0)

out.model6 = as.data.frame(lsoda(start, times, model6, parms))

points(out.model6$time, out.model6$N, type='l', col='blue', lty=3, lwd=4)

# It works for some parameters....
# red line is model 1
# blue line is model5 using my analytical solution of ODE
# thick blue dashes is model5 solved numerically
# model4 should equal model1 when a=0