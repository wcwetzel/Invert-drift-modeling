## likelihood for simple linear drift simulation
# assumes food density doesn't change
# Pr(drift) depends on initial food and current density
library(bbmle)
# 1st run 'simple drift simulation.R' to generate data


E.Nt = function(N0, b, f, t){
	N0 * exp(b * t / f)
}
E.Nt(113, -0.1, 1, 10)


Vt = function(N0, b, f, t, V0){
	E.Nt(N0, b, f, t) / 2 + (V0 - N0/2) * exp(-2 * b * t / f)
}
SDt = function(N0, b, f, t, V0){
	sqrt( E.Nt(N0, b, f, t) / 2 + (V0 - N0/2) * exp(-2 * b * t / f) )
}
# WHAT IS V0????? it gives a negative variance when = 0
# simulated variance is about 22 at t=10 with b=-0.1
# a V0 of about 56.9 will give approx this variance

SDt(113, -0.1, 1, 10, 56.9)

nLL = function(N, b, f, N0, V0, t){
	-sum( dnorm(x = N, 
		mean = E.Nt(N0, b, f, t), 
		sd = SDt(N0, b, f, t, V0),
		log = TRUE) 
	)
}

nLL.simp = function(N, m, s){
	-sum( dnorm(x = N, 
		mean = m, 
		sd = s,
		log = TRUE) 
	)
}


model = mle2(minuslogl = nLL, 
	start = list(b = -0.11, V0 = 56.9),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']))

model.simp = mle2(minuslogl = nLL.simp,
	start = list(m = mean(data.sim[,'N']), 
		s=sd(data.sim[,'N'])),
	data = list(N = data.sim[,'N']))


AICtab(model, model.simp)

