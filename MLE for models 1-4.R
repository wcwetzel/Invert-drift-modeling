## likelihood for simple linear drift simulation
# assumes food density doesn't change
# Pr(drift) depends on initial food and current density
library(bbmle)
library(deSolve)
# 1st run 'Model1 stoch sim.R', 'Model2 stoch sim.R',
# or 'Model3 stoch sim.R' to generate data
# or 'Model4 stoch sim.R'

# Analytical models and NLLs for models1-4
# V(t) for model 4 is numerical simulation


# Model1
E.Nt1 = function(N0, m, t){
	N0 * exp(-m * t)
}
E.Nt1(113, 0.1, 10)
E.SDt1 = function(N0, m, t, V0){
	sqrt( exp(2 * m * t) * (V0 + N0 * ({exp(-3 * m * t) - 1} / 3)) )
}
E.SDt1(113, 0.1, 10, 56.9)
nll1 = function(N, m, N0, V0, t){
	-sum( dnorm(x = N, 
		mean = E.Nt1(N0, m, t), 
		sd = E.SDt1(N0, m, t, V0),
		log = TRUE) 
	)
}

# Model2
E.Nt2 = function(N0, m, b, f, t){
	N0 * exp((b*f - m) * t)
}
E.Nt2(113, 0.1, 0.01, 1:5, 10)

E.SDt2 = function(N0, m, b, f, t, V0){
	sqrt( exp(2 * t * (m - b * f)) * (V0 + N0 * {exp(3 * t * (b * f - m)) - 1} / 3) )
}
E.SDt2(113, 0.1, 0.01, 1:5, 10, 56.9)
nll2 = function(N, m, b, f, N0, V0, t){
	-sum( dnorm(x = N, 
		mean = E.Nt2(N0, m, b, f, t), 
		sd = E.SDt2(N0, m, b, f, t, V0),
		log = TRUE) 
	)
}

# Model3
E.Nt3 = function(N0, m, c, f, t){
	N0 * exp(-m * t / (c + f))
}
E.Nt3(113, 0.1, 0.01, 1:5, 10)

E.SDt3 = function(N0, m, c, f, t, V0){
	sqrt( exp(2 * m * t / (c + f)) * (V0 + N0 * {exp(3 * m * t / (c + f)) - 1} / 3) )
}
E.SDt3(113, 0.1, 0.01, 1:5, 10, 56.9)
nll3 = function(N, m, c, f, N0, V0, t){
	-sum( dnorm(x = N, 
		mean = E.Nt3(N0, m, c, f, t), 
		sd = E.SDt3(N0, m, c, f, t, V0),
		log = TRUE) 
	)
}

# Model4
E.Nt4 = function(N0, m, a, f, t){
	m / ((m/N0 + a) * exp(t*m) - a)
}
E.Nt4(113:114, 0.1, 0.01, 1:5, 10)

model4.sim = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -m*N - a*N^2					# dN/dt
		dv.dt = -2 * (- m - 2 * a * N) * V + (- m*N - a * N^2) # dV/dt
		res = c(dn.dt, dv.dt)
		list(res)
	})
}
	
E.SDt4 = function(N0, m, a, f, t, V0){
	Vvec = vector()
	sdvec = vector()
	Nvec = vector()
	for(i in 1:length(N0)){
		N0i = N0[i]
		initial = c(N = N0i, V = V0)
		parms = c(m = m, f=f, a=a)
		out.model4 = as.data.frame(lsoda(initial, times=c(0,t), model4.sim, parms))
		V = out.model4$V[length(out.model4$V)]
		sd = sqrt(out.model4$V[length(out.model4$V)])
		Vvec[i] = V
		sdvec[i] = sd
		Nvec[i] = out.model4$N[length(out.model4$N)]
	}
	return(data.frame(f=f, N0=N0, N = Nvec, V = Vvec, sd = sdvec))
}
E.SDt4(113:200, 0.15, 0, 1, 10, 500/3)

E.SDt4(data.sim[,'N0'][1:30], 0.15, 0, data.sim[,'f'][1:30], 10, 500/3)

nll4 = function(N, m, a, f, N0, V0, t, trace=TRUE){
	nll = -sum( dnorm(x = N,
		mean = E.Nt4(N0, m, a, f, t),
		sd = E.SDt4(N0, m, a, f, t, V0)$sd,
		log = TRUE
		))
	if(trace) {cat('m:', m, 'V0:', V0, 'a:', a, 'nll:', nll, '\n')}
	return(nll)
}
nll4(N=data.sim[,'N'], m=0.15, a=0, f=data.sim[,'f'], N0=data.sim[,'N0'], V0=500/3, t=10, trace=TRUE)

# choose a simulated dataset
data.sim=data.sim4

# run MLE
model1 = mle2(minuslogl = nll1, 
	start = list(m = 0.15, V0 = 113/3),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

model2 = mle2(minuslogl = nll2, 
	start = list(m = 0.15, b = 0.01, V0 = 113/3),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

model3 = mle2(minuslogl = nll3, 
	start = list(m = 0.15, c = 0.01, V0 = 113/3),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

model4 = mle2(minuslogl = nll4,
	start = list(m=0.15, V0 = 100, a = 0),
	data = list(N0 = data.sim[,'N0'], t = 10,
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='SANN')#, lower=c(m=0.0001, V0=213/3, a=0))

AICctab(model1, model2, model3, model4, nobs=nrow(data.sim))



# make model predictions
pfood = seq(1, 9.5, by=0.25)
p1 = E.Nt1(N0=113, t=10, m=coef(model1)['m'])
p2 = E.Nt2(N0=113, t=10, m=coef(model2)['m'], 
	b=coef(model2)['b'], f=pfood)
p3 = E.Nt3(N0=113, t=10, m=coef(model3)['m'], 
	c=coef(model3)['c'], f=pfood)


# plot data and predictions
plot(data.sim$N ~ data.sim$f, las=1, ylab='N', xlab='Food')
abline(h=p1, col='red')
points(p2 ~ pfood, col='blue', type='l')
points(p3 ~ pfood, col='orange', type='l')

