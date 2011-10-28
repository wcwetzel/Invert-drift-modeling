## likelihood for simple linear drift simulation
# assumes food density doesn't change
# Pr(drift) depends on initial food and current density
library(bbmle)
# 1st run 'simple drift stochastic simulation.R' to generate data

# Analytical models and NLLs for models


# Model1
E.Nt1 = function(N0, m, t){
	N0 * exp(-m * t)
}
E.Nt1(113, 0.1, 10)
E.SDt1 = function(N0, m, t, V0){
	sqrt( exp(2 * m * t) * (V0 + N0 * {exp(-3 * m * t) - 1} / 3) )
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






# choose a simulated dataset
data.sim=data.sim3

# run MLE
model1 = mle2(minuslogl = nll1, 
	start = list(m = 0.15, V0 = 56.9),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

model2 = mle2(minuslogl = nll2, 
	start = list(m = 0.15, b = 0.01, V0 = 56.9),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

model3 = mle2(minuslogl = nll3, 
	start = list(m = 0.15, c = 0.01, V0 = 56.9),
	data = list(N0 = 113, t = 10, 
	f = data.sim[,'f'], N = data.sim[,'N']),
	method='Nelder-Mead')

AICctab(model1, model2, model3)



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

