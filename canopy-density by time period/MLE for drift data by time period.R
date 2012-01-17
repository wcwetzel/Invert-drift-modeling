#### MLE for excess food drift with time.csv ####
# 12 Dec 2011
library(bbmle)

# FIRST run 'prob models for drift.R'

data = read.csv("/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$rate = data$total / data$N0

datacp = data
dataexcessfood = data[which(data$experiment=='excess_food'),]
data = data[-which(data$experiment=='pred_canopy'),]


# neg log L for model 1
nll1 = function(m){
	with(data,
	-sum( c(dbinom(t1, N0, pm1(m), log=TRUE),
		dbinom(t2, N0 - t1, pm1(m), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm1(m), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm1(m), log=TRUE)))
	)
}

nll1(.1)

# neg log L for model 2
nll2 = function(m, b){
	with(data,
	-sum( c(dbinom(t1, N0, pm2(m, b), log=TRUE),
		dbinom(t2, N0 - t1, pm2(m, b), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm2(m, b), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm2(m, b), log=TRUE)))
	)
}

nll2(.1, 0.01)

# neg log L for model 3
nll3 = function(m, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm3(m, c), log=TRUE),
		dbinom(t2, N0 - t1, pm3(m, c), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm3(m, c), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm3(m, c), log=TRUE)))
	)
}

nll3(.1, 1)

# neg log L for model 3 WITH k!!!
nll3k = function(m, k, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm3k(m, k, c), log=TRUE),
		dbinom(t2, N0 - t1, pm3k(m, k, c), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm3k(m, k, c), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm3k(m, k, c), log=TRUE)))
	)
}

nll3k(.1, 0.01, 1)


# neg log L for model 4
nll4 = function(m, a){
	with(data,
	-sum( c(dbinom(t1, N0, pm4(m, a, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm4(m, a, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm4(m, a, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm4(m, a, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll4(.1, 0.001)

# neg log L for model 5
nll5 = function(m, a, b){
	with(data,
	-sum( c(dbinom(t1, N0, pm5(m, a, b, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm5(m, a, b, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm5(m, a, b, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm5(m, a, b, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll5(.1, 0.001, 0)

# neg log L for model 6
nll6 = function(m, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6(m, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6(m, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6(m, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6(m, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll6(.1, 0.01, 1)

# neg log L for model 6 WITH k!!!
nll6k = function(m, k, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6k(m, k, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6k(m, k, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6k(m, k, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6k(m, k, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll6k(.1, 0.01, 0.0001, 1)

# neg log L for model 6 WITH k!!! and no intercept
nll6kn = function(m, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6kn(m, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6kn(m, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6kn(m, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6kn(m, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll6kn(-2, 0.0001, 1)

# neg log L for model 6 without parm m, ie, no inadvert drift
nll6n = function(a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6n(a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6n(a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6n(a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6n(a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

# neg log L for model 7
nll7 = function(m, a, c, alpha){
	with(data,
	-sum( c(dbinom(t1, N0, pm7(m, a, c, alpha, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm7(m, a, c, alpha, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm7(m, a, c, alpha, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm7(m, a, c, alpha, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

# neg log L for model 7
nll7k = function(m, k, a, c, alpha){
	with(data,
	-sum( c(dbinom(t1, N0, pm7k(m, k, a, c, alpha, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm7k(m, k, a, c, alpha, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm7k(m, k, a, c, alpha, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm7k(m, k, a, c, alpha, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}


# neg log L for model 8
nll8 = function(m, a, b, alpha){
	with(data,
	-sum( c(dbinom(t1, N0, pm8(m, a, b, alpha, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm8(m, a, b, alpha, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm8(m, a, b, alpha, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm8(m, a, b, alpha, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

# neg log L for model 9
nll9 = function(m, a, alpha){
	with(data,
	-sum( c(dbinom(t1, N0, pm9(m, a, alpha, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm9(m, a, alpha, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm9(m, a, alpha, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm9(m, a, alpha, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

# neg log L for model 1V
nll1v = function(m, m1){
	with(data,
	-sum( c(dbinom(t1, N0, pm1v(m, m1), log=TRUE),
		dbinom(t2, N0 - t1, pm1v(m, m1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm1v(m, m1), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm1v(m, m1), log=TRUE)))
	)
}

nll1(.1)

# neg log L for model 2V
nll2v = function(m, m1, b){
	with(data,
	-sum( c(dbinom(t1, N0, pm2v(m, m1, b), log=TRUE),
		dbinom(t2, N0 - t1, pm2v(m, m1, b), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm2v(m, m1, b), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm2v(m, m1, b), log=TRUE)))
	)
}

nll2(.1, 0.01)

# neg log L for model 3V
nll3v = function(m, m1, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm3v(m, m1, c), log=TRUE),
		dbinom(t2, N0 - t1, pm3v(m, m1, c), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm3v(m, m1, c), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm3v(m, m1, c), log=TRUE)))
	)
}

nll3(.1, 0.01)

# neg log L for model 4V
nll4v = function(m, m1, a){
	with(data,
	-sum( c(dbinom(t1, N0, pm4v(m, m1, a, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm4v(m, m1, a, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm4v(m, m1, a, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm4v(m, m1, a, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll4(.1, 0.001)

# neg log L for model 5V
nll5v = function(m, m1, a, b){
	with(data,
	-sum( c(dbinom(t1, N0, pm5v(m, m1, a, b, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm5v(m, m1, a, b, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm5v(m, m1, a, b, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm5v(m, m1, a, b, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

nll5(.1, 0.001, 0)

# neg log L for model 6V
nll6v = function(m, m1, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6v(m, m1, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6v(m, m1, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6v(m, m1, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6v(m, m1, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}

# neg log L for model 6kV
nll6kv = function(m, m1, k, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6kv(m, m1, k, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6kv(m, m1, k, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6kv(m, m1, k, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6kv(m, m1, k, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}


# neg log L for model 3 negexp
nll3negexp = function(m, a, b){
	with(data,
	-sum( c(dbinom(t1, N0, pm3negexp(m, a, b), log=TRUE),
		dbinom(t2, N0 - t1, pm3negexp(m, a, b), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm3negexp(m, a, b), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm3negexp(m, a, b), log=TRUE)))
	)
}

# neg log L for model 6p
nll6p = function(p, m, a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6p(p, m, a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6p(p, m, a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6p(p, m, a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6p(p, m, a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}
nll6p(0, 0.1, 0.0001, 1)


# set simulated data set, comment out for real data loaded above
#data=data.m6

# initial guesses for mle
mguess = logit(0.07)
bguess = 0.0002
cguess = 1.2
aguess = 0.00001


model1 = mle2(nll1, start = list(m = -2.966), data=data)
model1s = mle2(total ~ dbinom(size=N0, prob = pm1(m)), 
	start = list(m = -1.50397), data=data)


model2 = mle2(nll2, start = list(m = -2.6083, b=0.009810448), 
	data=data, method='Nelder-Mead')
model2s = mle2(total ~ dbinom(size=N0, prob = pm2(m, b)), 
	start = list(m = -1.09635883, b=0.03515887), data=data)

model3 = mle2(nll3, start = list(m = -2.3281463, c=0.4366301), 
	data=data, method='Nelder-Mead')
model3k = mle2(nll3k, start = list(m = -2.3281463, k = 0.01, c=0.4366301), 
	data=data, method='Nelder-Mead')

model3s = mle2(total ~ dbinom(size=N0, prob = pm3(m, c)), 
	start = list(m = -0.5995892, c= 0.5724655), 
	method = 'Nelder-Mead', data=data)

model4 = mle2(nll4, start = list(m = -2.672963, a=-5.55209e-06), 
	data=data, method='Nelder-Mead')
model4s = mle2(total ~ dbinom(size=N0, prob = pm4(m, a, N=N0)), 
	start = list(m = -2.7181420912, a = 0.0003156788), 
	data=data, method='Nelder-Mead')


model5 = mle2(nll5, start = list(m = -2.586689, a=0.000100, b=0.0002), 
	data=data, method='Nelder-Mead')
model5s = mle2(total ~ dbinom(size=N0, prob = pm5(m, a, b, N=N0)), 
	start = list(m = -2.6148440037, a = 0.0005968744, b = 0.0254360974), 
	data=data, method='Nelder-Mead')


model6 = mle2(nll6, start = list(m = -2.6535, a=0.0001886, c=1.35064), 
	data=data, method='Nelder-Mead')
model6k = mle2(nll6k, start = list(m = -6.0876078754, k=-2.6427495764, a=0.0001177, c=0.9006994661), 
	data=data, method='Nelder-Mead')
model6kn = mle2(nll6kn, start = list(m = -3.8214296773, a=0.0000358755, c=-0.639032), 
	data=data, method='SANN')
model6s = mle2(total ~ dbinom(size=N0, prob = pm6(m, a, c, N=N0)), 
	start = list(m = -2.3497323009, a = 0.0009231747, c = 1.1825592849), 
	data=data, method='Nelder-Mead')
model6n = mle2(nll6n, start = list(a= 0.0003762844, c= 1.2377230986),
	data=data, method='Nelder-Mead') # no inadvertant drift


model6nnc = with(data, 
	 optimize(nll6n, lower = 0.00000000001, upper = 0.001, c=0)
)

model7 = mle2(nll7, start = list(m = -2.4898977683, 
	a= 0.0001072355, c= 1.1130006792, alpha= 0.919995347), 
	data=data, method='SANN')

model7k = mle2(nll7k, start = list(m = -7.0522809542, k=-2.7231198215, a= 0.0000683319, c= 1.0152550536, alpha= 0.7726926992), 
	data=data, method='SANN')

model8 = mle2(nll8, start = list(m = mguess, a=aguess, b= bguess, alpha=1), 
	data=data, method='Nelder-Mead')

model9 = mle2(nll9, start = list(m = mguess, a=aguess, alpha=1), 
	data=data, method='Nelder-Mead')

model1v = mle2(nll1v, start = list(m = -2.966, m1=0), data=data)
model2v = mle2(nll2v, start = list(m = -2.6083, m1=0, b=0.009810448), 
	data=data, method='Nelder-Mead')
model3v = mle2(nll3v, start = list(m = -2.3281463, m1=0, c=0.4366301), 
	data=data, method='Nelder-Mead')
model4v = mle2(nll4v, start = list(m = -2.672963, m1=0, a=-5.55209e-06), 
	data=data, method='Nelder-Mead')
model5v = mle2(nll5v, start = list(m = -2.586689, m1=0, a=0.000100, b=0.0002), 
	data=data, method='Nelder-Mead')
model6v = mle2(nll6v, start = list(m = -2.6535, m1=0, a=0.0001886, c=1.35064), 
	data=data, method='Nelder-Mead')

model6kv = mle2(nll6kv, start = list(m = -8.2618617852, m1= 2.7754386851, k=-2.5009778563, a= 0.0001355539, c= 1.1824372505), 
	data=data, method='Nelder-Mead')

model3ne = mle2(nll3negexp, start = list(m = -2, a = 0.5, b = 0.25), 
	data=data, method='Nelder-Mead')

model3ne.EF = mle2(nll3negexp, start = list(m = -2, a = 0.5, b = 0.25), 
	data=dataexcessfood, method='Nelder-Mead')

model3k.EF = mle2(nll3k, start = list(m = -2.3281463, k = 0.01, c=0.4366301), 
	data=dataexcessfood, method='Nelder-Mead')
	
model3.EF = mle2(nll3, start = list(m = -2.3281463, c=0.4366301), 
	data=dataexcessfood, method='Nelder-Mead')

model6p = mle2(nll6p, start = list(p = 0, m = -2.6535, a=0.0001886, c=1.35064), 
	data=data, method='Nelder-Mead')

# model3 with just excess food


AICtab(model1, model2, model3, model3k, model4, model5, model6, model6k, model6kn, model6n, model7, model7k, model8, model9, model3ne, weights=TRUE)
AICtab(model1, model2, model3, model3k, model4, model5, model6, model6k, model6kn, model6n, model7, model7k, model8, model9, model1v, model2v, model3v, model4v, model5v, model6v, model6kv, weights=TRUE)
AICtab(model3ne.EF, model3k.EF)
BICtab(model1, model2, model3, model3k, model4, model5, model6, model6k, model6kn, model6n, model7, model7k, model8, model9, weights=TRUE)



AICctab(model1s, model2s, model3s, model4s, model5s, model6s, nobs=nrow(data), weights=TRUE)
AICtab(model1s, model2s, model3s, model4s, model5s, model6s)
BICtab(model1s, model2s, model3s, model4s, model5s, model6s)

#prom6 = profile(model6)
prom7k = profile(model7k)

