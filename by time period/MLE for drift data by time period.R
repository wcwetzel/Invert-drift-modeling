#### MLE for excess food drift with time.csv ####
# 12 Dec 2011
library(bbmle)

# FIRST run 'prob models for drift.R'

data = read.csv("/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
#data$N0 = 213
data$rate = data$total / data$N0

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

nll3(.1, 0.01)

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

# neg log L for model 6 without parm m
nll6n = function(a, c){
	with(data,
	-sum( c(dbinom(t1, N0, pm6n(a, c, N0), log=TRUE),
		dbinom(t2, N0 - t1, pm6n(a, c, N0 - t1), log=TRUE),
		dbinom(t3, N0 - apply(cbind(t1, t2), 1, sum), pm6n(a, c, N0 - apply(cbind(t1, t2), 1, sum)), log=TRUE),
		dbinom(t4, N0 - apply(cbind(t1, t2, t3), 1, sum), pm6n(a, c, N0 - apply(cbind(t1, t2, t3), 1, sum)), log=TRUE)))
	)
}


# set simulate data set, comment out for real data loaded above
#data=data.m6

# initial guesses for mle
mguess = logit(0.07)
bguess = 0.002
cguess = 1.2
aguess = 0.0001


model1 = mle2(nll1, start = list(m = mguess), data=data)
model1s = mle2(total ~ dbinom(size=N0, prob = pm1(m)), 
	start = list(m = mguess), data=data)


model2 = mle2(nll2, start = list(m = mguess, b=bguess), data=data)
model2s = mle2(total ~ dbinom(size=N0, prob = pm2(m, b)), 
	start = list(m = mguess, b=bguess), data=data)

model3 = mle2(nll3, start = list(m = mguess, c=cguess), data=data, method='SANN')
model3s = mle2(total ~ dbinom(size=N0, prob = pm3(m, c)), 
	start = list(m = mguess, c= cguess), data=data)

model4 = mle2(nll4, start = list(m = mguess, a=aguess), data=data, method='SANN')
model4s = mle2(total ~ dbinom(size=N0, prob = pm4(m, a, N=N0)), 
	start = list(m = mguess, a = aguess), data=data, method='SANN')


model5 = mle2(nll5, start = list(m = mguess, a=aguess, b=bguess), 
	data=data, method='SANN')
model5s = mle2(total ~ dbinom(size=N0, prob = pm5(m, a, b, N=N0)), 
	start = list(m = mguess, a = aguess, b = bguess), data=data, method='SANN')


model6NM = mle2(nll6, start = list(m = -2.6535, a=0.0001886, c=1.35064), 
	data=data, method='Nelder-Mead')
model6s = mle2(total ~ dbinom(size=N0, prob = pm6(m, a, c, N=N0)), 
	start = list(m = mguess, a = aguess, c = cguess), data=data, method='SANN')
model6n = mle2(nll6n, start = list(a=aguess, c=cguess),
	data=data, method='SANN')


AICctab(model1, model2, model3, model4, model5, model6, model6n,
	nobs=nrow(data), weights=TRUE)
AICtab(model1, model2, model3, model4, model5, model6, model6n)
BICtab(model1, model2, model3, model4, model5, model6, model6n)

AICctab(model1s, model2s, model3s, model4s, model5s, model6s, nobs=30, weights=TRUE)
AICtab(model1s, model2s, model3s, model4s, model5s, model6s)
BICtab(model1s, model2s, model3s, model4s, model5s, model6s)