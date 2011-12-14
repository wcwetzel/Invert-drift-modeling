#### MLE for excess food drift with time.csv ####
# 12 Dec 2011
library(bbmle)

#data = read.csv("/Users/will/Documents/Analysis for colleagues/bruce/Drift data/excess food drift with time.csv")
#data$N0 = 113

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


data=data.m6

model1 = mle2(nll1, start = list(m = 0.2), data=data)

model2 = mle2(nll2, start = list(m = 0.2, b=0), data=data)

model3 = mle2(nll3, start = list(m = 0.2, c=2), data=data)

model4 = mle2(nll4, start = list(m = 0.15, a=0.001), data=data)

model5 = mle2(nll5, start = list(m = 0.15, a=0.001, b=0), data=data)

model6 = mle2(nll6, start = list(m = 0.15, a=0.001, c=1), data=data)

AICctab(model1, model2, model3, model4, model5, model6)
