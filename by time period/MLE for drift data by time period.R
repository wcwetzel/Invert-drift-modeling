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










model1 = mle2(nll1, start = list(m = 0.2), data=data)

model2 = mle2(nll2, start = list(m = 0.2, b=0), data=data)





