### MLE for drift by time period ###
# 12 Dec 2011

# 1ST RUN 'drift simulaiton by time period.R'

library(bbmle)

N0 = 113
nll1 = function(m){
	-sum( c(dbinom(N[,1], N0, pm1(m), log=TRUE),
		dbinom(N[,2], N0 - N[,1], pm1(m), log=TRUE),
		dbinom(N[,3], N0 - apply(N[,1:2], 1, sum), pm1(m), log=TRUE),
		dbinom(N[,4], N0 - apply(N[,1:3], 1, sum), pm1(m), log=TRUE)))
}

nll1(.1)

model1 = mle2(nll1, start = list(m = 0.2))
summary(model1)
confint(model1)



