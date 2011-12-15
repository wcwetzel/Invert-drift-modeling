#### simulate binomial drift by 4 time periods ####
# 12 Dec 2011


# FIRST run 'prob models for drift.R'

n = 6 # sample size
foodlevels = 6 # number of food levels
N0 = 113 # initial number at beginning of time period
m = logit(0.1)
b = 0.2
c = 1
a = 0.01





# inital values and food:
data = data.frame(N0 = sort(rep(seq(113, 213, length=n), foodlevels)), food = rep(1:foodlevels,n))

# simulate model 1:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm1(m))
}
data.m1 = cbind(data, N, total=apply(N, 1, sum))
data.m1 = cbind(data.m1, rate = data.m1$total/data.m1$N0)
names(data.m1) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')


# simulate model 2:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm2(m, b))
}
data.m2 = cbind(data, N, total=apply(N, 1, sum))
data.m2 = cbind(data.m2, rate = data.m2$total/data.m2$N0)
names(data.m2) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')


# simulate model 3:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm3(m, c))
}
data.m3 = cbind(data, N, total=apply(N, 1, sum))
data.m3 = cbind(data.m3, rate = data.m3$total/data.m3$N0)
names(data.m3) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')

# simulate model 4:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm4(m, a, N0 - apply(N, 1, sum)))
}
data.m4 = cbind(data, N, total=apply(N, 1, sum))
data.m4 = cbind(data.m4, rate = data.m4$total/data.m4$N0)
names(data.m4) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')

# simulate model 5:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm5(m, a, b, N0 - apply(N, 1, sum)))
}
data.m5 = cbind(data, N, total=apply(N, 1, sum))
data.m5 = cbind(data.m5, rate = data.m5$total/data.m5$N0)
names(data.m5) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')

# simulate model 6:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm6(m, a, c, N0 - apply(N, 1, sum)))
}
data.m6 = cbind(data, N, total=apply(N, 1, sum))
data.m6 = cbind(data.m6, rate = data.m6$total/data.m6$N0)
names(data.m6) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total', 'rate')
