#### simulate binomial drift by 4 time periods ####
# 12 Dec 2011

n = 30 # sample size
N0 = 113 # initial number at beginning of time period
m = 0.1
b = 0.01
c = 1

# models:
pm1 = function(m) {m} 
pm2 = function(m, b) {m - b * data$food}
pm3 = function(m, c) {m / (c + data$food)}
pm4 = function(m){m + a * data$}

# inital values and food:
data = data.frame(N0 = 113, food = rep(1:6,4))

# simulate model 1:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm1(m))
}
data.m1 = cbind(data, N, total=apply(N, 1, sum))
names(data.m1) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total')


# simulate model 2:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm2(m, b))
}
data.m2 = cbind(data, N, total=apply(N, 1, sum))
names(data.m2) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total')


# simulate model 3:
N = matrix(0, nrow=nrow(data), ncol=4) # vector of numbers drifting each t
for(i in 1:4){
	N[,i] = rbinom(nrow(data), N0 - apply(N, 1, sum), pm3(m, c))
}
data.m3 = cbind(data, N, total=apply(N, 1, sum))
names(data.m3) = c('N0', 'food', 't1', 't2', 't3', 't4', 'total')




