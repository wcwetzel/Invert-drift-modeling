### simulation of model predictions ###
# 17 Dec 2011


# FIRST run 'prob models for drift.R'
# then run 'MLE for drift data by time period.R'

###############################
## models with logisitic(m): ##
## predicted with predicted food:

# constant rate
ppm1 = function(m) {logistic( m )} 

# linear decrease with food
ppm2 = function(m, b) {logistic( m ) - b * pdata$food } 

# rational decrease with food
ppm3 = function(m, c) {logistic( m ) / (c + pdata$food) } 

# linear increase with density
ppm4 = function(m, a, N) {logistic( m ) + a * N } 

# linear decrease with food, increase with density
ppm5 = function(m, a, b, N) {logistic( m ) + a * N - b * pdata$food }
 
# rational decrease with food, linear increase with density
ppm6 = function(m, a, c, N) { (logistic(m) + a * N) / (c + pdata$food) }

# rational decrease with food, linear increase with density
ppm6n = function(a, c, N) { (a * N) / (c + pdata$food) }

# rational decrease with food, linear increase with density WITH k!!!
ppm6k = function(m, k, a, c, N) { logistic(m) + {(logistic(k) + a * N) / (c + pdata$food)} }



# inital values and food:
pdata = data.frame(N0 = sort(rep(c(113, 163, 213, 263, 313, 363, 413, 463, 513, 563, 613), 21)), food = rep(seq(1,6,by=0.25),11))

# simulate model 1:
N = matrix(0, nrow=nrow(pdata), ncol=4) # vector of numbers drifting each t

N[,1] = pdata$N0 * ppm1(coef(model1)['m'])
N[,2] = (pdata$N0 - N[,1]) * ppm1(coef(model1)['m'])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm1(coef(model1)['m'])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm1(coef(model1)['m'])
pdata$totalm1 = apply(N, 1, sum)
pdata$ratem1 = pdata$totalm1/pdata$N0


# simulate model 2:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm2(coef(model2)['m'], coef(model2)['b'])
N[,2] = (pdata$N0 - N[,1]) * ppm2(coef(model2)['m'], coef(model2)['b'])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm2(coef(model2)['m'], coef(model2)['b'])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm2(coef(model2)['m'], coef(model2)['b'])
pdata$totalm2 = apply(N, 1, sum)
pdata$ratem2 = pdata$totalm2/pdata$N0

# simulate model 3:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm3(coef(model3)['m'], coef(model3)['c'])
N[,2] = (pdata$N0 - N[,1]) * ppm3(coef(model3)['m'], coef(model3)['c'])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm3(coef(model3)['m'], coef(model3)['c'])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm3(coef(model3)['m'], coef(model3)['c'])
pdata$totalm3 = apply(N, 1, sum)
pdata$ratem3 = pdata$totalm3/pdata$N0


# simulate model 4:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm4(coef(model4)['m'], coef(model4)['a'], pdata$N0)
N[,2] = (pdata$N0 - N[,1]) * ppm4(coef(model4)['m'], coef(model4)['a'], pdata$N0 - N[,1])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm4(coef(model4)['m'], coef(model4)['a'], pdata$N0 - N[,1] - N[,2])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm4(coef(model4)['m'], coef(model4)['a'], pdata$N0 - N[,1] - N[,2] - N[,3])
pdata$totalm4 = apply(N, 1, sum)
pdata$ratem4 = pdata$totalm4/pdata$N0

# simulate model 5:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm5(coef(model5)['m'], coef(model5)['a'], coef(model5)['b'], pdata$N0)
N[,2] = (pdata$N0 - N[,1]) * ppm5(coef(model5)['m'], coef(model5)['a'], coef(model5)['b'], pdata$N0 - N[,1])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm5(coef(model5)['m'], coef(model5)['a'], coef(model5)['b'], pdata$N0 - N[,1] - N[,2])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm5(coef(model5)['m'], coef(model5)['a'], coef(model5)['b'], pdata$N0 - N[,1] - N[,2] - N[,3])
pdata$totalm5 = apply(N, 1, sum)
pdata$ratem5 = pdata$totalm5/pdata$N0


# simulate model 6:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm6(coef(model6)['m'], coef(model6)['a'], coef(model6)['c'], pdata$N0)
N[,2] = (pdata$N0 - N[,1]) * ppm6(coef(model6)['m'], coef(model6)['a'], coef(model6)['c'], pdata$N0 - N[,1])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm6(coef(model6)['m'], coef(model6)['a'], coef(model6)['c'], pdata$N0 - N[,1] - N[,2])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm6(coef(model6)['m'], coef(model6)['a'], coef(model6)['c'], pdata$N0 - N[,1] - N[,2] - N[,3])
pdata$totalm6 = apply(N, 1, sum)
pdata$ratem6 = pdata$totalm6/pdata$N0


# simulate model 6k:
N = matrix(0, nrow=nrow(pdata), ncol=4) 
N[,1] = pdata$N0 * ppm6k(coef(model6k)['m'], coef(model6k)['k'], coef(model6k)['a'], coef(model6k)['c'], pdata$N0)
N[,2] = (pdata$N0 - N[,1]) * ppm6k(coef(model6k)['m'], coef(model6k)['k'], coef(model6k)['a'], coef(model6k)['c'], pdata$N0 - N[,1])
N[,3] = (pdata$N0 - N[,1] - N[,2]) * ppm6k(coef(model6k)['m'], coef(model6k)['k'], coef(model6k)['a'], coef(model6k)['c'], pdata$N0 - N[,1] - N[,2])
N[,4] = (pdata$N0 - N[,1] - N[,2] - N[,3]) * ppm6k(coef(model6k)['m'], coef(model6k)['k'], coef(model6k)['a'], coef(model6k)['c'], pdata$N0 - N[,1] - N[,2] - N[,3])
pdata$totalm6k = apply(N, 1, sum)
pdata$ratem6k = pdata$totalm6k/pdata$N0
