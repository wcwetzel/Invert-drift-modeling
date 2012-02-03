# Beta regression
# CANOPY AND DENSITY EXPERIMENTS
# 27 Jan 2012

library(bbmle)
library(ggplot2)

data = read.csv(
"/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$rate = data$total / data$N0
datacd = data[-which(data$experiment=='excess_food'),]
datac = data[data$experiment=='pred_canopy',]
datad = data[data$experiment=='density',]

hist(data$rate, xlim=c(0,1), prob=TRUE)
lines(density(data$rate))

logistic = function(x) {1 / (1 + exp(-x))}

# mu functions
mu0 = function(b0){logistic(b0)}

mu1 = function(b0, b1, N0){logistic(b0 + b1 * N0)}

mu2 = function(b0, b2, canopy){logistic(b0 + b2 * canopy)}

mu3 = function(b0, b1, N0, b2, canopy){logistic(b0 + b1 * N0 + b2 * canopy)}

mu0p = function(b0, p, preds){logistic(b0 + p * preds)}

mu1p = function(b0, b1, N0, p, preds){logistic(b0 + b1 * N0 + p * preds)}

mu2p = function(b0, b2, canopy, p, preds){logistic(b0 + b2 * canopy + p * preds)}

mu3p = function(b0, b1, N0, b2, canopy, p, preds){logistic(b0 + b1 * N0 + b2 * canopy + p * preds)}

mu2pint = function(b0, b2, canopy, p, preds, b3){logistic(b0 + b2 * canopy + p * preds + preds * b3 * canopy)}

mu4 = function(b0, b4, food){logistic(b0 + b4 * food)}

mu5 = function(b0, a, c, food){logistic(b0 + a/(c + food) )}

# also add models without intercept

plot(rate ~ N0, data=datacd)

# MODELS WITH COMBINED DENSITY/CANOPY DATASET - bad idea
m0 = mle2(rate ~ dbeta(shape1=mu0(b0) * theta, 
	shape2=(1-mu0(b0)) * theta), 
	start=list(b0=0.5, theta=10), 
	data=datacd, method='Nelder-Mead')
abline(h = logistic(coef(m0)))

m1 = mle2(rate ~ dbeta(shape1=mu1(b0, b1, N0) * theta, 
	shape2=(1-mu1(b0, b1, N0)) * theta), 
	start=list(b0=0.5, theta=10, b1=0), 
	data=datacd, method='Nelder-Mead')
m1.predicted = logistic(coef(m1)['b0'] + coef(m1)['b1'] * 100:600)
points(100:600, m1.predicted, type='l', col='red')

m2 = mle2(rate ~ dbeta(shape1=mu2(b0, b2, canopy) * theta, 
	shape2=(1-mu2(b0, b2, canopy)) * theta), 
	start=list(b0=0.5, theta=10, b2=0), 
	data=datacd, method='Nelder-Mead')
m2.predicted = logistic(coef(m2)['b0'] + coef(m2)['b2'] * 0:100)
plot(rate ~ canopy, data=datacd)
points(0:100, m2.predicted, type='l', col='red')

m3 = mle2(rate ~ dbeta(shape1=mu3(b0, b1, N0, b2, canopy) * theta, 
	shape2=(1-mu3(b0, b1, N0, b2, canopy)) * theta), 
	start=list(b0=0.5, theta=10, b2=0, b1=0), 
	data=datacd, method='Nelder-Mead')
#m3.predicted = logistic(coef(m3)['b0'] + coef(m3)['b1'] * 
	#rep(seq(100,600,length=100),100) + 
	#coef(m3)['b2'] * sort(rep(1:100, 100)))
#plot(rate ~ canopy, data=datacd)
#points(0:100, m3.predicted, type='l', col='red')

AICtab(m0, m1, m2, m3)


# models for CANOPY dataset
m0 = mle2(rate ~ dbeta(shape1=mu0(b0) * theta, 
	shape2=(1-mu0(b0)) * theta), 
	start=list(b0=0.5, theta=10), 
	data=datac, method='Nelder-Mead')

m0p = mle2(rate ~ dbeta(shape1=mu0p(b0, p, preds) * theta, 
	shape2=(1-mu0p(b0, p, preds)) * theta), 
	start=list(b0=0.5, theta=10, p=0), 
	data=datac, method='Nelder-Mead')

m1 = mle2(rate ~ dbeta(shape1=mu2(b0, b2, canopy) * theta, 
	shape2=(1-mu2(b0, b2, canopy)) * theta), 
	start=list(b0=0.5, theta=10, b2=0), 
	data=datac, method='Nelder-Mead')
m1.predicted = logistic(coef(m1)['b0'] + coef(m1)['b2'] * 0:100)

m1p = mle2(rate ~ dbeta(shape1=mu2p(b0, b2, canopy, p, preds) * theta, 
	shape2=(1-mu2p(b0, b2, canopy, p, preds)) * theta), 
	start=list(b0=0.5, theta=10, b2=0, p=0), 
	data=datac, method='Nelder-Mead')
m1p.nopred.predicted = logistic(coef(m1p)['b0'] + coef(m1p)['b2'] * 0:100)
m1p.pred.predicted = logistic(coef(m1p)['b0'] + coef(m1p)['b2'] * 0:100 + coef(m1p)['p'])

m2p = mle2(rate ~ dbeta(shape1=mu2pint(b0, b2, canopy, p, preds, b3) * theta, 
	shape2=(1-mu2pint(b0, b2, canopy, p, preds, b3)) * theta), 
	start=list(b0=-0.77, theta=10, b2=-0.009, p=0.65, b3=0.013), 
	data=datac, method='Nelder-Mead')


AICtab(m0, m0p, m1, m1p, m2p)
anova(m0, m0p)
anova(m0p, m1p)
anova(m1, m0)


#m1p.predicted = logistic(coef(m2)['b0'] + coef(m2)['b2'] * 0:100)
#points(0:100, m2.predicted, type='l', col='red')


plot(rate ~ canopy, data=datac, col=datac$preds+1, pch=20)
#abline(h = logistic(coef(m0)))
abline(h = logistic(coef(m0p)['b0']), lty=2)
abline(h = logistic(coef(m0p)['b0'] + coef(m0p)['p']), col='red', lty=2)
points(0:100, m1.predicted, type='l', col='grey')
points(0:100, m1p.nopred.predicted, type='l')
points(0:100, m1p.pred.predicted, type='l', col='red')






# models for DENSITY dataset
m0 = mle2(rate ~ dbeta(shape1=mu0(b0) * theta, 
	shape2=(1-mu0(b0)) * theta), 
	start=list(b0=0.5, theta=10), 
	data=datad, method='Nelder-Mead')

m1 = mle2(rate ~ dbeta(shape1=mu1(b0, b1, N0) * theta, 
	shape2=(1-mu1(b0, b1, N0)) * theta), 
	start=list(b0=0.5, theta=10, b1=0), 
	data=datad, method='Nelder-Mead')
m1.predicted = logistic(coef(m1)['b0'] + coef(m1)['b1'] * 100:600)


plot(rate ~ N0, data=datad, pch=20)
points(100:600, m1.predicted, type='l', col='red')

anova(m0,m1)
AICtab(m0, m1)


# models for FOOD dataset