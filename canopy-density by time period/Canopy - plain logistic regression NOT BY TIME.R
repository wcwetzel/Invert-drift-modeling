# plain logistic/binomial regression
# CANOPY AND DENSITY EXPERIMENTS
# CANOPY, density, predation, velocity
# 16 Jan 2012

library(bbmle)

data = read.csv("/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$rate = data$total / data$N0
datacd = data[-which(data$experiment=='excess_food'),]

# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

# logit transform:
logit = function(x){
	log( x / (1 - x))
}

# the models

m1 = mle2(total ~ dbinom(size=N0, prob=logistic(m)), start=list(m=0), data=datacd)

m2 = mle2(total ~ dbinom(size=N0, prob=logistic(m + b1 * N0)), 
	start=list(m=0, b1=0), data=datacd)

m3 = mle2(total ~ dbinom(size=N0, prob=logistic(b1 * N0)), 
	start=list(b1=0), data=datacd)

m4 = mle2(total ~ dbinom(size=N0, prob=logistic(m + b2 * canopy)), 
	start=list(m=0, b2=0), data=datacd)
	
m5 = mle2(total ~ dbinom(size=N0, prob=logistic(b2 * canopy)), 
	start=list(b2=0), data=datacd)

m6 = mle2(total ~ dbinom(size=N0, prob=logistic(b1 * N0 + b2 * canopy)), 
	start=list(b1=0, b2=0), data=datacd)
	
m7 = mle2(total ~ dbinom(size=N0, prob=logistic(m + b1 * N0 + b2 * canopy)), 
	start=list(m=0, b1=0, b2=0), data=datacd)

m7velocity = mle2(total ~ dbinom(size=N0, 
prob=logistic(m + b1 * N0 + b2 * canopy + b3 * velocity)), 
	start=list(m=0, b1=0, b2=0, b3=0), data=datacd)

m7VP = mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b1 * N0 + b2 * canopy + b3 * velocity + b4 * preds)), 
	start=list(m=0, b1=0, b2=0, b3=0, b4=0), data=datacd)

m7P  = mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b1 * N0 + b2 * canopy + b4 * preds)), 
	start=list(m=0, b1=0, b2=0, b4=0), data=datacd)


mP = mle2(total ~ dbinom(size=N0, 
	prob=logistic( m + b4 * preds)), 
	start=list(m=0, b4=0), data=datacd)

m2P  = mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b1 * N0 + b4 * preds)), 
	start=list(m=0, b1=0, b4=0), data=datacd)

m6P  = mle2(total ~ dbinom(size=N0, 
	prob=logistic(b1 * N0 + b2 * canopy + b4 * preds)), 
	start=list(b1=0, b2=0, b4=0), data=datacd)





AICtab(m1, m2, m3, m4, m5, m6, m7, m7velocity, m7VP, m7P, mP, m2P, m6P, weights=TRUE)
anova(m7,m6)
anova(m7,m4)
anova(m7,m2)
anova(m7, m7velocity)
anova(m7VP, m7P)
anova(m7, m7P)
anova(m7P, m2P)

# Ploting predictions from m7P
newN0 = rep(seq(0, 600, length = 10),10)
newcanopy = sort(rep(seq(0, 100, length = 10),10))

pppreds = logistic(coef(m7P)['m'] + coef(m7P)['b1'] * new )
ppnopreds
