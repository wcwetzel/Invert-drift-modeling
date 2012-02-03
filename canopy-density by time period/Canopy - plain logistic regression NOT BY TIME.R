# plain logistic/binomial regression
# CANOPY AND DENSITY EXPERIMENTS
# CANOPY, density, predation, velocity
# 16 Jan 2012

# with pred_canopy data only, model with canopy and predation 
# is significantly better than model with predation only
# when we include density data, model with predation, density, canopy
# is only 0.1 AIC units better than predation and density (w/o canopy).

library(bbmle)
library(ggplot2)

data = read.csv(
"/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$rate = data$total / data$N0
datacd = data[-which(data$experiment=='excess_food'),]
datac = data[data$experiment=='pred_canopy',]
datad = data[data$experiment=='density',]
datad$stay = datad$N0 - datad$total


plot(rate ~ canopy, data=datac)
ggplot(data=data, aes(y=rate, x=canopy)) +
	geom_point() +
	stat_smooth() + stat_smooth(method='lm')

# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

# logit transform:
logit = function(x){
	log( x / (1 - x))
}

# the models
datacd = datacd
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

m7P.nopred = mle2(total ~ dbinom(size=N0,
	prob=logistic(m + b1 * N0 + b2 * canopy)), 
	start=list(m=0, b1=0, b2=0), data=datacd[datacd$preds==0,])

m7P.nopred.nocanopy = mle2(total ~ dbinom(size=N0,
	prob=logistic(m + b1 * N0)), 
	start=list(m=0, b1=0), data=datacd[datacd$preds==0,])
	
m7P.nopred.intercept = mle2(total ~ dbinom(size=N0,
	prob=logistic(m)), 
	start=list(m=0), data=datacd[datacd$preds==0,])

m7P.nopred.nodensity = mle2(total ~ dbinom(size=N0,
	prob=logistic(m + b2 * canopy)), 
	start=list(m=0, b2=0), data=datacd[datacd$preds==0,])


anova(m7P.nopred, m7P.nopred.nocanopy)
AICtab(m7P.nopred, m7P.nopred.nocanopy, m7P.nopred.intercept, m7P.nopred.nodensity)
anova(m7P.nopred.intercept, m7P.nopred.nodensity)



mP = mle2(total ~ dbinom(size=N0, 
	prob=logistic( m + b4 * preds)), 
	start=list(m=0, b4=0), data=datacd)

m2P  = mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b1 * N0 + b4 * preds)), 
	start=list(m=0, b1=0, b4=0), data=datacd)

m6P  = mle2(total ~ dbinom(size=N0, 
	prob=logistic(b1 * N0 + b2 * canopy + b4 * preds)), 
	start=list(b1=0, b2=0, b4=0), data=datacd)

mCanopyP =  mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b2 * canopy + b4 * preds)), 
	start=list(m=0, b2=0, b4=0), data=datacd)

mNOcanopyP =  mle2(total ~ dbinom(size=N0, 
	prob=logistic(m + b4 * preds)), 
	start=list(m=0, b4=0), data=datacd)




AICtab(m1, m2, m3, m4, m5, m6, m7, m7velocity, m7VP,
	m7P, mP, m2P, m6P, mCanopyP, mNOcanopyP, weights=TRUE)
anova(m7,m6)
anova(m7,m4)
anova(m7,m2)
anova(m7, m7velocity)
anova(m7VP, m7P)
anova(m7, m7P)
anova(m7P, m2P)
anova(mCanopyP, mNOcanopyP)




ggplot(data=datacd, aes(x=canopy, y=rate, colour=factor(preds))) + 
	geom_point() +
	stat_smooth(method='lm') 

quartz()
ggplot(data=datac, aes(x=canopy, y=rate, colour=factor(preds))) + 
	geom_point() +
	stat_smooth(method='lm') 
	
quartz()
ggplot(data=data, aes(x=canopy, y=rate, colour=factor(preds))) + 
	geom_point()  +
	stat_smooth(method='lm') 

quartz()
data.food1= data[data$food==1,]
ggplot(data=data.food1, aes(x=canopy, y=rate, colour=factor(preds))) + 
	geom_point()  +
	stat_smooth(method='lm') 

data.food1.nodensity= data.food1[!data$experiment=='density',]
ggplot(data=data.food1.nodensity, aes(x=canopy, y=rate, colour=factor(preds))) + 
	geom_point()  +
	stat_smooth(method='lm') 
	
ggplot(data=datacd[datacd$preds==0,], aes(x=canopy, y=rate, colour=factor(N0))) + 
	geom_point() +
	stat_smooth(method='lm') 




quarz()
plot(rate ~ canopy, pch=preds+10, data=datac)
newcan = 0:100
points()






# Ploting predictions from m7P
newN0 = rep(seq(0, 600, length = 7),6)
newcanopy = sort(rep(seq(0, 100, length = 6),7))

pppreds = logistic(coef(m7P)['m'] + coef(m7P)['b1'] * newN0 + 
	coef(m7P)['b2'] * newcanopy + coef(m7P)['b4'])
ppnopreds = logistic(coef(m7P)['m'] + coef(m7P)['b1'] * newN0 + 
	coef(m7P)['b2'] * newcanopy)
	
par(mfrow=c(1,1))
trellis.par.set("axis.line", list(col="transparent"))
print( wireframe(ppnopreds ~ newN0 + newcanopy,
				drape=TRUE,
				layout=c(0,1),
				xlab = list('Initial density', rot=-48),
				ylab='% Canopy',
				zlab=list('Drift rate', rot=90),
				col.regions = rev( gray(seq(0, 1, length=20))),
				at = seq(0, 0.253, length=20),
				aspect=c(1,1),
				screen = list(z = -245, x=-75),
				scales = list(arrows=FALSE,
				cex=0.5, col='black',
				font=3, tck=1),
				colorkey=FALSE
	))


ppnopreds.nopred = logistic(coef(m7P.nopred)['m'] + coef(m7P.nopred)['b1'] * newN0 + 
	coef(m7P.nopred)['b2'] * newcanopy)


trellis.par.set("axis.line", list(col="transparent"))
print( wireframe(ppnopreds.nopred ~ newN0 + newcanopy,
				drape=TRUE,
				layout=c(0,1),
				xlab = list('Initial density', rot=-48),
				ylab='% Canopy',
				zlab=list('Drift rate', rot=90),
				col.regions = 0,#rev( gray(seq(0, 1, length=20))),
				at = seq(0, 0.253, length=20),
				aspect=c(1,1),
				screen = list(z = -245, x=-75),
				scales = list(arrows=FALSE,
				cex=0.5, col='black',
				font=3, tck=1),
				colorkey=FALSE,
	))


###
ggplot(data=data[data$experiment=='density',], aes(x=N0, y=rate))+
	geom_point()+ stat_smooth(method='lm')
	
ggplot(data=data[data$experiment=='excess_food',], aes(x=food, y=rate))+
	geom_point() + stat_smooth()

