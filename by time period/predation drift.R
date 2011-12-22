### separate drift analysis for predation
# 19 Dec 2011

# first run 'prob model for drift.R'

datacp = read.csv("/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
datacp$rate = datacp$total / datacp$N0
datacp = datacp[which(datacp$experiment=='pred_canopy'),]
datacp$canopy = datacp$canopy/100
datacp$revcanopy = abs(1-datacp$canopy)


# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

modelp0 = mle2(total ~ dbinom(size = N0, prob = logistic(m)),
	start=list(m = -2), data=datacp)


modelp1 = mle2(total ~ dbinom(size = N0, prob = pmp1(p, m)),
	start=list(p=-10, m = -2), data=datacp)


modelp2 = mle2(total ~ dbinom(size = N0, prob = pmp2(c, m)),
	start=list(c=0, m = -2), data=datacp)


modelp3 = mle2(total ~ dbinom(size = N0, prob = pmp3(p, c, m)),
	start=list(p=-10, c=0, m = -2), data=datacp)

modelp4 = mle2(total ~ dbinom(size = N0, prob = pmp4(p, c, m, cp)),
	start= list(p=-10, c=0, m = -2, cp=0), data=datacp)

modelp5 = mle2(total ~ dbinom(size = N0, prob = pmp5(p, m, c)),
	start = list (p=0.1, m=0.3, c = 1), data=datacp)

AICtab(modelp0, modelp1, modelp2, modelp3, modelp4, modelp5)

anova(modelp4, modelp3)
anova(modelp3, modelp2)

#ploting

plot(rate ~ canopy, data=datacp, col=datacp$preds+1)
pcanopy = rep(seq(0,1, length=20), 2)
ppred = sort(rep(0:1, 20))

pcanopyrate1 = logistic(ppred * coef(modelp1)['p'] + coef(modelp1)['m'])
points(pcanopyrate1[ppred==0] ~ pcanopy[ppred==0], type='l', lwd=1.5, lty=3)
points(pcanopyrate1[ppred==1] ~ pcanopy[ppred==1], type='l', lwd=1.5, lty=3, col=2)


pcanopyrate2 = logistic(coef(modelp2)['c'] * pcanopy + coef(modelp2)['m'])
points(pcanopyrate2 ~ pcanopy, type='l', lwd=1.5, lty=1)

pcanopyrate3 = logistic(coef(modelp3)['p'] * ppred + coef(modelp3)['c'] * pcanopy + coef(modelp3)['m'])
points(pcanopyrate3[ppred==0] ~ pcanopy[ppred==0], type='l', col=1, lty=2)
points(pcanopyrate3[ppred==1] ~ pcanopy[ppred==1], type='l', col=2, lty=2)

pcanopyrate4 = logistic(coef(modelp4)['m'] + coef(modelp4)['p'] * 
	ppred + coef(modelp4)['c'] * pcanopy + coef(modelp4)['cp'] * ppred * pcanopy)
points(pcanopyrate4[ppred==0] ~ pcanopy[ppred==0], type='l', col=3, lty=4)
points(pcanopyrate4[ppred==1] ~ pcanopy[ppred==1], type='l', col=3, lty=4)



