## long format data binomial GLMM
########## separate analyses by species ###########
# 13 April 2012
# 1st run 'drift data reshape.R' to put the data
# into the long format, in which each row is a spot check
# 2nd run 'drift data reshape by species.R'
# to do the same by species

library(lme4)
library(ggplot2)
library(rethinking)
library(bbmle)

# logistic transform:
logistic = function(x){1 / (1 + exp(-x))}

logit<-function(x) {exp(x)/(1+exp(x))}

ldata$rate = ldata$drift / ldata$initial
lhdataspp$rate = lhdataspp$drift / lhdataspp$initial
lbdataspp$rate = lbdataspp$drift / lbdataspp$initial





####--------------------density--------------------####

### both species combined ###
dldata = ldata[ldata$experiment=='density',]

d0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dldata, family='binomial')

d1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dldata, family='binomial')

AICctab(d0, d1, weights=TRUE, nobs=28)
anova(d0, d1)
precis(d1)
logistic(precis(d1))
post.d1 = sample.naive.posterior(d1)
HPDI(post.d1)
HPDI(logistic(post.d1))
logistic(HPDI(post.d1))

### Hept (h) ###
dlhdataspp = lhdataspp[lhdataspp$experiment=='density',]

d0h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dlhdataspp, family='binomial')

d1h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dlhdataspp, family='binomial')

AICctab(d0h, d1h, weights=TRUE, nobs=28)
anova(d0h, d1h)
precis(d1h)

### Baet (b) ###
dlbdataspp = lbdataspp[lbdataspp$experiment=='density',]

d0b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dlbdataspp, family='binomial')

d1b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dlbdataspp, family='binomial')

AICctab(d0b, d1b, weights=TRUE, nobs=28)
anova(d0b, d1b)
precis(d1b)

##################



####---------canopy + predation + Herbivore density---------####

## both species combinded ##
cpldata = ldata[ldata$experiment=='pred_canopy',]

cp0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=cpldata, family='binomial')

cp1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds, data=cpldata, family='binomial')

cp2 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy, data=cpldata, family='binomial')

cp3 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds, data=cpldata, family='binomial')

cp4 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + canopy:preds, data=cpldata, family='binomial')

cp0i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	initial, data=cpldata, family='binomial')

cp1i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial, data=cpldata, family='binomial')

cp2i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial, data=cpldata, family='binomial')

cp3i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial, data=cpldata, family='binomial')

cp4i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + canopy:preds, data=cpldata, family='binomial')

cp1id = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens, data=cpldata, family='binomial')

cp1id.int = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens + Hdens:preds, data=cpldata, family='binomial')

cp2id = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cpldata, family='binomial')
	
cp2id.int = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cpldata, family='binomial')

cp3id = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens, data=cpldata, family='binomial')

cp3id.int = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens + Hdens:preds,
	 data=cpldata, family='binomial')


plot(rate ~ Hdens, data=cpldata[cpldata$time==2,], pch=preds*2, col=preds+1)
newHdens = 200:1600
Hdens.nopred.pred = logistic(fixef(cp1id)[1] + fixef(cp1id)[2]*2 +fixef(cp1id)[3]*2^2 +
	fixef(cp1id)[4] * 0 + fixef(cp1id)[5] * 212.075 + fixef(cp1id)[6] * newHdens)

Hdens.pred.pred = logistic(fixef(cp1id)[1] + fixef(cp1id)[2]*2 +fixef(cp1id)[3]*2^2 +
	fixef(cp1id)[4] * 1 + fixef(cp1id)[5] * 212.075 + fixef(cp1id)[6] * newHdens)

lines(newHdens, Hdens.nopred.pred, col=1)
lines(newHdens, Hdens.pred.pred, col=2)


AICctab(cp0, cp1, cp2, cp3, cp4, cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE, nobs=40)
AICctab(cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE, nobs=40)
AICctab(cp0i, cp1i, cp2i, cp3i, cp4i, cp1id, cp1id.int, cp2id, cp2id.int,
	cp3id, cp3id.int, weights=TRUE, nobs=40)

anova(cp1id, cp1i)
anova(cp1i, cp3i)
anova(cp1i, cp0i)


## Hept (h) ##
cplhdataspp = lhdataspp[lhdataspp$experiment=='pred_canopy',]

cp0h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=cplhdataspp, family='binomial')

cp1h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds, data=cplhdataspp, family='binomial')

cp2h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy, data=cplhdataspp, family='binomial')

cp3h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds, data=cplhdataspp, family='binomial')

cp4h = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + canopy:preds, data=cplhdataspp, family='binomial')

cp0ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	initial, data=cplhdataspp, family='binomial')

cp1ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial, data=cplhdataspp, family='binomial')

cp2ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial, data=cplhdataspp, family='binomial')

cp3ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial, data=cplhdataspp, family='binomial')

cp4ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + canopy:preds, data=cplhdataspp, family='binomial')

cp1idh = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens, data=cplhdataspp, family='binomial')

cp1id.inth = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens + Hdens:preds, data=cplhdataspp, family='binomial')

cp2idh = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cplhdataspp, family='binomial')
	
cp2id.inth = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cplhdataspp, family='binomial')

cp3idh = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens, data=cplhdataspp, family='binomial')

cp3id.inth = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens + Hdens:preds,
	 data=cplhdataspp, family='binomial')


plot(rate ~ Hdens, data=cplhdataspp[cplhdataspp$time==2,], pch=preds*2, col=preds+1)
newHdens = 200:1600
Hdens.nopred.predh = logistic(fixef(cp1idh)[1] + fixef(cp1idh)[2]*2 +fixef(cp1idh)[3]*2^2 +
	fixef(cp1idh)[4] * 0 + fixef(cp1idh)[5] * 212.075 + fixef(cp1idh)[6] * newHdens)

Hdens.pred.predh = logistic(fixef(cp1idh)[1] + fixef(cp1idh)[2]*2 +fixef(cp1idh)[3]*2^2 +
	fixef(cp1idh)[4] * 1 + fixef(cp1idh)[5] * 212.075 + fixef(cp1idh)[6] * newHdens)

lines(newHdens, Hdens.nopred.predh, col=1)
lines(newHdens, Hdens.pred.predh, col=2)


AICctab(cp0h, cp1h, cp2h, cp3h, cp4h, cp0ih, cp1ih, cp2ih, cp3ih, cp4ih, weights=TRUE, nobs=40)
AICctab(cp0ih, cp1ih, cp2ih, cp3ih, cp4ih, weights=TRUE, nobs=40)
AICctab(cp0ih, cp1ih, cp2ih, cp3ih, cp4ih, cp1idh, cp1id.inth, cp2idh, cp2id.inth,
	cp3idh, cp3id.inth, weights=TRUE, nobs=40)

anova(cp1idh, cp1ih)
anova(cp1ih, cp3ih)
anova(cp1ih, cp0ih)


## Baet (b) ##
cplbdataspp = lbdataspp[lbdataspp$experiment=='pred_canopy',]

cp0b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=cplbdataspp, family='binomial')

cp1b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds, data=cplbdataspp, family='binomial')

cp2b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy, data=cplbdataspp, family='binomial')

cp3b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds, data=cplbdataspp, family='binomial')

cp4b = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + canopy:preds, data=cplbdataspp, family='binomial')

cp0ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	initial, data=cplbdataspp, family='binomial')

cp1ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial, data=cplbdataspp, family='binomial')

cp2ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial, data=cplbdataspp, family='binomial')

cp3ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial, data=cplbdataspp, family='binomial')

cp4ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + canopy:preds, data=cplbdataspp, family='binomial')

cp1idb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens, data=cplbdataspp, family='binomial')

cp1id.intb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial + Hdens + Hdens:preds, data=cplbdataspp, family='binomial')

cp2idb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cplbdataspp, family='binomial')
	
cp2id.intb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial + Hdens, data=cplbdataspp, family='binomial')

cp3idb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens, data=cplbdataspp, family='binomial')

cp3id.intb = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + Hdens + Hdens:preds,
	 data=cplbdataspp, family='binomial')


plot(rate ~ Hdens, data=cplbdataspp[cplbdataspp$time==2,], pch=preds*2, col=preds+1)
newHdens = 200:1600
Hdens.nopred.predb = logistic(fixef(cp1idb)[1] + fixef(cp1idb)[2]*2 +fixef(cp1idb)[3]*2^2 +
	fixef(cp1idb)[4] * 0 + fixef(cp1idb)[5] * 212.075 + fixef(cp1idb)[6] * newHdens)

Hdens.pred.predb = logistic(fixef(cp1idb)[1] + fixef(cp1idb)[2]*2 +fixef(cp1idb)[3]*2^2 +
	fixef(cp1idb)[4] * 1 + fixef(cp1idb)[5] * 212.075 + fixef(cp1idb)[6] * newHdens)

lines(newHdens, Hdens.nopred.predb, col=1)
lines(newHdens, Hdens.pred.predb, col=2)


AICctab(cp0b, cp1b, cp2b, cp3b, cp4b, cp0ib, cp1ib, cp2ib, cp3ib, cp4ib, weights=TRUE, nobs=40)
AICctab(cp0ib, cp1ib, cp2ib, cp3ib, cp4ib, weights=TRUE, nobs=40)
AICctab(cp0ib, cp1ib, cp2ib, cp3ib, cp4ib, cp1idb, cp1id.intb, cp2idb, cp2id.intb,
	cp3idb, cp3id.intb, weights=TRUE, nobs=40)

anova(cp1idb, cp1ib)
anova(cp1ib, cp3ib)
anova(cp1ib, cp0ib)


#############################







####---------------excess food------------------------####
# no food, food linear, food curved
# initial vs. no initial
# food + initial vs. foodperinitial


### for both species combined ###
fldata = ldata[ldata$experiment=='excess_food',]

f0 = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2), data=fldata, family='binomial')

f1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food, data=fldata, family='binomial')

f2 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + I(food^2), data=fldata, family='binomial')

f0i = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2) + initial, data=fldata, family='binomial')

f1i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + initial, data=fldata, family='binomial')

f2i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + I(food^2) + initial, data=fldata, family='binomial')

# f3 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	# foodperinitial, data=fldata, family='binomial')

# f4 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	# foodperinitial + I(foodperinitial^2), 
	# data=fldata, family='binomial')


AICctab(f0, f1, f2, f0i, f1i, f2i, weights=TRUE, nobs=nrow(fldata)/4)
AICctab(f0i, f1i, f2i, weights=TRUE, nobs=nrow(fldata)/4)


anova(f2i, f1i)

curve(logistic(-8.69862 + x * -1.05402 + x^2 * 0.1055), 0, 6)

ranef(f2i)
dotplot(ranef(f2i))

y = logistic(coef(f2i)[[1]][1])
x = 69:98
plot(y[,1] ~ x)
abline(logistic(fixef(f2i)[1]), 0, col=2, lty=2)


### Hept (h) ###
flhdataspp = lhdataspp[lhdataspp$experiment=='excess_food',]

f0ih = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2) + initial, data=flhdataspp, family='binomial')

f1ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + initial, data=flhdataspp, family='binomial')

f2ih = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + I(food^2) + initial, data=flhdataspp, family='binomial')

AICctab(f0ih, f1ih, f2ih, weights=TRUE, nobs=nrow(flhdataspp)/4)
anova(f2ih, f1ih)

### Baet (b) ###
flbdataspp = lbdataspp[lbdataspp$experiment=='excess_food',]

f0ib = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2) + initial, data=flbdataspp, family='binomial')

f1ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + initial, data=flbdataspp, family='binomial')

f2ib = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + I(food^2) + initial, data=flbdataspp, family='binomial')

AICctab(f0ib, f1ib, f2ib, weights=TRUE, nobs=nrow(flbdataspp)/4)
anova(f2ib, f1ib)