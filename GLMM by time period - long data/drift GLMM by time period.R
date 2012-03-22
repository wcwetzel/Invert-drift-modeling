## long format data binomial GLMM
# 8 Feb 2012
# 1st run 'drift data reshape.R' to put the data
# into the long format, in which each row is a spot check

library(lme4)
library(ggplot2)
library(rethinking)
library(bbmle)

# logistic transform:
logistic = function(x){1 / (1 + exp(-x))}

logit<-function(x) {exp(x)/(1+exp(x))}

ldata$rate = ldata$drift / ldata$initial

#### do I include a term for time or not? ####
ggplot( data=ldata, aes(y=ldata$drift, x=ldata$time)) +
	geom_point()+
	stat_smooth()
# YES, use a polynomial. time + time^2
# take it for granted, don't test its significance


#### do I include velocity or not? ####
ggplot( data=ldata[ldata$time==2,], aes(y=drift, x=velocity)) +
	geom_point()+
	stat_smooth()
# NO, looks like no relationship. I think it makes sense
# to just leave it out


#### density ####
dldata = ldata[ldata$experiment=='density',]

d0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dldata, family='binomial')

d1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dldata, family='binomial')

AIC(d0, d1)
AICctab(d0, d1, weights=TRUE, nobs=28)
anova(d0, d1)
precis(d1)

##################
#### looking at time effect ####

plot(drift/initial ~ time, data=dldata, xlim=c(0, 5))
newt = seq(0,5, length=100)
newtd = logistic(fixef(d0)[1] + fixef(d0)[2] * newt + fixef(d0)[3] * newt^2)
points(newtd ~ newt, type='l')

##################

#### canopy + predation + Herbivore density ####
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

#############################

#### excess food ####
# no food, food linear, food curved
# initial vs. no initial
# food + initial vs. foodperinitial

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
