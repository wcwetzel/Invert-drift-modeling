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
ggplot( data=ldata, aes(y=ldata$drift, x=ldata$velocity)) +
	geom_point()+
	stat_smooth()
# NO, looks like no relationship. I think it makes sense
# to just leave it out


#### density ####
dldata = ldata[ldata$experiment=='density',]

m0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dldata, family='binomial')

m1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dldata, family='binomial')

AIC(m0, m1)
AICctab(m0,m1, nobs=28)
anova(m0, m1)


#### canopy + predation ####
cpldata = ldata[ldata$experiment=='pred_canopy',]

m0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) + 
	velocity, data=cpldata, family='binomial')

m1 = lmer(cbind(drift, stay) ~ (1|rep) + preds + velocity + 
	time + I(time^2), data=cpldata, family='binomial')

m2 = lmer(cbind(drift, stay) ~ (1|rep) + canopy + 
	velocity + time + I(time^2), data=cpldata, family='binomial')

m3 = lmer(cbind(drift, stay) ~ (1|rep) + canopy + preds + 
	velocity + time + I(time^2), data=cpldata, family='binomial')

m4 = lmer(cbind(drift, stay) ~ (1|rep) + canopy + preds + 
	velocity + time + I(time^2) + initial, data=cpldata, family='binomial')

m4.5 = lmer(cbind(drift, stay) ~ (1|rep) + preds + 
	velocity + time + I(time^2) + initial, data=cpldata, family='binomial')

m4.7 = lmer(cbind(drift, stay) ~ (1|rep) + preds + 
	time + I(time^2) + initial, data=cpldata, family='binomial')


m5 = lmer(cbind(drift, stay) ~ (1|rep) + canopy + preds + 
	velocity + time + I(time^2) + (1 + canopy|preds), 
	data=cpldata, family='binomial')


AIC(m0, m1, m2, m3, m4, m5, m4.5, m4.7)
AICctab(logLik(m0), logLik(m1), logLik(m2), logLik(m3), logLik(m4), logLik(m4.5),
	logLik(m5), logLik(m4.7), weights=TRUE, nobs=nrow(cpldata)/4)

anova(m0, m1, m3)
anova(m2, m3)
anova(m2, m0)
anova(m1,m4)


# excess food
fldata = ldata[ldata$experiment=='excess_food',]

m0 = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2), data=fldata, family='binomial')

m1 = lmer(cbind(drift, stay) ~ (1|rep) + food + 
	time + I(time^2), data=fldata, family='binomial')

m2 = lmer(cbind(drift, stay) ~ (1|rep) + food + I(food^2) +
	time + I(time^2), data=fldata, family='binomial')

m3 = lmer(cbind(drift, stay) ~ (1|rep) + food + I(food^2) +
	time + I(time^2) + I(food^3), data=fldata, family='binomial')

m4 = lmer(cbind(drift, stay) ~ (1|rep) + food + I(food^2) +
	time + I(time^2) + velocity, data=fldata, family='binomial')

m4.5 = lmer(cbind(drift, stay) ~ (1|rep) + I(food/initial)  +
	time + I(time^2), data=fldata, family='binomial')

m4.7 = lmer(cbind(drift, stay) ~ (1|rep) + I(food/initial)  +
	time + I(time^2) + initial, data=fldata, family='binomial')

m4.8 = lmer(cbind(drift, stay) ~ (1|rep) + initial +
	time + I(time^2), data=fldata, family='binomial')
	
m4.9 = lmer(cbind(drift, stay) ~ (1|rep) + initial + food + 
	I(food^2) + time + I(time^2), data=fldata, family='binomial')


m5 = lmer(cbind(drift, stay) ~ (1|rep) + I(food/N0) +
	time + I(time^2), data=fldata, family='binomial')


AICctab(logLik(m0), logLik(m1), logLik(m2), logLik(m3), logLik(m4), logLik(m4.5),
	logLik(m5), logLik(m4.7), logLik(m4.8), logLik(m4.9), 
	weights=TRUE, nobs=nrow(fldata)/4)

anova(m0, m1, m2, m3)

curve(logistic(-8.69862 + x * -1.05402 + x^2 * 0.1055), 0, 6)

ranef(m2)
dotplot(ranef(m2))

y = logistic(coef(m2)[[1]][1])
x = 69:98
plot(y[,1] ~ x)
abline(logistic(fixef(m2)[1]), 0, col=2, lty=2)
