## long format data binomial GLMM
# 9 Feb 2012
# 1st run 'drift data reshape.R' to put the data
# into the long format, in which each row is a spot check

library(lme4)
library(ggplot2)
# canopy + predation

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
	velocity + time + I(time^2) + (1 + canopy|preds), data=cpldata, family='binomial')


AICctab(logLik(m0), logLik(m1), logLik(m2), logLik(m3), logLik(m4), logLik(m4.5),
	logLik(m5), logLik(m4.7), weights=TRUE, nobs=nrow(cpldata)/4)


# plotting

plot(resid(m4.7) ~ cpldata$canopy)
anova(m4.7)
newdata = data.frame(preds=0, canopy=1:99, time=3, initial=200)
p = predict(m4.7, newdata=newdata, se.fit=TRUE)

ggplot(data = cpldata, )

