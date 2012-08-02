## long format data binomial GLMM
# 8 Feb 2012
# 1st run 'drift data reshape.R' to put the data
# into the long format, in which each row is a spot check

library(lme4)
library(bbmle)

# logistic transform:
logistic = function(x){1 / (1 + exp(-x))}


ldata$rate = ldata$drift / ldata$initial


#--------------------- density ---------------------------#
# data for density
dldata = ldata[ldata$experiment=='density',]

# fit models
d0 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2),
	data=dldata, family='binomial')

d1 = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2)
	+ initial, data=dldata, family='binomial')

# compare models, test specification of number of obs
AICctab(d0, d1, weights=TRUE, nobs=14)
AICctab(d0, d1, weights=TRUE, nobs=28)
AICctab(d0, d1, weights=TRUE, nobs=100000)
AICtab(d0, d1, weights=TRUE)

############################################################


#--------------------- canopy + predation ---------------#
# data for canopy + predation
cpldata = ldata[ldata$experiment=='pred_canopy',]

# fit models
cp0i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	initial, data=cpldata, family='binomial')

cp1i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	preds + initial, data=cpldata, family='binomial')

cp2i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + initial, data=cpldata, family='binomial')

cp3i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial, data=cpldata, family='binomial')

cp4i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	canopy + preds + initial + canopy:preds, data=cpldata, 
	family='binomial')

# compare models, test specification of number of obs
AICctab(cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE, nobs=20)
AICctab(cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE, nobs=40)
AICctab(cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE, nobs=10000)
AICtab(cp0i, cp1i, cp2i, cp3i, cp4i, weights=TRUE)

##############################################################


#---------------------- excess food -------------------------#
# data for food
fldata = ldata[ldata$experiment=='excess_food',]

# fit models
f0i = lmer(cbind(drift, stay) ~ (1|rep) + time +
	 I(time^2) + initial, data=fldata, family='binomial')

f1i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + initial, data=fldata, family='binomial')

f2i = lmer(cbind(drift, stay) ~ (1|rep) + time + I(time^2) +
	food + I(food^2) + initial, data=fldata, family='binomial')

# compare models, test specification of number of obs
AICctab(f0i, f1i, f2i, weights=TRUE, nobs=15)
AICctab(f0i, f1i, f2i, weights=TRUE, nobs=30)
AICctab(f0i, f1i, f2i, weights=TRUE, nobs=1000000)
AICtab(f0i, f1i, f2i, weights=TRUE)


