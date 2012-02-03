# Binomial GLMM!!!!
# All invert DRIFT experiments separately
# canopy+predation, density, food
# 31 Jan 2012

library(bbmle)
library(ggplot2)
library(lme4)

data = read.csv(
"/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time.csv")
data$rate = data$total / data$N0
data$rep = 1:nrow(data)
data$stay = data$N0 - data$total
datacd = data[-which(data$experiment=='excess_food'),]
datac = data[data$experiment=='pred_canopy',]
datad = data[data$experiment=='density',]
dataf = data[data$experiment=='excess_food',]



# logistic transform:
logistic = function(x){
	1 / (1 + exp(-x))
}

# canopy+predation
m0 = lmer(cbind(total, stay) ~ (1|rep) + 1, data=datac, family=binomial)

m1 = lmer(cbind(total, stay) ~ (1|rep) + canopy, data=datac, family=binomial)

m1p = lmer(cbind(total, stay) ~ (1|rep) + canopy + preds, data=datac, family=binomial)

m0p = lmer(cbind(total, stay) ~ (1|rep) + preds, data=datac, family=binomial)

AICtab(m0, m1, m0p, m1p)
anova(m0, m1)
anova(m0, m0p)

# density
m0 = lmer(cbind(total, stay) ~ (1|rep) + 1, data=datad, family=binomial)

m1 = lmer(cbind(total, stay) ~ (1|rep) + N0, data=datad, family=binomial)

AICtab(m0, m1)
anova(m0, m1)

# food


