### plotting drift binomial GLMMs
# 26 Feb 2012


# 1st run 'drift data reshape.R'
# 2nd run 'drift GLMM by time period.R' to fit models
# 3rd plot models below...

library(rethinking)
library(ggplot2)
library(bbmle)

#### ploting density experiment ####
post.d1 = sample.naive.posterior(d1) # sample posterior of d1
# add random effect noise to posterior table
post.d1$ran = rnorm(nrow(post.d1), mean=0, sd=0.43574)

new.initial = seq(0, 606, by=0.5) # fake density values

## Plot just for time period 2
## first ignoring random effect
# calcuate mu
# calculate predicted value for ALL posterior sampled parms
# for each and every fake density value using sapply
# then take the mean across sampled parms for each fake density
mu = sapply(new.initial, 
	function(z) logistic(mean(post.d1[,1] + post.d1[,2] *2 + 
	post.d1[,3] * 2^2 + post.d1[,4] * z)))

# calculate ci by taking middle 95% of posterior
mu.ci = sapply(new.initial, 
	function(z) logistic(HPDI(post.d1[,1] + post.d1[,2] *2 + 
	post.d1[,3] * 2^2 + post.d1[,4] * z)))

# plot data, predictions, and CIs
plot(drift/initial ~ initial, 
	data=dldata[dldata$time==2,], pch=20,
	las=1, xlab='Density', ylab='Proportion drifting')
lines(new.initial, mu)
lines(new.initial, mu.ci[1,], lty=2)
lines(new.initial, mu.ci[2,], lty=2)


## now with random effect
mu.ran = sapply(new.initial, 
	function(z) logistic(mean(post.d1[,1] + post.d1[,2] *2 + 
	post.d1[,3] * 2^2 + post.d1[,4] * z + post.d1[,5])))

# calculate ci by taking middle 95% of posterior
mu.ci.ran = sapply(new.initial, 
	function(z) logistic(HPDI(post.d1[,1] + post.d1[,2] *2 + 
	post.d1[,3] * 2^2 + post.d1[,4] * z + post.d1[,5])))

# plot data, predictions, and CIs
plot(drift/initial ~ initial, 
	data=dldata[dldata$time==2,], pch=20,
	las=1, xlab='Density', ylab='Proportion drifting')
lines(new.initial, mu.ran)
lines(new.initial, mu.ci.ran[1,], lty=2)
lines(new.initial, mu.ci.ran[2,], lty=2)
# I don't think we should do it with the random effect because
# we are really just interested in what happens within a channel
# how the density affects things in a channel




#######################################

#### plotting canopy + predation ####
# I'm just going to do preds vs. no preds
# No model averaging. The AIC values of canopy vs. 
# no canopy were close, but the coef for canopy was
# so tiny as to be ecologically unimportant.

# I'm going to do this for the mean initial for time
# period 2 (212.075) at time period 2

post.cp1i = sample.naive.posterior(cp1i)

new.preds = 0:1
# new.canopy = seq(0, 100, by = 0.1) not going to plot canopy

# calculate mean for preds
mu.cp = sapply(new.preds, 
	function(z) logistic(mean(post.cp1i[,1] + post.cp1i[,2] * 
	2 +	post.cp1i[,3] * 2^2 + post.cp1i[,4] * z +
	post.cp1i[,5] * 212.075)))

# calculate ci by taking middle 95% of posterior
mu.ci.cp = sapply(new.preds, 
	function(z) logistic(HPDI(post.cp1i[,1] + post.cp1i[,2] * 
	2 +	post.cp1i[,3] * 2^2 + post.cp1i[,4] * z +
	post.cp1i[,5] * 212.075)))

# plot data, predictions, and CIs
boxplot(drift/initial ~ as.factor(preds),
	data=cpldata[cpldata$time==2,], border='white',
	las=1, xlab='', ylab='Proportion drifting',
	names=c('No predators', 'Predators'),
	ylim=c(0,0.6))
points(new.preds+1, mu.cp, pch=5)
arrows(1, mu.ci.cp[1,1], 1, mu.ci.cp[2,1], length=0.05, 
	angle=90, code=3)
arrows(2, mu.ci.cp[1,2], 2, mu.ci.cp[2,2], length=0.05, 
	angle=90, code=3)

##################################################

#### plotting excess food ####

# I'm going to do this for the mean initial for time
# period 2 (211.8333) at time period 2

post.f2i = sample.naive.posterior(f2i) # sample posterior of d1
# add random effect noise to posterior table

new.food = seq(0, 6, by=0.01) # fake density values


## Plot just for time period 2

# calcuate mu
mu.f2i = sapply(new.food, 
	function(z) logistic(mean(post.f2i[,1] + post.f2i[,2] * 2 + 
	post.f2i[,3] * 2^2 + post.f2i[,4] * z + 
	post.f2i[,5] * z^2 + post.f2i[,6] * 211.8333)))

# calculate ci by taking middle 95% of posterior
mu.ci.f2i = sapply(new.food, 
	function(z) logistic(HPDI(post.f2i[,1] + post.f2i[,2] * 2 + 
	post.f2i[,3] * 2^2 + post.f2i[,4] * z + 
	post.f2i[,5] * z^2 + post.f2i[,6] * 211.8333)))

# plot data, predictions, and CIs
plot(drift/initial ~ food, 
	data=fldata[fldata$time==2,], pch=20,
	las=1, xlab='Food isolation days', 
	ylab='Proportion drifting')
lines(new.food, mu.f2i)
lines(new.food, mu.ci.f2i[1,], lty=2)
lines(new.food, mu.ci.f2i[2,], lty=2)



