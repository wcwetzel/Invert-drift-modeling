### plotting drift binomial GLMMs
# 26 Feb 2012


# 1st run 'drift data reshape.R'
# 2nd run 'drift GLMM by time period.R' to fit models
# 3rd plot models below...

library(rethinking)
library(ggplot2)
library(bbmle)




#------------plotting canopy + predation---------------------####
# model 3 canopy + predation
# I'm going to do this for the mean initial for time
# period 2 (212.075) at time period 2

### for both species combined ###
post.cp3i = sample.naive.posterior(cp3i)

new.canopy = seq(0, 100, by = 0.1)

# calculate mean for preds
mu.cp3i.preds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3i[,1] + post.cp3i[,2] * 
	2 +	post.cp3i[,3] * 2^2 + post.cp3i[,4] * z + post.cp3i[,5] * 1 +
	post.cp3i[,6] * 212.075)))

# calculate ci for preds
ci.cp3i.preds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3i[,1] + post.cp3i[,2] * 
	2 +	post.cp3i[,3] * 2^2 + post.cp3i[,4] * z + post.cp3i[,5] * 1 +
	post.cp3i[,6] * 212.075)))

# calculate mean for no preds
mu.cp3i.nopreds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3i[,1] + post.cp3i[,2] * 
	2 +	post.cp3i[,3] * 2^2 + post.cp3i[,4] * z + post.cp3i[,5] * 0 +
	post.cp3i[,6] * 212.075)))

# calculate ci for no preds
ci.cp3i.nopreds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3i[,1] + post.cp3i[,2] * 
	2 +	post.cp3i[,3] * 2^2 + post.cp3i[,4] * z + post.cp3i[,5] * 0 +
	post.cp3i[,6] * 212.075)))

plot(drift/initial ~ canopy, data=cpldata[cpldata$time==2 & cpldata$preds==1,],
	pch=2, col='red', ylim=c(0,0.7), main='both spp')
points(rate ~ canopy, data=cpldata[cpldata$time==2 & cpldata$preds==0,],
	pch=1, col='black')
lines(new.canopy, mu.cp3i.preds, col='red')
lines(new.canopy, ci.cp3i.preds[1,], col='red', lty=2)
lines(new.canopy, ci.cp3i.preds[2,], col='red', lty=2)
lines(new.canopy, mu.cp3i.nopreds)
lines(new.canopy, ci.cp3i.nopreds[1,], lty=2)
lines(new.canopy, ci.cp3i.nopreds[2,], lty=2)


### for Hept (h) ###
# time period 2
# mean initial Hept at time period 2 = 88.75

post.cp3ih = sample.naive.posterior(cp3ih)

new.canopy = seq(0, 100, by = 0.1)

# calculate mean for preds
mu.cp3ih.preds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3ih[,1] + post.cp3ih[,2] * 
	2 +	post.cp3ih[,3] * 2^2 + post.cp3ih[,4] * z + post.cp3ih[,5] * 1 +
	post.cp3ih[,6] * 88.75)))

# calculate ci for preds
ci.cp3ih.preds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3ih[,1] + post.cp3ih[,2] * 
	2 +	post.cp3ih[,3] * 2^2 + post.cp3ih[,4] * z + post.cp3ih[,5] * 1 +
	post.cp3ih[,6] * 88.75)))

# calculate mean for no preds
mu.cp3ih.nopreds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3ih[,1] + post.cp3ih[,2] * 
	2 +	post.cp3ih[,3] * 2^2 + post.cp3ih[,4] * z + post.cp3ih[,5] * 0 +
	post.cp3ih[,6] * 88.75)))

# calculate ci for no preds
ci.cp3ih.nopreds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3ih[,1] + post.cp3ih[,2] * 
	2 +	post.cp3ih[,3] * 2^2 + post.cp3ih[,4] * z + post.cp3ih[,5] * 0 +
	post.cp3ih[,6] * 88.75)))

plot(drift/initial ~ canopy, data=cplhdataspp[cplhdataspp$time==2 & cplhdataspp$preds==1,],
	pch=2, col='red', ylim=c(0,0.7), main='Hept')
points(rate ~ canopy, data=cplhdataspp[cplhdataspp$time==2 & cplhdataspp$preds==0,],
	pch=1, col='black')
lines(new.canopy, mu.cp3ih.preds, col='red')
lines(new.canopy, ci.cp3ih.preds[1,], col='red', lty=2)
lines(new.canopy, ci.cp3ih.preds[2,], col='red', lty=2)
lines(new.canopy, mu.cp3ih.nopreds)
lines(new.canopy, ci.cp3ih.nopreds[1,], lty=2)
lines(new.canopy, ci.cp3ih.nopreds[2,], lty=2)

### for Baet (b) ###
# time period 2
# mean initial Baet at time period 2 = 123.325
post.cp3ib = sample.naive.posterior(cp3ib)

new.canopy = seq(0, 100, by = 0.1)

# calculate mean for preds
mu.cp3ib.preds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3ib[,1] + post.cp3ib[,2] * 
	2 +	post.cp3ib[,3] * 2^2 + post.cp3ib[,4] * z + post.cp3ib[,5] * 1 +
	post.cp3ib[,6] * 123.325)))

# calculate ci for preds
ci.cp3ib.preds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3ib[,1] + post.cp3ib[,2] * 
	2 +	post.cp3ib[,3] * 2^2 + post.cp3ib[,4] * z + post.cp3ib[,5] * 1 +
	post.cp3ib[,6] * 123.325)))

# calculate mean for no preds
mu.cp3ib.nopreds = sapply(new.canopy, 
	function(z) mean(logistic(post.cp3ib[,1] + post.cp3ib[,2] * 
	2 +	post.cp3ib[,3] * 2^2 + post.cp3ib[,4] * z + post.cp3ib[,5] * 0 +
	post.cp3ib[,6] * 123.325)))

# calculate ci for no preds
ci.cp3ib.nopreds = sapply(new.canopy, 
	function(z) HPDI(logistic(post.cp3ib[,1] + post.cp3ib[,2] * 
	2 +	post.cp3ib[,3] * 2^2 + post.cp3ib[,4] * z + post.cp3ib[,5] * 0 +
	post.cp3ib[,6] * 123.325)))

plot(drift/initial ~ canopy, data=cplbdataspp[cplbdataspp$time==2 & cplbdataspp$preds==1,],
	pch=2, col='red', ylim=c(0,0.7), main='Baet')
points(rate ~ canopy, data=cplbdataspp[cplhdataspp$time==2 & cplbdataspp$preds==0,],
	pch=1, col='black')
lines(new.canopy, mu.cp3ib.preds, col='red')
lines(new.canopy, ci.cp3ib.preds[1,], col='red', lty=2)
lines(new.canopy, ci.cp3ib.preds[2,], col='red', lty=2)
lines(new.canopy, mu.cp3ib.nopreds)
lines(new.canopy, ci.cp3ib.nopreds[1,], lty=2)
lines(new.canopy, ci.cp3ib.nopreds[2,], lty=2)





##################################################

#-------------------plotting excess food--------------------------#

#### for both spp combined ####
# I'm going to do this for the mean initial for time
# period 2 (211.8333) at time period 2

post.f2i = sample.naive.posterior(f2i) # sample posterior of d1
# is the posterior approximately MVN?
# parm values are not near the boundaries except for intercept
# which is close to zero logistic(-14)=7e-7 and time which is
# close to 1, logistic(6.3)=0.998

new.food = seq(1, 6, by=0.01) # fake density values


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
quartz(width=1.2*3.5, height=1*3.5)
par(mar=c(5,4,4,5)+0.1)
plot(drift/initial ~ I(food-1), 
	data=fldata[fldata$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='Proportion drifting', main='both spp')
lines(new.food-1, mu.f2i)
lines(new.food-1, mu.ci.f2i[1,], lty=2)
lines(new.food-1, mu.ci.f2i[2,], lty=2)

par(new=TRUE)
plot(afdm ~ I(food-1), data=fldata[fldata$time==2,],
	col='white', xaxt='n', yaxt='n', xlab='', ylab='')
axis(4, las=1)
mtext(expression(paste('Remaining biofilm (g m'^{-2},')')), side=4, line=3)
mafdm = lm(afdm ~ food, data=data[data$experiment=='excess_food',])
mafdm0 = lm(afdm ~ 1, data=data[data$experiment=='excess_food',])
mafdm2 = mle2(afdm ~ dnorm(mean = a + b * food, sd=s), start=list(a=1, b=0, s=2), 
	data=data[data$experiment=='excess_food',])
#AICtab(mafdm, mafdm0)
#anova(mafdm, mafdm0)
#pmafdm = profile(mafdm2)
#plot(pmafdm)
post.mafdm = sample.naive.posterior(mafdm)

# calcuate mu
mu.afdm = sapply(new.food, 
	function(z) mean(post.mafdm[,1] + post.mafdm[,2] * z))

# calculate ci by taking middle 95% of posterior
ci.afdm = sapply(new.food, 
	function(z) HPDI(post.mafdm[,1] + post.mafdm[,2] * z))

lines(new.food-1, mu.afdm, col='grey')
lines(new.food-1, ci.afdm[1,], lty=2, col='grey')
lines(new.food-1, ci.afdm[2,], lty=2, col='grey')



#### for Hept (h) ####
# I'm going to do this for the mean initial for time
# period 2 (88.8) at time period 2

post.f2ih = sample.naive.posterior(f2ih) # sample posterior of d1
# is the posterior approximately MVN?
# parm values are not near the boundaries except for intercept
# which is close to zero logistic(-14)=7e-7 and time which is
# close to 1, logistic(6.3)=0.998

new.food = seq(1, 6, by=0.01) # fake density values


## Plot just for time period 2

# calcuate mu
mu.f2ih = sapply(new.food, 
	function(z) logistic(mean(post.f2ih[,1] + post.f2ih[,2] * 2 + 
	post.f2ih[,3] * 2^2 + post.f2ih[,4] * z + 
	post.f2ih[,5] * z^2 + post.f2ih[,6] * 88.8)))

# calculate ci by taking middle 95% of posterior
mu.ci.f2ih = sapply(new.food, 
	function(z) logistic(HPDI(post.f2ih[,1] + post.f2ih[,2] * 2 + 
	post.f2ih[,3] * 2^2 + post.f2ih[,4] * z + 
	post.f2ih[,5] * z^2 + post.f2ih[,6] * 88.8)))

# plot data, predictions, and CIs
quartz(width=1.2*3.5, height=1*3.5)
par(mar=c(5,4,4,5)+0.1)
plot(drift/initial ~ I(food-1), 
	data=flhdataspp[flhdataspp$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='Proportion drifting', main='Hept')
lines(new.food-1, mu.f2ih)
lines(new.food-1, mu.ci.f2ih[1,], lty=2)
lines(new.food-1, mu.ci.f2ih[2,], lty=2)


#### for Baet (b) ####
# I'm going to do this for the mean initial for time
# period 2 (123.0333) at time period 2
post.f2ib = sample.naive.posterior(f2ib) 

new.food = seq(1, 6, by=0.01) # fake density values


## Plot just for time period 2

# calcuate mu
mu.f2ib = sapply(new.food, 
	function(z) logistic(mean(post.f2ib[,1] + post.f2ib[,2] * 2 + 
	post.f2ib[,3] * 2^2 + post.f2ib[,4] * z + 
	post.f2ib[,5] * z^2 + post.f2ib[,6] * 123.0333)))

# calculate ci by taking middle 95% of posterior
mu.ci.f2ib = sapply(new.food, 
	function(z) logistic(HPDI(post.f2ib[,1] + post.f2ib[,2] * 2 + 
	post.f2ib[,3] * 2^2 + post.f2ib[,4] * z + 
	post.f2ib[,5] * z^2 + post.f2ib[,6] * 123.0333)))

# plot data, predictions, and CIs

plot(drift/initial ~ I(food-1), 
	data=flbdataspp[flbdataspp$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='Proportion drifting', main='Baet')
lines(new.food-1, mu.f2ib)
lines(new.food-1, mu.ci.f2ib[1,], lty=2)
lines(new.food-1, mu.ci.f2ib[2,], lty=2)

### for f1ib ###
post.f1ib = sample.naive.posterior(f1ib)

# calcuate mu
mu.f1ib = sapply(new.food, 
	function(z) logistic(mean(post.f1ib[,1] + post.f1ib[,2] * 2 + 
	post.f1ib[,3] * 2^2 + post.f1ib[,4] * z + 
	post.f1ib[,5] * 123.0333)))

# calculate ci by taking middle 95% of posterior
mu.ci.f1ib = sapply(new.food, 
	function(z) logistic(HPDI(post.f1ib[,1] + post.f1ib[,2] * 2 + 
	post.f1ib[,3] * 2^2 + post.f1ib[,4] * z + 
	post.f1ib[,5] * 123.0333)))

# plot data, predictions, and CIs

lines(new.food-1, mu.f1ib, col=3)
lines(new.food-1, mu.ci.f1ib[1,], lty=2, col=3)
lines(new.food-1, mu.ci.f1ib[2,], lty=2, col=3)















#----------------------fancier plots----------------------------------#
# plotting combined excess food plus predators
# first with preds on right
par(fig=c(0, 0.6, 0, 1))
plot(drift/initial ~ I(food-1), 
	data=fldata[fldata$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='Proportion drifting', ylim=c(0,0.5))
lines(new.food-1, mu.f2i)
lines(new.food-1, mu.ci.f2i[1,], lty=2)
lines(new.food-1, mu.ci.f2i[2,], lty=2)
# adding AFDM
par(new=TRUE)
plot(afdm ~ I(food-1), data=fldata[fldata$time==2,],
	col='white', xaxt='n', yaxt='n', xlab='', ylab='')
axis(4, las=1)
mtext(expression(paste('Remaining biofilm (g m'^{-2},')')), side=4, line=3)

lines(new.food-1, mu.afdm, col='grey')
lines(new.food-1, ci.afdm[1,], lty=2, col='grey')
lines(new.food-1, ci.afdm[2,], lty=2, col='grey')

# adding the predator vs. no predator part
par(fig=c(0.6,1,0,1), new=TRUE)
boxplot(1:2 ~ as.factor(1:2),
	data=cpldata[cpldata$time==2,], border='white',
	las=1, xlab='', ylab='Proportion drifting',
	names=c('No preds', 'Preds'),
	ylim=c(0,0.5))
pty = 20
ps = 0.9
points(new.preds+1, mu.cp, pch=pty, cex=ps)
arrows(1, mu.ci.cp[1,1], 1, mu.ci.cp[2,1], length=0.05, 
	angle=90, code=3)
arrows(2, mu.ci.cp[1,2], 2, mu.ci.cp[2,2], length=0.05, 
	angle=90, code=3)


# second with preds on left
par(fig=c(0.25, 1, 0, 1))
plot(drift/initial ~ I(food-1), 
	data=fldata[fldata$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='', yaxt='n', ylim=c(0,0.5))
lines(new.food-1, mu.f2i)
lines(new.food-1, mu.ci.f2i[1,], lty=2)
lines(new.food-1, mu.ci.f2i[2,], lty=2)
# adding AFDM
par(new=TRUE)
plot(afdm ~ I(food-1), data=fldata[fldata$time==2,],
	col='white', xaxt='n', yaxt='n', xlab='', ylab='')
axis(4, las=1, col.ticks='gray', col.axis='gray')
mtext(expression(paste('Remaining biofilm (g m'^{-2},')')), side=4, line=3, col='gray')
lines(new.food-1, mu.afdm, col='grey')
lines(new.food-1, ci.afdm[1,], lty=2, col='grey')
lines(new.food-1, ci.afdm[2,], lty=2, col='grey')

# adding the predator vs. no predator part
par(fig=c(0.01,0.5,0,1), new=TRUE)
boxplot(1:2 ~ as.factor(1:2),
	data=cpldata[cpldata$time==2,], border='white',
	las=1, xlab='Predator treatment', ylab='Proportion drifting',
	names=c('No preds', 'Preds'),
	ylim=c(0,0.5))
pty = 20
ps = 0.9
points(new.preds+1, mu.cp, pch=pty, cex=ps)
arrows(1, mu.ci.cp[1,1], 1, mu.ci.cp[2,1], length=0.05, 
	angle=90, code=3)
arrows(2, mu.ci.cp[1,2], 2, mu.ci.cp[2,2], length=0.05, 
	angle=90, code=3)


par(fig=c(0.25,1,0,1))
plot(1)


# third with preds on left and AFDM on top
topofdrift = 0.8
bottomofAFDM = 0.45
leftofAFDM = 0.2325

par(fig=c(leftofAFDM, 1, 0, topofdrift))
plot(drift/initial ~ I(food-1), 
	data=fldata[fldata$time==2,], pch=20,
	las=1, xlab='Substrate isolation (days)', 
	ylab='', yaxt='n', ylim=c(0,0.5))
lines(new.food-1, mu.f2i)
lines(new.food-1, mu.ci.f2i[1,], lty=2)
lines(new.food-1, mu.ci.f2i[2,], lty=2)

# adding AFDM
par(new=TRUE, fig=c(leftofAFDM, 1, bottomofAFDM, 1))
plot(afdm ~ I(food-1), data=fldata[fldata$time==2,],
	pch=20, xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(0,3))
axis(1, labels=FALSE)
axis(4, las=1, at=c(0,1,2,3))
mtext(expression(paste('Remaining biofilm (g m'^{-2},')')), side=4, line=3)
lines(new.food-1, mu.afdm)
lines(new.food-1, ci.afdm[1,], lty=2)
lines(new.food-1, ci.afdm[2,], lty=2)

# adding the predator vs. no predator part
par(fig=c(0.01,0.5,0,topofdrift), new=TRUE)
boxplot(1:2 ~ as.factor(1:2),
	data=cpldata[cpldata$time==2,], border='white',
	las=1, xlab='Predator treatment', ylab='Proportion drifting',
	names=c('No preds', 'Preds'),
	ylim=c(0,0.5))
pty = 20
ps = 0.9
points(new.preds+1, mu.cp, pch=pty, cex=ps)
arrows(1, mu.ci.cp[1,1], 1, mu.ci.cp[2,1], length=0.05, 
	angle=90, code=3)
arrows(2, mu.ci.cp[1,2], 2, mu.ci.cp[2,2], length=0.05, 
	angle=90, code=3)
