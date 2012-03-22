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
### without Hdens ###
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
# run code for excess food first
boxplot(1:4 ~ as.factor(1:4),
	data=cpldata[cpldata$time==2,], border='white',
	las=1, xlab='', ylab='Proportion drifting',
	names=c('No predators', 'Predators', 'Excess food', 'Natural food'),
	ylim=c(0,0.6))
pty = 20
ps = 0.9
points(new.preds+1, mu.cp, pch=pty, cex=ps)
arrows(1, mu.ci.cp[1,1], 1, mu.ci.cp[2,1], length=0.05, 
	angle=90, code=3)
arrows(2, mu.ci.cp[1,2], 2, mu.ci.cp[2,2], length=0.05, 
	angle=90, code=3)
# add excess food points
points(c(4,3), mu.f2i[new.food==6 | new.food==1], pch=pty, cex=ps)
arrows(3, mu.ci.f2i[1,new.food==6], 3, mu.ci.f2i[2,new.food==6],
	length=0.05, angle=90, code=3)


### with Hdens ###
post.cp1id = sample.naive.posterior(cp1id)
newHdens = 200:1600

Hdens.mu.nopred = sapply(newHdens,
	function(z) logistic( mean(post.cp1id[,1] + post.cp1id[,2] * 
	2 +	post.cp1id[,3] * 4 + post.cp1id[,4] * 0 +
	post.cp1id[,5] * 212.075 + post.cp1id[,6] * z))
)

Hdens.mu.pred = sapply(newHdens,
	function(z) logistic( mean(post.cp1id[,1] + post.cp1id[,2] * 
	2 +	post.cp1id[,3] * 2^2 + post.cp1id[,4] * 1 +
	post.cp1id[,5] * 212.075 + post.cp1id[,6] * z))
)



Hdens.ci.nopred = sapply(newHdens,
	function(z) logistic( HPDI(post.cp1id[,1] + post.cp1id[,2] * 
	2 +	post.cp1id[,3] * 2^2 + post.cp1id[,4] * 0 +
	post.cp1id[,5] * 212.075 + post.cp1id[,6] * z))
)

Hdens.ci.pred = sapply(newHdens,
	function(z) logistic( HPDI(post.cp1id[,1] + post.cp1id[,2] * 
	2 +	post.cp1id[,3] * 2^2 + post.cp1id[,4] * 1 +
	post.cp1id[,5] * 212.075 + post.cp1id[,6] * z))
)

plot(rate ~ Hdens, data=cpldata[cpldata$time==2,], pch=preds*2, col=preds+1, ylim=c(0,1))
lines(newHdens, Hdens.mu.nopred)
lines(newHdens, Hdens.ci.nopred[1,], lty=2)
lines(newHdens, Hdens.ci.nopred[2,], lty=2)

lines(newHdens, Hdens.mu.pred, col=2)
lines(newHdens, Hdens.ci.pred[1,], lty=2, col=2)
lines(newHdens, Hdens.ci.pred[2,], lty=2, col=2)


##################################################

#### plotting excess food ####

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
	ylab='Proportion drifting')
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
AICtab(mafdm, mafdm0)
anova(mafdm, mafdm0)
pmafdm = profile(mafdm2)
plot(pmafdm)
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
