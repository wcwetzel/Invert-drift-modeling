#### plotting drift model predictions
# 18 Dec 2011

# first run 'prob models for drift.R'
# second run 'MLE for drift data by time period.R'
# third run 'simulation of drift model predictions by time period.R'


plot(N0 ~ food, data=data)

par(mfrow=c(2,2))
plot(rate ~ food, data=data)
points(ratem6k ~ food, data=pdata[pdata$N0==213,], type='l', col='green')
points(ratem6 ~ food, data=pdata[pdata$N0==213,], type='l', col='red')
points(ratem5 ~ food, data=pdata[pdata$N0==213,], type='l', col='blue')
points(ratem4 ~ food, data=pdata[pdata$N0==213,], type='l', col='orange')
points(ratem3 ~ food, data=pdata[pdata$N0==213,], type='l', col='red', lty=2)
points(ratem2 ~ food, data=pdata[pdata$N0==213,], type='l', col='blue', lty=2)
points(ratem1 ~ food, data=pdata[pdata$N0==213,], type='l', col='orange', lty=2)

plot(rate ~ food, data=data)
points(ratem6 ~ food, data=pdata[pdata$N0==113,], type='l', col='red')
points(ratem6 ~ food, data=pdata[pdata$N0==213,], type='l', col='red')
points(ratem6 ~ food, data=pdata[pdata$N0==413,], type='l', col='red')
#points(ratem3 ~ food, data=pdata[pdata$N0==113,], type='l', col='red', lty=2)
 points(ratem6k ~ food, data=pdata[pdata$N0==113,], type='l', col='green')
 points(ratem6k ~ food, data=pdata[pdata$N0==213,], type='l', col='green')
 points(ratem6k ~ food, data=pdata[pdata$N0==413,], type='l', col='green')

plot(rate ~ N0, data=data)
points(ratem6k ~ N0, data=pdata[pdata$food==1,], type='l', col='green')
points(ratem6 ~ N0, data=pdata[pdata$food==1,], type='l', col='red')
points(ratem5 ~ N0, data=pdata[pdata$food==1,], type='l', col='blue')
points(ratem4 ~ N0, data=pdata[pdata$food==1,], type='l', col='orange')
points(ratem3 ~ N0, data=pdata[pdata$food==1,], type='l', col='red', lty=2)
points(ratem2 ~ N0, data=pdata[pdata$food==1,], type='l', col='blue', lty=2)
points(ratem1 ~ N0, data=pdata[pdata$food==1,], type='l', col='orange', lty=2)

plot(rate ~ N0, data=data)
points(ratem6 ~ N0, data=pdata[pdata$food==1,], type='l', col='red')
points(ratem6 ~ N0, data=pdata[pdata$food==3,], type='l', col='red')
points(ratem6 ~ N0, data=pdata[pdata$food==6,], type='l', col='red')


points(ratem6k ~ N0, data=pdata[pdata$food==1,], type='l', col='green')
points(ratem6k ~ N0, data=pdata[pdata$food==3,], type='l', col='green')
points(ratem6k ~ N0, data=pdata[pdata$food==6,], type='l', col='green')

par(mfrow=c(1,1))
trellis.par.set("axis.line", list(col="transparent"))
print( wireframe(ratem6 ~ N0 + food, data=pdata,
				drape=TRUE,
				layout=c(0,1),
				xlab = list('Initial density', rot=-48),
				ylab='Food',
				zlab=list('Drift rate', rot=90),
				col.regions = rev( gray(seq(0, 1, length=20))),
				at = seq(0, 0.253, length=20),
				aspect=c(1,1),
				screen = list(z = -245, x=-75),
				scales = list(arrows=FALSE,
				cex=0.5, col='black',
				font=3, tck=1),
				colorkey=FALSE
	))

