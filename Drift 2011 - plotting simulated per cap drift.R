#### plotting per capita drift rates for varied densities and food abundances
#21 Sep 2011
library(lattice)
d = expand.grid(N = seq(1,113,by=2), f = seq(0,8, length=15))

m = .001
b = 1

d$NN = m

d$LRpgd = m * d$N / (b + d$f)
d$LRtot = m * d$N / (b + d$f) * d$N

trellis.par.set("axis.line", list(col="transparent"))
par(mfrow=c(1,2))
print(wireframe(d$LRpgd ~ d$N + d$f,
	main = "Hypothesis 6",
	zlab = list(label = "Per capita drift rate", rot = 90),
	xlab = 'Density',
	ylab = 'Food abundance',
	screen = c(z=-245, x=-75),
	scales = list(arrows=FALSE, cex=.6, col=1, tck=0.75)

))

print(wireframe(d$LRtot ~ d$N + d$f,
	main = "Hypothesis 6",
	zlab = list(label = "Population drift rate", rot = 90),
	xlab = 'Density',
	ylab = 'Food abundance',
	screen = c(z=-245, x=-75),
	scales = list(arrows=FALSE, cex=.6, col=1, tck=0.75)
))
