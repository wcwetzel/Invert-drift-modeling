### DRIFT
### 16 Sep 2011
# exploring determinisitc models for drift/emigration

require(deSolve)

## fixed number emigrating
# model
fixed = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -E
		res = dn.dt
		list(res)
		})
}

# parms
parms = c( E = 10 )

# define timepoints to be recorded
tspan = 10
tsteps = 500
times = seq(0, tspan, length = tsteps)

# define initial values
xstart = c(n = 113)

# run lsoda
out.fixed = as.data.frame(lsoda(xstart, times, fixed, parms))

## fixed proportion emigrating
# model
prop = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -d * n
		res = dn.dt
		list(res)
		})
}

# parms
parms = c( d = 0.15 )

# define timepoints to be recorded
tspan = 10
tsteps = 500
times = seq(0, tspan, length = tsteps)

# define initial values
xstart = c(n = 113)

# run lsoda
out.prop = as.data.frame(lsoda(xstart, times, prop, parms))

## dd emigrating
# model
dd = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = -d * n^2
		res = dn.dt
		list(res)
		})
}

# parms
parms = c( d = 0.002 )

# define timepoints to be recorded
tspan = 10
tsteps = 500
times = seq(0, tspan, length = tsteps)

# define initial values
xstart = c(n = 113)

# run lsoda
out.dd = as.data.frame(lsoda(xstart, times, dd, parms))

## scaled exponential emigrating
# model
se = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = - exp(s * n) / (10 + exp(s * n)) * n
		res = dn.dt
		list(res)
		})
}

# parms
parms = c( s = 0.001 )

# define timepoints to be recorded
tspan = 10
tsteps = 500
times = seq(0, tspan, length = tsteps)

# define initial values
xstart = c(n = 113)

# run lsoda
out.se = as.data.frame(lsoda(xstart, times, se, parms))

### HOW ABOUT A TWO PART E RATE? ONE CONSTANT PROPORTION. ONE DD.
## two part rate:
# model
tp = function(t, x, parms) {
	with( as.list( c(parms, x)), {
		dn.dt = (-d * n) - E
		res = dn.dt
		list(res)
		})
}

# parms
parms = c( d = 0.1, E = 5 )

# define timepoints to be recorded
tspan = 10
tsteps = 500
times = seq(0, tspan, length = tsteps)

# define initial values
xstart = c(n = 113)

# run lsoda
out.tp = as.data.frame(lsoda(xstart, times, tp, parms))





# plot
plot(out.fixed$n ~ out.fixed$time, type = 'l', ylim=c(0,113))
points(out.prop$n ~ out.prop$time, type = 'l', lty=2)
points(out.dd$n ~ out.dd$time, type = 'l', lty=3, col='red')
points(out.se$n ~ out.se$time, type = 'l', lty=4)
points(out.tp$n ~ out.tp$time, type = 'l', lty=5)

