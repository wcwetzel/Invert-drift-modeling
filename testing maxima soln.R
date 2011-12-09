# testing maxima's solution for model4's N(t)

mine = function(t, N0=100, m=0.1, a=0.1){
	m / ((m/N0 + a)*exp(t*m)-a)
}

maxima = function(t, N0=100, m=0.1, a=0.1){
	(m * N0) / ((m + N0 * a) * exp(m * t) - N0 * a)
}

t = seq(0,10, by=0.1)


plot(t, mine(t), type='l')
points(t, maxima(t), type='l', col='red')