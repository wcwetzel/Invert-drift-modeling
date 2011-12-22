# what functional form do I want for models that have
# rational decrease with food???
# I want them to asymptote at a low of m.


m = 0.1
a = 0.0001
c = 1
k = 0.1

par(mfrow=c(1,2))

## 1st for when denity matters (model6) ##
N = 600
curve(m + (a * N)/(c + x), 0, 10, ylim=c(0.09, 0.3))

N = 0
curve(m + (a * N)/(c + x), 0, 10, add=TRUE, col='blue')

N = 50
curve(m + (a * N)/(c + x), 0, 10, add=TRUE)

N = 250
curve(m + (a * N)/(c + x), 0, 10, add=TRUE)

N = 100
curve(m + (a * N)/(c + x), 0, 10, add=TRUE)

# with extra parm for food intercept
N = 600
curve(m + (k + a * N)/(c + x), 0, 10, add=TRUE, lty=2)

N = 0
curve(m + (k + a * N)/(c + x), 0, 10, add=TRUE, col='blue', lty=2)

N = 50
curve(m + (k + a * N)/(c + x), 0, 10, add=TRUE, lty=2)

N = 250
curve(m + (k + a * N)/(c + x), 0, 10, add=TRUE, lty=2)

N = 100
curve(m + (k + a * N)/(c + x), 0, 10, add=TRUE, lty=2)

# when density doesn't matter a or N = 0