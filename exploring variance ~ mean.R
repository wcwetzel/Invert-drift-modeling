### DRIFT
### 12 Sep 2011
# exploring the relationship between the variance and the mean

V0 = 56.9
N0 = 113
b = 0.1
t = 10
f = 1
N = 50
c = 0.4

curve(x/2 + (V0 - N0/2) * exp(-2 * b * t / f), 0, 113)


curve(x/2 + c * exp(-2 * b * t / f), 0 , 113)