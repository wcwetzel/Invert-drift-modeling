# by species
# Hept and Baet separate
# reshaping bruce's data from wide to long
# I want each row to be a spot check 
# (spot check = counts from one time period)
# 8 Feb 2012

#------- read data ----------------------------
dataspp = read.csv(
"/Users/will/Documents/Analysis for colleagues/bruce/Drift data/drift experiments with time by species.csv")

#------- calculate number in channel ---------

## for Hept (h)
dataspp$N1h = dataspp$N0h - dataspp$t1h
dataspp$N2h = dataspp$N1h - dataspp$t2h
dataspp$N3h = dataspp$N2h - dataspp$t3h
dataspp$N4h = dataspp$N3h - dataspp$t4h
dataspp$N0th = dataspp$N0

## for Baet (b)
dataspp$N1b = dataspp$N0b - dataspp$t1b
dataspp$N2b = dataspp$N1b - dataspp$t2b
dataspp$N3b = dataspp$N2b - dataspp$t3b
dataspp$N4b = dataspp$N3b - dataspp$t4b
dataspp$N0tb = dataspp$N0b



# first get rid of block, 
# because i don't much care about it
dataspp = dataspp[,-c(2)]




#------------------- reshape data to long format ----------------

## for both species
ldataspp = reshape(dataspp, varying=list(c('t1h','t2h','t3h','t4h'),
	c('N0th', 'N1h', 'N2h', 'N3h'), c('N1h', 'N2h', 'N3h', 'N4h'), 
	c('t1b','t2b','t3b','t4b'), c('N0tb', 'N1b', 'N2b', 'N3b'), 
	c('N1b', 'N2b', 'N3b', 'N4b')),
	v.names=c('drifth', 'initialh', 'stayh', 'driftb', 'initialb', 'stayh'),
	timevar='time', direction='long', idvar='rep')
# ldataspp$foodperinitial = ldataspp$food / ldataspp$initial
# each row is the count from one time one channel was checked
# initial is the number in the channel at the beginning of
# that time period.
# drift is number drifted
# stay is number stayed

## for Hept (h)
lhdataspp = reshape(dataspp, varying=list(c('t1h','t2h','t3h','t4h'),
	c('N0th', 'N1h', 'N2h', 'N3h'), c('N1h', 'N2h', 'N3h', 'N4h')),
	v.names=c('drift', 'initial', 'stay'), timevar='time', 
	direction='long', idvar='rep')
lhdataspp$foodperinitial = lhdataspp$food / lhdataspp$initial


## for Baet (b)
lbdataspp = reshape(dataspp, varying=list(c('t1b','t2b','t3b','t4b'),
	c('N0tb', 'N1b', 'N2b', 'N3b'), c('N1b', 'N2b', 'N3b', 'N4b')),
	v.names=c('drift', 'initial', 'stay'), timevar='time', 
	direction='long', idvar='rep')
lbdataspp$foodperinitial = lbdataspp$food / lbdataspp$initial
