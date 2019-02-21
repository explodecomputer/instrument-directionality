library(TwoSampleMR)
library(simulateGP)

filename <- commandArgs(T)[1]
outfile <- commandArgs(T)[2]
jid <- as.numeric(commandArgs(T)[3])

load(filename)

set.seed(jid)

sims <- length(l)
for(i in 1:sims)
{
	message(i)
	id <- names(l)[[i]]
	message(id)
	l[[i]]$estimates <- try(test_system(l[[i]], id))
}

save(l, file=filename)
