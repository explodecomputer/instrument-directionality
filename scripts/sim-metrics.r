library(TwoSampleMR)
library(simulateGP)

infile <- commandArgs(T)[1]
outfile <- commandArgs(T)[2]
jid <- as.numeric(commandArgs(T)[3])

load(infile)

set.seed(jid)

sims <- length(l)
for(i in 1:sims)
{
	message(i)
	id <- runif(1) %>% as.character %>% gsub("0.","",.) %>% paste0(jid, ":", .)
	message(id)
	l[[i]]$estimates <- try(test_system(l[[i]], id))
	l[[i]]$id <- id
	names(l)[i] <- id
}

save(l, file=outfile)
