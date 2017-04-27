library(TwoSampleMR)
library(dplyr)

try_extract_outcome <- function(SNP, id, tries=3)
{
	i <- 1
	out <- 1
	class(out) <- 'try-error'
	while(class(out) == 'try-error' & i <= tries)
	{
		message("trying iteration ", i)
		out <- try(extract_outcome_data(SNP, id))
		i <- i + 1
	}
	return(out)
}

try_extract_instrument <- function(id, pval, tries=3)
{
	i <- 1
	out <- 1
	class(out) <- 'try-error'
	while(class(out) == 'try-error' & i <= tries)
	{
		message("trying iteration ", i)
		out <- try(extract_instruments(id, pval))
		i <- i + 1
	}
	return(out)
}
toggle_dev("test")

ao <- available_outcomes()
traits <- subset(ao, category=="Metabolites" & author=="Shin")

selected_ids <- traits$id


l <- list()
m <- list()
for(i in 1:length(selected_ids))
{
	message(selected_ids[i])
	l[[i]] <- try_extract_instrument(selected_ids[i], pval=5e-8)
	m[[i]] <- try_extract_outcome(l[[i]]$SNP, selected_ids)
}


save(l, m, traits, file="../data/shin_pairwise.rdata")
