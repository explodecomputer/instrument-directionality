library(tidyverse)
library(TwoSampleMR)

ao <- available_outcomes()


# Select most reliable studies for diseases and risk factors
# Exclude lipids
traits <- filter(ao, 
		category %in% c("Disease", "Risk factor"),
		(grepl("Mixed", population) | grepl("European", population)),
		mr == 1,
		sample_size > 1000,
		nsnp > 95000,
		access %in% c("immunobase_users", "public")
	) %>%
	arrange(desc(year), desc(nsnp)) %>%
	distinct(trait, .keep_all=TRUE) %>%
	filter(! id %in% c(299:302))


load("../data/extract_everything.rdata")

# Identify NULL elements in list

i <- which(sapply(l, is.null))
class(m[[i]])
j <- which(sapply(l, function(x) class(x)=="try-error"))

which(sapply(m, function(x) class(x)=="try-error"))
m[[36]] <- try(extract_outcome_data(l[[36]]$SNP, traits$id))