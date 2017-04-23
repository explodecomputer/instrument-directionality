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

load("../data/extract_everything.rdata")

ao <- available_outcomes()
selected_ids <- scan("../data/selected_ids.txt", what=numeric())
traits <- subset(ao, id %in% selected_ids)

# Select most reliable studies for diseases and risk factors
# Exclude lipids
traits_orig <- filter(ao, 
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

selected_ids <- traits$id
original_ids <- traits_orig$id

remaining_ids <- selected_ids[!selected_ids %in% original_ids]


l1 <- list()
m1 <- list()
for(i in 1:length(selected_ids))
{
	index <- which(original_ids %in% selected_ids[i])
	if(length(index) == 1)
	{
		message(selected_ids[i], " already present")
		l1[[i]] <- l[[index]]
		m1[[i]] <- m[[index]]
		temp <- try_extract_outcome(l1[[i]]$SNP, remaining_ids)
		if(class(temp) != "try-error")
		{
			m1[[i]] <- subset(m1[[i]], id.outcome %in% selected_ids)
			m1[[i]] <- rbind(m1[[i]], temp)
		}
	} else {
		message(selected_ids[i], " not present")
		l1[[i]] <- try_extract_instrument(selected_ids[i], pval=5e-6)
		m1[[i]] <- try_extract_outcome(l1[[i]]$SNP, selected_ids)
	}
}

l <- l1
m <- m1
save(l, m, traits, file="../data/extract_everything2.rdata")
