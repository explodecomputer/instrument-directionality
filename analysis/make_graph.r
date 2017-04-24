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


load("../data/extract_everything2.rdata")

# Identify NULL elements in list

which(sapply(m, function(x) class(x)=="try-error"))
which(sapply(l, function(x) class(x)=="try-error"))

sapply(l, nrow)
table(sapply(l, function(x) sum(x$pval.exposure < 5e-8)) == 0)

d <- list()
for(i in 1:length(l))
{
	message(i)
	a <- subset(l[[i]], pval.exposure < 5e-8)
	d[[i]] <- harmonise_data(a, m[[i]])
}

m_ivw <- lapply(d, function(x) {
	mr(x, method_list=c("mr_wald_ratio", "mr_ivw"))
})

m_weighted_median <- lapply(d, function(x) {
	mr(x, method_list=c("mr_wald_ratio", "mr_weighted_median"))
})

m_mode <- lapply(d, function(x) {
	mr(x, method_list=c("mr_wald_ratio", "mr_mode"))
})

mr(d[[4]])

index <- which(sapply(m_ivw, function(x) nrow(x)==0))

length(m_ivw)
dim(traits)

temp <- m_ivw[-index]
traits2 <- traits[-index, ]
for(i in 1:nrow(traits2))
{
	temp[[i]] <- subset(temp[[i]], !id.outcome == traits2$id[i])
}

m_ivwd <- bind_rows(temp)
dim(m_ivwd)
head(m_ivwd)
(m_ivwd$pval < 0.05/nrow(m_ivwd), na.rm=TRUE) 

a <- subset(m_ivwd, pval < (0.05/nrow(m_ivwd)))
m_ivwd$fdr <- p.adjust(m_ivwd$pval, "fdr")

head(m_ivwd)
a <- subset(m_ivwd, fdr < 0.05)

dim(a)
table(a$exposure)

subset(a, exposure == "Intracranial volume || ENIGMA || 2015 || mm3")

