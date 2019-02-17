library(dplyr)
library(multidplyr)
library(pROC)

load("../results/sim3/agg-estimates.rdata")
load("../results/sim3/agg-param.rdata")
load("../results/sim3/agg-heterogeneity.rdata")
load("../results/sim3/agg-plei_summary.rdata")

plei <- group_by(plei_summary, id) %>%
	summarise(nonpl = sum(nonpl * effd / sum(effd)))

res <- inner_join(estimates, param, by="id")
res <- subset(res, nsnp >= 5)
res$code <- paste0(res$id, res$hypothesis)
res <- inner_join(res, plei, by=c("code"="id"))
res <- subset(res, !is.na(res$nonpl))
res$strategy <- paste(res$selection, as.numeric(res$steiger_filtered), as.numeric(res$outlier_filtered))
res$n <- res$nidx
res$n[res$hypothesis == "y"] <- res$nidy[res$hypothesis == "y"]
res$nbin <- cut(res$n, breaks=10, labels=FALSE)
res$nbin <- res$nbin * 50000
res$pbin <- cut(res$nonpl, breaks=quantile(res$nonpl, probs = seq(0,1,0.2), include.lowest = T), labels=FALSE)
res$nsnpbin <- cut(res$nsnp, breaks=quantile(res$nsnp, probs = seq(0,1,0.2), include.lowest = T), labels=FALSE)
res <- filter(res, !is.na(nsnpbin), !is.na(pbin), !is.na(nbin))
res$causal <- as.numeric(res$hypothesis == "x" & res$eff_x.y != 0)
save(res, file="res.rdata")

##

cluster <- create_cluster(11)
set_default_cluster(cluster)

rocs <- res %>% filter(selection == "e") %>% 
partition(method, strategy) %>%
cluster_library("pROC") %>%
cluster_library("dplyr") %>%
do({
	r <- pROC::roc(.$causal, abs(.$b/.$se), ci=TRUE)
	tibble(
		n=nrow(.), 
		auc=r$auc[1],
		ci_lo=r$ci[1],
		ci_up=r$ci[3]
	)
})

rocs <- as_data_frame(rocs)
rocs$strategy[rocs$strategy == "e 0 0"] <- "Tophits"
rocs$strategy[rocs$strategy == "e 1 0"] <- "DF"
rocs$strategy[rocs$strategy == "e 0 1"] <- "HF"
rocs$strategy[rocs$strategy == "e 1 1"] <- "DF + HF"
rocs$method <- factor(rocs$method)
rocs$method <- factor(rocs$method, levels=levels(rocs$method)[c(7,2,6,1,5,8,10,3,9,11,4)])

rocs2 <- res %>% filter(selection == "e") %>%
partition(method, strategy, nsnpbin, pbin) %>%
cluster_library("pROC") %>%
do({
	r <- pROC::roc(.$causal, abs(.$b/.$se), ci=TRUE)
	tibble(
		n=nrow(.), 
		auc=r$auc[1],
		ci_lo=r$ci[1],
		ci_up=r$ci[3]
	)
})
rocs2 <- as_data_frame(rocs2)
rocs2$strategy[rocs2$strategy == "e 0 0"] <- "Tophits"
rocs2$strategy[rocs2$strategy == "e 1 0"] <- "DF"
rocs2$strategy[rocs2$strategy == "e 0 1"] <- "HF"
rocs2$strategy[rocs2$strategy == "e 1 1"] <- "DF + HF"
rocs2$method <- factor(rocs2$method)
rocs2$method <- factor(rocs2$method, levels=levels(rocs2$method)[c(7,2,6,1,5,8,10,3,9,11,4)])

save(rocs, rocs2, file="rocs.rdata")

