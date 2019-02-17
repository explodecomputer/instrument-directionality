make_resnull <- function(res)
{
	bind_rows(
		filter(res, hypothesis=="x", eff_x.y == 0),
		filter(res, hypothesis=="y")
	)	
}

make_resxy <- function(res)
{
	filter(res, hypothesis == "x", eff_x.y != 0)

}


simeval <- function(res)
{

	require(tidyverse)
	require(ggrepel)

	# Organise data
	resnull_s_d <- group_by(resnull, method, strategy, hypothesis) %>%
		summarise(fdr = sum(pval < 0.05)/n())
	levels(resnull_s_d$hypothesis) <- c("No causal effect", "Reverse cause")

	resnull_s <- group_by(resnull, method, strategy) %>%
		summarise(fdr = sum(pval < 0.05)/n())

	resxy$eff_bin <- cut(resxy$eff_x.y, 3)

	resxy_s <- group_by(resxy, method, strategy, eff_bin) %>%
		summarise(tdr=sum(pval < 0.05)/n(), bias=mean(b - eff_x.y), n=n(), bias_se=sd(b - eff_x.y)/sqrt(n))



	temp <- resxy %>%
		group_by(method, strategy) %>%
		summarise(power=sum(pval < 0.01)/n())
	# temp <- inner_join(temp, meth, by="Method")
	# temp$Method <- as.factor(temp$Method)
	# temp$Method <- factor(temp$Method, levels=as.character(meth$Method))

	p_power <- ggplot(temp %>% filter(!grepl("^o", strategy)), aes(y=power, x=method)) +
		geom_point(position=position_dodge(width=0.3), aes(colour=strategy), size=3) +
		theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
		scale_colour_brewer(type="qual") +
		labs(colour="Instrument\nselection")



	temp <- resnull_s
	# temp <- inner_join(temp, meth, by="Method")
	# temp$Method <- as.factor(temp$Method)
	# temp$Method <- factor(temp$Method, levels=as.character(meth$Method))

	p_fdr <- ggplot(temp %>% filter(!grepl("^o", strategy)), aes(y=fdr, x=method)) +
		geom_point(position=position_dodge(width=0.3), aes(colour=strategy), size=3) +
		theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
		scale_colour_brewer(type="qual") +
		labs(colour="Instrument\nselection")


	temp1 <- resxy %>%
		group_by(method, strategy) %>%
		summarise(power=sum(pval < 0.01)/n())
	# temp1 <- inner_join(temp1, meth, by="Method")
	# temp1$Method <- as.factor(temp1$Method)
	# temp1$Method <- factor(temp1$Method, levels=as.character(meth$Method))

	temp2 <- resnull_s
	# temp2 <- inner_join(temp2, meth, by="Method")
	# temp2$Method <- as.factor(temp2$Method)
	# temp2$Method <- factor(temp2$Method, levels=as.character(meth$Method))

	temp <- merge(temp1, temp2, by=c("method", "strategy"))
	temp$meth2 <- paste0(temp$method, " - ", temp$strategy)
	p_powerfdr <- ggplot(subset(temp, !grepl("^o", strategy)), aes(x=power, y=fdr)) +
		geom_point(aes(colour=strategy), size=3) +
		geom_text_repel(aes(label=method)) +
		scale_colour_brewer(type="qual") +
		labs(x="Power (higher is better)", y="FDR (lower is better)", colour="Instrument\nselection")

	temp <- resxy %>%
		group_by(method, strategy) %>%
		summarise(bias=mean(b - eff_x.y), se=sd(b-eff_x.y))
	# temp <- inner_join(temp, meth, by="Method")
	# temp$Method <- as.factor(temp$Method)
	# temp$Method <- factor(temp$Method, levels=as.character(meth$Method))

	p_bias <- ggplot(temp %>% filter(!method %in% c("Steiger null", "Wald ratio"), !grepl("^o", strategy)), aes(y=bias, x=method)) +
		geom_point(position=position_dodge(width=0.3), aes(colour=strategy), size=3) +
		geom_errorbar(position=position_dodge(width=0.3), aes(ymin=bias-se, ymax=bias+se, colour=strategy), width=0) +
		theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
		scale_colour_brewer(type="qual") +
		labs(colour="Instrument\nselection")

	return(list(p_power=p_power, p_fdr=p_fdr, p_powerfdr=p_powerfdr, p_bias=p_bias))

}



make_optim_dataset <- function(optvec, res, metrics, ncores)
{
	# For each method calculate the abs(bias) for each simulation
	# Fit RF of metrics against abs(bias)

	require(tidyverse)
	require(randomForest)
	require(parallel)

	res$optvec <- optvec

	temp <- subset(res, strategy != "oracle", select=c(optvec, Method, strategy, hypothesis, sim))
	temp$method <- paste0(temp$Method, " - ", temp$strategy)
	# temp$optvec <- abs(temp$Estimate - temp$eff_x.y)
	temp <- subset(temp, select=-c(Method, strategy))

	sp <- subset(temp, method == method[1], select=c(sim, hypothesis))
	sp$tt <- runif(1:nrow(sp))

	temp <- inner_join(temp, sp, by=c("hypothesis", "sim"))
	ind <- temp$tt > 0.33

	l <- split(temp[ind,], temp$method[ind])

	rf <- mclapply(l, function(x){
		message(x$method[1])
		x <- inner_join(x, metrics, by=c("hypothesis", "sim")) %>%
			dplyr::select(-c(hypothesis, sim, method, tt))

		return(randomForest(optvec ~ ., x, ntree=60, mtry=15, importance=TRUE, do.trace=TRUE))
	}, mc.cores=ncores)

	out <- list(rf=rf, sp=sp)
	return(out)
}




make_optim_dataset2 <- function(optvec, res, metrics, ncores)
{
	# For each method calculate the abs(bias) for each simulation
	# Fit RF of metrics against abs(bias)

	require(tidyverse)
	require(MASS)
	require(parallel)

	res$optvec <- optvec

	temp <- subset(res, strategy != "oracle", select=c(optvec, Method, strategy, hypothesis, sim))
	temp$method <- paste0(temp$Method, " - ", temp$strategy)
	# temp$optvec <- abs(temp$Estimate - temp$eff_x.y)
	temp <- subset(temp, select=-c(Method, strategy))

	sp <- subset(temp, method == method[1], select=c(sim, hypothesis))
	sp$tt <- runif(1:nrow(sp))

	temp <- inner_join(temp, sp, by=c("hypothesis", "sim"))
	ind <- temp$tt > 0.33

	l <- split(temp[ind,], temp$method[ind])

	rf <- mclapply(l, function(x){
		message(x$method[1])
		x <- inner_join(x, metrics, by=c("hypothesis", "sim")) %>%
			dplyr::select(-c(hypothesis, sim, method, tt))

		return(lda(optvec ~ ., x))
	}, mc.cores=ncores)

	out <- list(rf=rf, sp=sp)
	return(out)
}




test_predictor <- function(rf, sp, res, metrics)
{
	temp <- subset(res, strategy != "oracle", select=c(Estimate, eff_x.y, Method, strategy, hypothesis, sim))
	temp$method <- paste0(temp$Method, " - ", temp$strategy)
	temp$bias <- abs(temp$Estimate - temp$eff_x.y)
	temp <- subset(temp, select=-c(Method, strategy, Estimate, eff_x.y))

	met <- inner_join(metrics, sp, by=c("hypothesis", "sim")) %>%
		filter(tt <= 0.33)
	nom <- names(rf)

	d <- data.frame(method=nom, rsq=NA)

	for(i in 1:length(nom))
	{
		message(nom[i])
		pr <- predict(rf[[nom[i]]], subset(met, select=-c(hypothesis, sim)))
		if(is.list(pr)) 
		{
			pr <- as.numeric(pr$x)
		}
		x <- tibble(pr=pr, hypothesis=met$hypothesis, sim=met$sim)
		y <- subset(temp, method == nom[i])
		xy <- inner_join(x, y, by=c("hypothesis", "sim"))
		d$rsq[d$method == nom[i]] <- cor(xy$pr, xy$bias)^2
	}
	return(d)
}

add_method <- function(rf, sp, res, metrics)
{
	# PREDICT BIAS FOR EVERY METHOD
	temp <- subset(res, strategy != "oracle", select=c(Estimate, eff_x.y, Method, strategy, hypothesis, sim))
	temp$method <- paste0(temp$Method, " - ", temp$strategy)
	temp$bias <- abs(temp$Estimate - temp$eff_x.y)
	temp <- subset(temp, select=-c(Method, strategy, Estimate, eff_x.y))

	met <- inner_join(metrics, sp, by=c("hypothesis", "sim")) %>%
		filter(tt <= 0.33)
	nom <- names(rf)

	d <- tibble(hypothesis=met$hypothesis, sim=met$sim)

	for(i in 1:length(nom))
	{
		message(nom[i])
		p <- predict(rf[[nom[i]]], subset(met, select=-c(hypothesis, sim)))
		if(is.list(p))
		{
			p <- as.numeric(p$x)
		}
		d[[nom[i]]] <- p
	}

	sel <- tibble(
		hypothesis=met$hypothesis, sim=met$sim,
		selmethod = nom[apply(d[,-c(1,2)], 1, function(x) which.max(x))]
	)

	res2 <- inner_join(res, sel, by=c("hypothesis", "sim")) 

	res3 <- res2 %>%
		filter(selmethod == method)
	res3$method <- "MoE - RF"
	res3$strategy <- "Mixed"
	res3 <- bind_rows(res2, res3)
	return(res3)
}
