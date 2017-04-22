library(tidyverse)
load("../results/simulate2.rdata")

# Validity

validity_s <- group_by(validity, hypothesis, strategy) %>%
	summarise(
		n=sum(value), 
		direct=sum(value[type=="valid"])/n,
		precursor=sum(value[type=="valid_precursor"])/n,
		confounder=sum(value[type=="confounder"])/n,
		reverse=sum(value[type=="reverse"])/n
	)



# FDR

res <- inner_join(res, param, by="sim")

resnull <- bind_rows(
	filter(res, hypothesis=="xy", eff_x.y == 0),
	filter(res, hypothesis=="yx", eff_x.y != 0)
)

resnull_s <- group_by(resnull, Method, strategy, hypothesis) %>%
	summarise(fdr = sum(P < 0.05)/n())
levels(resnull_s$hypothesis) <- c("No causal effect", "Reverse cause")

ggplot(resnull_s, aes(x=Method, y=fdr)) +
geom_bar(stat="identity", position="dodge", aes(fill=strategy)) +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
scale_fill_brewer(type="qual") +
geom_hline(yintercept=0.05, linetype="dotted") +
facet_grid(hypothesis ~ .)


# Bias

resxy <- filter(res, hypothesis == "xy", eff_x.y != 0)
resxy$eff_bin <- cut(resxy$eff_x.y, 10)

resxy_s <- group_by(resxy, Method, strategy, eff_bin) %>%
	summarise(tdr=sum(P < 0.05)/n(), bias=mean(Estimate - eff_x.y), n=n(), bias_se=sd(Estimate - eff_x.y)/sqrt(n))
resxy_s

ggplot(resxy_s, aes(x=eff_bin, y=bias, group=strategy)) +
geom_point(aes(colour=strategy)) +
geom_errorbar(aes(ymin=bias-bias_se*1.96, ymax=bias+bias_se*1.96, colour=strategy), width=0) +
geom_line(aes(colour=strategy)) +
facet_grid(. ~ Method) +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
scale_colour_brewer(type="qual")


# Power

ggplot(resxy_s, aes(x=eff_bin, y=tdr, group=strategy)) +
geom_point(aes(colour=strategy)) +
geom_line(aes(colour=strategy)) +
facet_grid(. ~ Method) +
scale_colour_brewer(type="qual") +
theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5))


# Rucker

table(ruck$model)
table(ruck$model, ruck$hypothesis, ruck$strategy)
