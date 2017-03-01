library(tidyverse)

load("../results/power_calc.rdata")

pow_var <- filter(res, model == "RF") %>%
	group_by(ncase, prev, h2phen, vpred) %>%
	dplyr::summarise(pow=sum(pval < 0.05) / n(), n=n())

pow_prs <- filter(res, model == "PRS of RF") %>%
	group_by(ncase, prev, h2phen, h2pred, vpred) %>%
	dplyr::summarise(pow=sum(pval < 0.05) / n(), n=n())


pow_var$lab1 <- paste0("prevalence = ", pow_var$prev)
pow_prs$lab1 <- paste0("prevalence = ", pow_prs$prev)
pow_var$lab2 <- paste0("# cases = ", pow_var$ncase)
pow_prs$lab2 <- paste0("# cases = ", pow_prs$ncase)
pow_var$lab3 <- paste0("r2_pred = ", pow_var$h2pred)
pow_prs$lab3 <- paste0("r2_pred = ", pow_prs$h2pred)


ggplot(pow_var, aes(x=as.factor(vpred), y=pow)) +
geom_bar(stat="identity", aes(fill=as.factor(h2phen)), position="dodge") +
facet_grid(lab2 ~ lab1) +
labs(x="Variance explained by risk factor", y="Power (100 simulations, alpha=0.05)", fill="Variance explained\nby genetic score")
ggsave("../images/power_rf.pdf")

ggplot(subset(pow_prs, ncase==50000), aes(x=as.factor(vpred), y=pow)) +
geom_bar(stat="identity", aes(fill=as.factor(h2phen)), position="dodge") +
facet_grid(lab3 ~ lab1) +
labs(x="Variance explained by risk factor", y="Power (100 simulations, alpha=0.05)", fill="Variance explained\nby genetic score")
ggsave("../images/power_prs_of_rf.pdf")
