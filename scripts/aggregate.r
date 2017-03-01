library(tidyverse)

l <- list()
for(i in 1:16)
{
	load(paste0("../results/power_calc_", i, ".rdata"))
	l[[i]] <- res
}
res <- bind_rows(res)
save(res, file="../results/power_calc.rdata")


