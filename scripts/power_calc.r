source("functions.r")

#####

param <- expand.grid(
	sim = 1:100,
	ncase = c(1000, 10000, 50000),
	prev = c(0.01, 0.1, 0.5),
	h2phen = c(0.01, 0.1),
	h2pred = c(0.01, 0.1),
	vpred = c(0.01, 0.05)
)
param$nid <- param$ncase / param$prev



arguments <- commandArgs(T)
jid <- as.numeric(arguments[1])
chunks <- as.numeric(arguments[2])
out <- arguments[3]

chunksize <- ceiling(nrow(param) / chunks)
t1 <- (jid - 1) * chunksize + 1
t2 <- min(jid * chunksize, nrow(param))

message("total size: ", nrow(param))
message("running: ", t1, " to ", t2)

param <- param[t1:t2, ]


res <- list()
for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	dat <- make_system(param$nid[i], param$prev[i], param$h2phen[i], param$h2pred[i], param$vpred[i])
	right <- do_test(dat)
	left <- param[rep(i, nrow(right)), ]
	res[[i]] <- cbind(left, right)
}

res <- bind_rows(res)

save(res, file=out)

