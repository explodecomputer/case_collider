# library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

make_geno <- function(nid, nsnp, maf)
{
	return(matrix(rbinom(nid * nsnp, 2, maf), nid, nsnp))
}

range01 <- function(x)
{
	(x-min(x))/(max(x)-min(x))
}

gx_to_gp <- function(gx, h2x, prev)
{
	x_prime <- qnorm(prev, 0, 1)
	p <- pnorm(x_prime, mean=gx, sd = sqrt(1 - h2x), lower.tail=FALSE)
	return(p)
}


# G should scaled
# var(eff) is the h2x
simulate <- function(G, eff, prevalence, prop_discovered)
{
	nid <- nrow(G)
	nsnp <- ncol(G)
	h2x <- var(eff)
	gx_true <- as.numeric(scale(G %*% eff)) * sqrt(h2x)
	prob_disease <- gx_to_gp(gx_true, h2x, 1-prevalence)
	disease <- rbinom(nid, 1, prob_disease)
	eff_pred <- eff
	eff_pred[sample(1:nsnp, nsnp * (1-prop_discovered))] <- 0
	gx_pred <- as.numeric(G %*% eff_pred / sqrt(nsnp))
	dat <- data.frame(gx_true=gx_true, gx_pred=gx_pred, prob_disease=prob_disease, disease=disease)
	return(dat)
}

make_system <- function(nid, prev, h2phen, h2pred, vpred)
{
	prs_phen <- rnorm(nid)
	prs_pred <- rnorm(nid)
	pred <- makePhen(sqrt(h2pred), prs_pred)
	px <- prs_phen * sqrt(h2phen) + pred * sqrt(h2pred)
	prob_disease <- gx_to_gp(px, var(px), 1-prev)
	disease <- rbinom(nid, 1, prob_disease)
	dat <- data.frame(
		disease=disease,
		prs_disease=prs_phen,
		pred=pred,
		prs_pred=prs_pred
	)
}


do_test <- function(dat)
{
	dat <- filter(dat, disease==1)
	a <- fastAssoc(dat$prs_disease, dat$pred)$pval
	b <- fastAssoc(dat$prs_disease, dat$prs_pred)$pval
	return(c(a, b))
}



cross_plot <- function(x, o, xlab="Values (low to high)", ylab="", title="", xname="GRS", oname="Disease")
{

	d <- data.frame(
		value = c(range01(o), range01(x)),
		key = c(rep(oname, length(o)), rep(xname, length(x))),
		gr = rep(1:length(x), times=2)
	)
	d$key <- factor(d$key, levels=c(xname, oname))
	ggplot(d, aes(x=value, y=key)) +
	geom_line(aes(group=gr), alpha=0.1) +
	geom_point(aes(colour=key)) +
	labs(x=xlab,y=ylab,title=title) +
	scale_colour_discrete(guide=FALSE) +
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
}


makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors^2) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

chooseEffects <- function(nsnp, totvar, sqrt=TRUE)
{
	eff <- rnorm(nsnp)
	eff <- sign(eff) * eff^2
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out)
}

fastAssoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 5)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}



# do_test(a)

# # simulate trait
# # cov2 = u + e
# # px = prs + cov1 + u + e 

# nid <- 100000
# prev <- 0.2
# h2 <- 0.1
# v1 <- 0.05
# v2 <- 0.05
# uv2 <- 0.05
# up <- 0.05

# prs <- rnorm(nid)
# cov1 <- rnorm(nid)
# u <- rnorm(nid)
# cov2 <- makePhen(sqrt(uv2), u)
# px <- prs * sqrt(h2) + cov1 * sqrt(v1) + u * sqrt(up)
# prob_disease <- gx_to_gp(px, var(px), 1-prev)
# disease <- rbinom(nid, 1, prob_disease)

# summary(lm(cov1[disease==1] ~ cov2[disease==1]))
# summary(lm(cov1[disease==1] ~ u[disease==1]))
# summary(lm(cov1[disease==1] ~ prs[disease==1]))




# a <- make_system(10000, 0.1, 0.1, 0.1, 0.1)
# table(a$disease)

param <- expand.grid(
	sim = 1:100,
	ncase = c(1000, 10000, 50000),
	prev = c(0.01, 0.1),
	h2phen = c(0.01, 0.1),
	h2pred = c(0.01, 0.1),
	vpred = c(0.01, 0.05),
	pval_var = NA,
	pval_prs = NA
)
param$nid <- param$ncase / param$prev

set.seed(100)
for(i in 1:nrow(param))
{
	message(i, " of ", nrow(param))
	dat <- make_system(param$nid[i], param$prev[i], param$h2phen[i], param$h2pred[i], param$vpred[i])
	res <- do_test(dat)
	param$pval_var[i] <- res[1]
	param$pval_prs[i] <- res[2]
}

save(param, file="power_calc.rdata")

