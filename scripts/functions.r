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

make_system <- function(nid, prev, h2disease, h2pred, vpred)
{
	prs_disease <- rnorm(nid)
	prs_pred <- rnorm(nid)
	pred <- makePhen(sqrt(h2pred), prs_pred)
	e <- rnorm(nid)
	ve <- 1 - h2disease - vpred

	px <- scale(prs_disease * sqrt(h2disease) + pred * sqrt(vpred) + e * sqrt(ve))
	prob_disease <- gx_to_gp(px, 1, 1-prev)
	disease <- rbinom(nid, 1, prob_disease)
	dat <- data.frame(
		disease=disease,
		prob_disease=prob_disease,
		prs_disease=prs_disease,
		pred=pred,
		prs_pred=prs_pred
	)
}


do_test <- function(dat)
{
	dat <- filter(dat, disease==1)
	a <- fastAssoc(dat$prs_disease, dat$pred)
	b <- fastAssoc(dat$prs_disease, dat$prs_pred)
	ab <- rbind(as.data.frame(a), as.data.frame(b))
	ab$model <- c("RF", "PRS of RF")
	return(ab)
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
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, n=n
	))
}
