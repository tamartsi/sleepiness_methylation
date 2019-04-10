
## betas and stdErrors are data frames or matrices, each with p rows n.strat columns. 
## Each row corresponds to a specific SNP/methylation sites. Each column corresponds to 
## a specific strata (males/females, African Americans, European Americans, etc.).
## The order has to match between betas and stdErrors!
## the rownames of betas (and stdErrors) should be probe IDs, or SNPs, or somthing similar. 

metaAnalysis <- function(betas, stdErrors, ns){
	
	betas <- as.matrix(betas)
	stdErrors <- as.matrix(stdErrors)
	if (!all(rownames(betas) == rownames(stdErrors))) stop("rownames of betas and stdErrors do not match")

  # get weights
  weights <- 1 / stdErrors^2
  sum.wt <- rowSums(weights, na.rm=T) # SNPs that are all NA become 0 here

  # number of contributing groups: non-NA betas minus 1
  n.strat = rowSums(!is.na(betas))
  df = as.integer(ifelse(n.strat-1 < 1, 1, n.strat-1))


  N.total <- rowSums(ns*(!is.na(betas)), na.rm = T)


  # calculate meta-analysis stats
  beta.meta <- rowSums(betas * weights, na.rm=T) / sum.wt
  se.meta <- sqrt(1/sum.wt)

  Z <- beta.meta / se.meta
  # is this correct?
  pval.Z <- 2*pnorm(abs(Z), lower.tail=FALSE) # lower.tail gives better precision than 1-pnorm(x)

  # direction
  direction <- apply(ifelse(is.na(betas), "?",
                            ifelse(betas == 0, "0",
                                   ifelse(betas < 0, "-", "+"))), 1, paste, collapse="")


  # heterogeneity test
  Q <- rowSums(weights * (betas - beta.meta)^2, na.rm=T)
  ISq.Q <- 100 * pmax((Q-df)/Q, 0) # i squared statistic
  pval.Q <- pchisq(Q, df=df, lower.tail=F) # lower.tail gives better precision than 1-pchisq(x)

  # meta data frame
  meta <- data.frame(probeID=rownames(betas), beta.meta=beta.meta, se.meta=se.meta, Z=Z, pval.Z=pval.Z, direction=direction, Q=Q, pval.Q=pval.Q, ISq.Q=ISq.Q, n.strat=n.strat, N.total = N.total, stringsAsFactors=F)
  
   # updates for NAs
  meta[meta$n.strat == 0, c("beta.meta", "se.meta", "Z", "pval.Z", "Q", "pval.Q", "ISq.Q")] <- NA
  meta
}



