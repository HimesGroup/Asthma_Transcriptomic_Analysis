#!/usr/bin/R

library(rmeta)

###
# Function Module
###

# effect size-based method
ranef <- function(d,se,method="random") {
  x <- meta.summaries(d=d, se=se, method=method, conf.level=0.95, logscale=FALSE)
  summ <- x$summary
  se.summ <- x$se.summary
  nstudy <- length(d)
  pval <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
  res=c(nstudy,summ,se.summ,pval)
  names(res)=c("N","es_method.eff.size","es_method.std.err","es_method.pval")
  return(res)
}

# fisher combined p-values: compute x2 and convert to p-values
fisherp <- function(p) {
  keep <- (p > 0) & (p <= 1) & (!is.na(p))
  lnp <- log(p[keep])
  chisq <- (-2) * sum(lnp)
  df <- 2 * length(lnp)
  fisherp = pchisq(chisq, df, lower.tail = FALSE)
  if (length(p) <= 1) {fisherp=NA}
  names(fisherp)="p_method.pval"
  return(fisherp)
}

# Rankprod method
rankprodt <- function(rank) {
    keep <- (rank > 0) & (!is.na(rank))
    rankprodt <- (prod(rank[keep],na.rm=T))^(1/(length(rank[keep])))
    if (length(rank[keep]) <= 1) {rankprodt <- NA}
    names(rankprodt)="Score"
    return(rankprodt)
}

# Permutation-based Rankprod method
permvec_func <- function(vec){vec[!is.na(vec)] <- sample(vec[!is.na(vec)]);return(vec)} # permute a vector but keep the NA position
rankprodt_perm <- function(dat.list,nperm,rand=123) {
  dat <- data.frame(lapply(dat.list,function(x){x$rank})) # combind within-study ranks in a data frame
  # compute experimental rank product
  Gene=dat.list[[1]]$Gene
  org.rankprod=apply(dat,1,rankprodt) # original real rank product
  # permute rank product
  set.seed(rand)
  perm.res=rep(0,nrow(dat)) # store number of events that experimental rank product < permutated rank product
  for (i in 1:nperm) {
    tmp.dat <- dat
    tmp.dat=apply(tmp.dat,2,permvec_func)
    tmp.rankprod=apply(tmp.dat,1,rankprodt) # a temporary rankprod after permutation
    #print(tmp.rankprod)
    tmp.res=as.numeric(org.rankprod>=tmp.rankprod) # compare if the permutated rankprod is larger than the original one
    #print(tmp.res)
    perm.res=sapply(1:nrow(dat),function(x){sum(perm.res[x]+tmp.res[x])})
    i=i+1
  }
  res=data.frame(Gene=Gene,RP=org.rankprod,count=perm.res,rank_method.pval=perm.res/nperm)
  res$rank_method.qval=res$rank_method.pval*length(Gene)/rank(res$rank_method.pval)
  return(res[,c("Gene","rank_method.pval","rank_method.qval")])
}

# convert one-sided p-value to two-sided
pone2two <- function(dt) { # input table
    ref_study <- dt[order(-dt$Total,dt$P.Value),"Unique_ID"][1] # define the study with largest sample size as the reference study
    ref_sign <- sign(dt$logFC[which(dt$Unique_ID==ref_study)]) # define reference direction
    dt$stats <- rep(NA,length(dt$Total)) # a new column for converted p-value
    dt$stats[which(dt$Unique_ID==ref_study)] <- dt$P.Value[which(dt$Unique_ID==ref_study)] # use the two-sided p-value for reference study
    dt$stats[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)==ref_sign))] <- dt$P.Value[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)==ref_sign))]/2 # if the direction is the same as the reference direction, convert two-sided p-value to one-sided p-value
    dt$stats[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)!=ref_sign))] <- 1-dt$P.Value[which((dt$Unique_ID!=ref_study)&(sign(dt$logFC)!=ref_sign))]/2 # if the direction is opposite to the reference direction, convert two-sided p-value to one-sided p-value, and then convert to the other tail
    return(dt)
}

# select genes that shared within a certain number of studies
sharegene_func <- function(dat.list, minnum) {
  if (missing(minnum)) {minnum=length(dat.list)}
  allgenes=table(unname(unlist(lapply(dat.list,function(x){as.character(x$Gene)})))) # compute gene frequency
  allgenes <- allgenes[!(is.na(allgenes) | allgenes=="NA")]
  genes=names(allgenes)[which(allgenes>=minnum)] # obtain genes appear in how many number of datasets
  return(genes)
}

# generate output table with converted p-values by gene
pconvtb_func <- function(dat.list,dat.info) { # dat.info: data information sheet with selected studies; minnum: genes shared within at leat the particular number of studies
  genes=sharegene_func(dat.list=dat.list,minnum=1) # genes from dat.list
  # obtain tables with converted p-values
  res_bygene=list()
  for (curr_gene in genes) {
    curr_gene <- as.character(curr_gene)
    output.table <- data.frame() # initiate output.table to combine results for curr_gene
    output.table <- do.call(rbind,lapply(dat_list,function(x){x[which(x$Gene==curr_gene),]}))
    output.table <- data.frame(Unique_ID=row.names(output.table),output.table)
    row.names(output.table) <- NULL
    output.table$Total=dat.info$Total[dat.info$Unique_ID%in%output.table$Unique_ID]
    if (nrow(output.table)>2) {
      output.table <- pone2two(output.table)
      res_bygene[[curr_gene]] <- output.table
    }
  }
  return(res_bygene)
}

# Integration via effect size-, p-value-, and rank-based methods
es_func <- function(stat.list) {
  es.res=lapply(stat.list,function(x){d=x$logFC;se=x$SD;ranef(d=d,se=se)})
  res=do.call(rbind,es.res)
  res=data.frame(Gene=names(stat.list),res)
  res$es_method.qval=p.adjust(res$es_method.pval,method="BH")
  row.names(res) <- NULL
  return(res)
}

p_func <- function(stat.list) {
  p.res=lapply(stat.list,function(x){fisherp(x$stats)})
  res=do.call(rbind,p.res)
  res=data.frame(Gene=names(stat.list),res)
  res$p_method.qval=p.adjust(res$p_method.pval,method="BH")
  row.names(res) <- NULL
  return(res)
}

summstat_func <- function(stat.list,dat.list) {
  es.res=es_func(stat.list)
  p.res=p_func(stat.list)
  rank.res=rankprodt_perm(dat.list,nperm=10000,rand=123)
  Reduce(function(x,y){merge(x,y,by="Gene",all=T)},list(es.res,p.res,rank.res))
}

###
# Integration Analysis
###

# read in data information
datinfo <- readRDS("data/datainfo.RDS")

# select glucocorticoid blood cell studies
datinfo.sub <- datinfo[datinfo$Asthma%in%"GC"&datinfo$Tissue%in%c("MCF10A-Myc","chALL","MACRO","PBMC","LCL"),]

# read in gene results
dat <- readRDS("data/example.RDS")
dat_list=list()
for (i in datinfo.sub$Unique_ID) {dat_list[[i]]=dat[[i]]}

# obtain shared genes within at least 2 studies
genes=sharegene_func(dat.list=dat_list,minnum=2)

# read in gene results
dat_list=list()
for (i in datinfo.sub$Unique_ID) {dat_list[[i]]=dat[[i]]}

# subset data with genes that appear in at least 2 studies
dat_list=lapply(dat_list,function(x){x[x$Gene%in%genes,]})

# compute integration results
summstat_func(stat.list=stat_list,dat.list=dat_list)
