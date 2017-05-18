### Jinliang Yang


library("data.table")
mlm <- fread("largedata/mlm_filter90.txt", data.table=FALSE)

ts <- unique(mlm$Trait)

### false discovery rate cutoff
CF = 0.05 # 0.1

out <- list()
idx_cf <- c() ## cutoff idx
for(i in 1:length(ts)){
  sub <- subset(mlm, Trait %in% ts[i])
  sub <- subset(sub, Chr %in% 1:10)
  names(sub)[3:4] <- c("chr", "pos")
  
  sub$qval <- p.adjust(sub$p, method ="fdr")
  
  message(sprintf("###>>> trait [ %s ]: [ %s/%s ] qval < [ %s ]", 
                  ts[i], sum(sub$qval < CF, na.rm=T), nrow(sub), CF))
  out[[ts[i]]] <- sub
  
  if(sum(sub$qval < CF, na.rm=T) > 0 ){
    

    idx <- which.min( abs(sub$qval - CF))
    idx_cf <- c(idx_cf, idx)
  }else{
    idx_cf <- c(idx_cf, 0)
  }
  
}



########
source("lib/quickMHTplot.R")
source("lib/newpos.R")

# location: 129.186.85.7

pdf("graphs/mht_plots_filter90_fdr0.05.pdf", width=10, height=4)

for(j in 1:length(out)){
  t <- names(out)[j]
  
  sub <- out[[t]]
  sub$log10p <- -log10(sub$p)
  if(idx_cf[j] > 0){
    
    sub1 <- subset(sub, qval > CF)
    sub2 <- subset(sub, qval < CF)
    
    ## determine the value of the cutoff
    mycf <- (min(sub2$log10p) + max(sub1$log10p))/2
  }
  
  
  sub <- as.data.frame(sub)
  sub <- subset(sub, log10p > 1)
  quickMHTplot(res=sub, cex=.6, pch=16, col=rep(c("slateblue", "cyan4"), 5), 
               GAP=5e+06, yaxis=NULL, main=t, ylab="-log10(p-value)",
               col2plot="log10p")
  if(idx_cf[j] > 0){
    abline(h=mycf, col="red", lty=2)  
  }
}

dev.off()


