library(tidyverse)

filter_relevant <- function(feature_table) {
  feature_table %>% filter(mzmed >= 70,mzmed <= 700, rtmed >= 0.4, rtmed <= 16)
}

add_blank_matches <- function(samplesFilePath,pblanksFilePath) {
  
  samples <- filter_relevant(read.table(samplesFilePath,sep = "\t",header=TRUE))
  pblanks <- filter_relevant(read.table(pblanksFilePath,sep = "\t",header=TRUE))
    
  blank_pvalue <- rep(0.0, nrow(samples))
  blank_pvalue2 <- rep(0.0, nrow(samples))
  blank_matches <- rep(0.0, nrow(samples))
  blank_max <- rep(0.0, nrow(samples))
  blank_max2 <- rep(0.0, nrow(samples))
  sample_min <- rep(0.0, nrow(samples))
  sample_min0 <- rep(0.0, nrow(samples))
  sample_median <- rep(0.0, nrow(samples))
  difmin <- rep(0.0, nrow(samples))
  difmedian <- rep(0.0, nrow(samples))
  difmin2 <- rep(0.0, nrow(samples))
  difmedian2 <- rep(0.0, nrow(samples))
  bf <- rep(0.0, nrow(samples))
  bf2 <- rep(0.0, nrow(samples))
  
  for (i in 1:nrow(samples)) {
    mzminSample <- samples$mzmin[i]
    mzmaxSample <- samples$mzmax[i]
    rtminSample <- samples$rtmin[i]
    rtmaxSample <- samples$rtmax[i]
      
    for (j in 1:nrow(pblanks)) {
      mzmedBlank <- pblanks$mzmed[j]
      rtmedBlank <- pblanks$rtmed[j]
      
      if (mzmedBlank < mzminSample || mzmedBlank > mzmaxSample) 
        next
      
      if (rtmedBlank < rtminSample || rtmedBlank > rtmaxSample) 
        next
      
      if (mean(unlist(samples[i,] %>% select(starts_with("QT_"))))>0 & max(pblanks[j,] %>% select(starts_with("QT_")))==0) 
        next
      
      #intensities <- pblanks[j,] %>% select(starts_with("QT_"))
      #average <- mean(unlist(intensities))
      #blank_intensity[i] <- average
      blank_matches[i] <- blank_matches[i] + 1
      bill=unlist(samples[i,] %>% select(starts_with("QT_")))
      sample_min[i] <- min(bill[bill>0])
      sample_min0[i] <- min(bill)
      sample_median[i] <- median(bill)
      if (blank_max[i]==0) {
        blank_max[i] <- max(pblanks[j,] %>% select(starts_with("QT_")))
        difmin[i] <- sample_min[i]-blank_max[i]
        difmedian[i] <- sample_median[i]-blank_max[i]
        bf[i] <- sample_min[i]-(3*blank_max[i])
      } else {
        blank_max2[i] <- max(pblanks[j,] %>% select(starts_with("QT_")))
        difmin2[i] <- sample_min[i]-blank_max2[i]
        difmedian2[i] <- sample_median[i]-blank_max2[i]
        bf2[i] <- sample_min[i]-(3*blank_max2[i])
      } 
      mitch=wilcox.test(x=unlist(samples[i,] %>% select(starts_with("QT_"))), mu=max(pblanks[j,] %>% select(starts_with("QT_"))), alternative = "greater", paired = F, var.equal = F, conf.level = 0.95)
      if (blank_matches[i]<2) {
        blank_pvalue[i] <- mitch$p.value
      } else {
        blank_pvalue2[i] <- mitch$p.value
      }
    }
  }
  
  matched <<- cbind(samples, sample_min0, sample_median, sample_min, blank_pvalue, blank_max, difmin, difmedian, blank_pvalue2, blank_max2, difmin2, difmedian2, blank_matches, bf, bf2)
  write.table(matched,file=str_replace(samplesFilePath,".tsv","_pbmatched.tsv"),sep="\t",row.names = FALSE)
}