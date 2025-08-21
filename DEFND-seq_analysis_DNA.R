#This is the file for processing #39 the DNA sample
#We mostly want to look at the insert size lengths and check whether there still 
#is a nucleosomal profile or not

##Libraries and environment
save.image("/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/02_analysis/DEFND_eff_env.RData")
load("/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/02_analysis/DEFND_eff_env.RData")
library(Rsamtools)
library(QuasR)


##We will analyze all four DNA files for sample #39
setwd('/omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq')
bams.files <- c('./AS-1617760-LR-80746/outs/possorted_bam.bam',
                './AS-1617768-LR-80746/outs/possorted_bam.bam',
                './AS-1617776-LR-80746/outs/possorted_bam.bam',
                './AS-1617784-LR-80746/outs/possorted_bam.bam')
scanBamHeader(bams.files[[1]])

#we need to set the parameters we're going to look for
#and extract them from all four files
para <- ScanBamParam(what = c("rname", "strand", "pos", "isize"))
bam.all <- lapply(bams.files, function(x){scanBam(x, param = para)})

#We create a histogram of all insert size lengths. To get a clear profile, 
#extreme insert sizes are excluded and the data is log2 transformed 
hist(log2(abs(bam.all[[1]][[1]]$isize[bam.all[[1]][[1]]$isize < 1000 
                                      & bam.all[[1]][[1]]$isize > -1000])), 
     breaks = 1000, 
     xlab = 'log2(insert size)', 
     main = 'insert length B6Cas neg_1',
     ylim = c(4, 10))
abline(v = c(log2(151), log2(301)), lty = 2)

#although two small peaks are visible at lengths of 151 and 301 bp, 
#the characteristic ATAC profile is gone
#there also appears to be no over-fragmentation
