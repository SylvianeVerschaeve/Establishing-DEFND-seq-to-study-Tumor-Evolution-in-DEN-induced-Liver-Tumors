##This is the file that includes the analysis of the ATAC, LAND 5, 10, 15 and 
#20 min samples. 
#The goal of this script is to figure out whether the LAND buffer is able 
#to dissolve the nucleosomes and cause a more equal distribution of the reads

##Environment and Libraries
setwd('/omics/groups/OE0538/internal/users')
save.image("./sylviane_B270/projects/LAND/04_analysis/LAND_environment_2.RData")
load("./sylviane_B270/projects/LAND/04_analysis/LAND_environment_2.RData")
library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer) 
library(readxl)
library(ggplot2)
library(GGally)
library(KernSmooth)
library(BSgenome.Mmusculus.UCSC.mm39)
library(QuasR)
library(txdbmaker)
library(tidyr)

###Peaks###----
#peaks are one of the strongest indicators of whether the LAND buffer worked. 
#A typical ATAC shows peaks in highly accessible regions such as promoter sites.
#LAND treated samples should loose these peaks, so have a smaller proportion of
#reads in those regions compared to the ATAC sample

##Peak Ranges
#get the peak files and convert them to a GRanges Object
peak.files <- list.files('./sylviane_B270/projects/LAND/03_peak', 
                         pattern='.narrowPeak$', 
                         recursive=TRUE, 
                         full.names = TRUE)
pk.gr <- lapply(peak.files, function(x){import(x)})

#we need to remove the blacklisted regions and the non genomic DNA regions 
#to create a GRanges file with all the peaks we want to include in the analysis
blckl <- import("../shared_data/excluderanges39_boyle.bed")
pk.bl.gr <- lapply(pk.gr, function(x){subsetByOverlaps(x, blckl, 
                                                       invert = TRUE)})
chr.gr <- 
  GRanges(seqnames = 
            Rle(names(seqlengths(BSgenome.Mmusculus.UCSC.mm39)[1:21])),
          strand = 
            Rle('*'),
          range = 
            IRanges(start = 1,
                    end = seqlengths(BSgenome.Mmusculus.UCSC.mm39)[1:21]))
pk.chr.gr <- lapply(pk.bl.gr, function(x){subsetByOverlaps(x, chr.gr)})


##Number of Reads in Peaks
#get the sequencing BAM files, combine them with sample names and save them 
#into a text file
bam.files <- list.files('./sylviane_B270/projects/LAND/02_align', 
                        pattern='.sorted.bam$', 
                        recursive=TRUE, 
                        full.names = TRUE)
atac.df <- data.frame('FileName'=bam.files, 
                      'SampleName'=c('ATAC','LAND5min','LAND10min','LAND15min',
                                     'LAND20min'), 
                      stringsAsFactors = FALSE)
write.table(atac.df, 
            file='sylviane_B270/projects/LAND/samples.txt', 
            sep='\t', 
            quote=FALSE, 
            row.names=FALSE)


#checks whether all reads in bam files are aligned to mm39 genome
#map.stats gives specific numbers
proj <- qAlign(sampleFile = 'sylviane_B270/projects/LAND/samples.txt',
               genome='BSgenome.Mmusculus.UCSC.mm39',
               paired='fr')
map.stats <- alignmentStats(proj)

#quantify reads of proj in query region 
#-> How many reads align to each noLAND peak in each of the samples?
atac.cts <- qCount(proj, query = pk.chr.gr[[1]])

#normalize the reads to the lowest total read count. 
#Additionally we normalize to the peak width
atac.norm <- 
  t(t(atac.cts[,-1]) * 
      min(map.stats[,'mapped'])/map.stats[,'mapped'])/atac.cts[,1]*1000
atac.norm.wo.length <- 
  t(t(atac.cts[,-1]) * min(map.stats[,'mapped'])/map.stats[,'mapped'])

##peak vizualization
#pairs plot of the read counts per peak region
colnames(atac.norm) <- 
  c('ATAC', 'LAND 5 min', 'LAND 10 min', 'LAND 15 min', 'LAND 20 min')
pairs(log2(atac.norm + 1), pch = ".", xlim = c(0, 20), ylim = c(0,20))

#scatterplots of LAND10 min vs ATAC and LAND5 min
smoothScatter(log2(atac.norm[,1]+1), log2(atac.norm[,3]+1), asp = 1, 
              xlab = "log2(#read ATAC)", ylab = "log2(#reads LAND 10 min)",
              xlim = c(0, 20), ylim = c(0,20)); abline(0,1)

smoothScatter(log2(atac.norm[,2]+1), log2(atac.norm[,3]+1), asp = 1, 
              xlab = "log2(#read LAND5 min)", ylab = "log2(#reads LAND10 min)",
              xlim = c(0, 20), ylim = c(0,20)); abline(0,1)

##Nr of peaks in peak regions vs library size
#total number of reads in peak-regions
atac.cts.pk <- colSums(atac.cts[,-1])
atac.cts.tot <- map.stats[, "mapped"]

#ratio  #reads in peaks region/ #reads in total
atac.cts.rat <- atac.cts.pk/atac.cts.tot
atac.cts.rat
atac.cts.ov <- data.frame(rbind(atac.cts.pk, atac.cts.tot, atac.cts.rat*100))
rownames(atac.cts.ov) = c("#read in peaks", '#reads in total', 'ratio [%]')

barplot(atac.cts.rat*100, ylab = 'fraction of reads in peaks [%]', las = 2,
        border = c("#F4A582","#B3CDE3", "#8C96C6", "#8856A7", "#810F7C"), 
        col = 'white',
        names.arg = c('ATAC', 'LAND\n5 min', 'LAND\n10 min', 'LAND\n15 min', 
                      'LAND\n20 min'))

###gaps between peaks----
#opposed to the peak regions we're also going to look at the gap regions 
#between those peaks. In the ideal case we will see an increase in read counts 
#in some of those regions although we have to consider that not all the gaps 
#between the peaks have low coverage, so some of them might still see a decrease
#in coverage after the LAND treatment

##gaps Ranges
#create GRanges Object with gap-regions and take out the blacklisted regions
gp.gr <- gaps(pk.gr[[1]])
gp.bl.gr <- setdiff(gp.gr, blckl)

##Read counts in gaps
#count read in gaps without blacklisted regions and normalize on library size 
#and normalization on gap length
gp.ct <- qCount(proj, query = gp.bl.gr)
gp.norm <- t(t(gp.ct[,-1])*min(map.stats[,'mapped'])/map.stats[,'mapped'])

##data visualization
#compute the delta between ATAC and  LAND
gp.delta <- log2(gp.norm[,c(2:5)]+1)-log2(gp.norm[,1] +1)
smoothScatter(gp.delta[,1], pch = '.', 
              main = "delta LAND5 - ATAC gap regions", 
              ylab = "delta")

#MA-plot
mean.cpm.atac.LAND10 <- rowMeans(gp.norm[,c(1, 3)])
fc.cpm.atac.LAND10 <- gp.norm[,3] / gp.norm[,1]
smoothScatter(log2(mean.cpm.atac.LAND10),log2(fc.cpm.atac.LAND10), pch = '.', 
              xlab = 'log2(mean CPM)', ylab = 'log2(LAND 10 min / ATAC)')
abline(h=0, lty = 2)

#MA of LAND samples
mean.cpm.LAND5.LAND10 <- rowMeans(gp.norm[,2:3])
fc.cpm.LAND5.LAND10 <- gp.norm[,2]/gp.norm[,3]
smoothScatter(log2(mean.cpm.LAND5.LAND10),log2(fc.cpm.LAND5.LAND10), 
              pch = '.', 
              main = 'MA-plot of LAND 5 and 10 min', 
              xlab = 'log2(mean CPM)', 
              ylab = 'log2(#read LAND5/ #read LAND10)')
abline(h=0, lty = 2)


###insert size lengths----
#We're mostly going to work with the RSamtools package. 
#Let's first use it to look at out BAM files
scanBamHeader(bam.files[[1]])

#Extract the important information (isize)
#this gives us the available fields
scanBamWhat()

##ATAC-sample
#Extract insert size lengths from BAM files
#let's select the params we want to look at. 
para <- ScanBamParam(what = c("rname", "strand", "pos", "isize"))
bam.atac <- scanBam(bam.files[[1]], param = para)

atac.sz <- data.frame(seqname = bam.atac[[1]]$rname,
                      strand = bam.atac[[1]]$strand,
                      pos5 = bam.atac[[1]]$pos,
                      size = abs(bam.atac[[1]]$isize))

#we need to take out the blacklisted regions too blckl
#remove all lines with NA
atac.sz.nna <- atac.sz[!is.na(atac.sz[,"pos5"]) & !is.na(atac.sz[,"size"]),]

#create GRanges Object
atac.sz.gr.blckl <- GRanges(seqnames = atac.sz.nna$seqname,
                            strand = atac.sz.nna$strand,
                            range = IRanges(start = atac.sz.nna$pos5, 
                                            end = atac.sz.nna$pos5 + 
                                              atac.sz.nna$size),
                            size = abs(atac.sz.nna$size))

atac.sz.gr <- subsetByOverlaps(atac.sz.gr.blckl, blckl, invert = TRUE)


#turns out we've got some non gDNA stuff in there - especially the mitochondrial 
#DNA is bad. 
atac.sz.gen.gr <- subsetByOverlaps(atac.sz.gr, chr.gr)


##LAND5 sample
#let's select the parameters we want to look at.
para <- ScanBamParam(what = c("rname", "strand", "pos", "isize"))
bam.LAND5 <- scanBam(bam.files[[2]], param = para)

LAND5.sz <- data.frame(seqname = bam.LAND5[[1]]$rname,
                       strand = bam.LAND5[[1]]$strand,
                       pos5 = bam.LAND5[[1]]$pos,
                       size = abs(bam.LAND5[[1]]$isize))

#we need to take out the blacklisted regions to remove all lines with NA
LAND5.sz.nna <- LAND5.sz[!is.na(LAND5.sz[,"pos5"]) & !is.na(LAND5.sz[,"size"]),]

#create GRanges Object
LAND5.sz.gr.blckl <- GRanges(seqnames = LAND5.sz.nna$seqname,
                             strand = LAND5.sz.nna$strand,
                             range = IRanges(start = LAND5.sz.nna$pos5, 
                                             end = LAND5.sz.nna$pos5 + 
                                               LAND5.sz.nna$size),
                             size = abs(LAND5.sz.nna$size))

LAND5.sz.gr <- subsetByOverlaps(LAND5.sz.gr.blckl, blckl, invert = TRUE)


#remove non genomic DNA
LAND5.sz.gen.gr <- subsetByOverlaps(LAND5.sz.gr, chr.gr)

##basic data visualization
#plot insert size frequency of the ATAC sample
hist(log2(atac.sz.gen.gr$size +1), 
     breaks = 3000, 
     xlim = c(4, 10), 
     xlab = "insert size", 
     main = "insert size distribution ATAC")
hist(log2(atac.sz.gen.gr$size[atac.sz.gen.gr$size > 100] +1), 
     breaks = 3000, 
     xlim = c(6.5, 8), 
     xlab = "insert size", 
     main = "insert size distribution ATAC")
abline(v = c(log2(150 +1), log2(300+1)))

#plot insert size frequency of LAND5 sample
hist(log2(LAND5.sz.gen.gr$size+1), 
     breaks = 5000, 
     xlim= c(4,10), 
     xlab = "insert size", 
     main = "insert size distribution LAND 5 min")
hist(log2(LAND5.sz.gen.gr$size[LAND5.sz.gen.gr$size > 100] +1), 
     breaks = 3000, 
     xlim = c(6.5, 8), 
     xlab = "insert size", 
     main = "insert size distribution LAND 5 min")
abline(v = c(log2(150+1), log2(300+1)))

##Nucleosomal profile plot
#Loess regression plots ATAC
#extracting values out of histogram  
ins.sz.hist.atac <- 
  hist(log2(atac.sz.gen.gr$size[atac.sz.gen.gr$size > (2^4.95 -1) 
                                & atac.sz.gen.gr$size < (2^9.2 -1)]), 
       breaks = 3000, 
       xlim = c(4.95, 9.2), 
       xlab = "insert size", 
       main = "insert size distribution ATAC")
lr.ct.atac <- ins.sz.hist$counts
lr.mid.atac <- ins.sz.hist$mids

mids.atac <- lr.mid.atac[which(lr.mid.atac > 5.05 & lr.mid.atac < 8)]
cts.atac <- lr.ct[which(lr.mid.atac > 5.05 & lr.mid.atac < 8)]
cts.0.atac <- cts.atac; cts.0.atac[which(cts.0.atac == 0)] <- NA
mids.0.atac <- mids.atac; mids.0.atac[which(cts.0.atac == 0)] <- NA

#create plot with loess regression
lr.df.atac <- data.frame(mids.0.atac = mids.0.atac,
                         cts.0.atac = cts.0.atac)
lr.ls.atac <- loess(cts.0.atac ~ mids.0.atac, df = lr.df.atac, span = 0.2)
lr.smth.atac <- predict(lr.ls.atac)

mids.atac.0.rm <- mids.0.atac[!is.na(cts.0.atac)]
plot(mids.0.atac, cts.0.atac, 
     pch = '.', 
     xlab = 'log2(framgent length)', 
     ylab = 'fragments')
lines(mids.atac.0.rm, lr.smth.atac, col = 'red')

#delta between regression and actual values
cts.atac.0.rm <- cts.0.atac[!is.na(cts.0.atac)]
lr.delta.atac <- cts.atac.0.rm - lr.smth.atac

plot(mids.atac.0.rm, lr.delta.atac, 
     pch = '.', 
     xlab = 'log2(framgent length)',
     ylab = 'log2(counts) - loess')
abline(h=0, lty = 2)

#fit a curve
df.deltadk.ls.atac <- data.frame(mids.atac.0.rm = mids.atac.0.rm,
                                 lr.delta.atac = lr.delta.atac)
delta.ls.atac <- smooth.spline(mids.atac.0.rm, lr.delta.atac)
delta.smth.atac <- predict(delta.ls.atac)

plot(mids.atac.0.rm, lr.delta.atac, 
     pch = '.', col = 'red', 
     xlab = 'log2(framgent length)', 
     ylab = 'log2(counts) - loess', 
     xlim = c(5.05, 8.0))
lines(mids.atac.0.rm, delta.smth.atac$y, col = 'red')

#loess regression plots LAND5
#extract values from histograms
ins.sz.hist.LAND5 <- 
  hist(log2(LAND5.sz.gen.gr$size[LAND5.sz.gen.gr$size > (2^4.95 -1) 
                                 & LAND5.sz.gen.gr$size < (2^9.2 -1)]), 
       breaks = 3000, 
       xlim = c(4.95, 9.2), 
       xlab = "insert size", 
       main = "insert size distribution LAND 5 min")

lr.cts.5 <- ins.sz.hist.LAND5$counts
lr.mid.5 <- ins.sz.hist.LAND5$mids

mids.5 <- lr.mid.5[which(lr.mid.5 > 5.05 & lr.mid.5 < 8)]
cts.5 <- lr.cts.5[which(lr.mid.5 > 5.05 & lr.mid.5 < 8)]
cts.5.0.rm <- cts.5[which(cts.5 !=0)]
mids.5.0.rm <- mids.5[which(cts.5 != 0)]

plot(mids.5.0.rm, cts.5.0.rm, 
     pch = '.', 
     xlab = 'log2(fragment length)', 
     ylab = '#fragments')

#do loess regression
lr.df.5 <- data.frame(mids.5.0.rm <- mids.5.0.rm,
                      cts.5.0.rm <- cts.5.0.rm)
lr.ls.5 <- loess(cts.5.0.rm ~ mids.5.0.rm, df = lr.df.5, span = 0.25)
lr.smth.5 <- predict(lr.ls.5)

plot(mids.5.0.rm, cts.5.0.rm, 
     pch = '.', 
     xlab = 'log2(fragment length)', 
     ylab = '#fragments')
lines(mids.5.0.rm, lr.smth.5, col = 'blue')

#compute deltas
lr.delta.5 <- cts.5.0.rm - lr.smth.5

plot(mids.5.0.rm, lr.delta.5, pch = '.', col = 'blue')

#fit a smooth line
df.delta.ls.5 <- data.frame(mids.5.0.rm = mids.5.0.rm, 
                            lr.delta.5 = lr.delta.5)
delta.ls.5 <- smooth.spline(mids.5.0.rm, lr.delta.5)
delta.smth.5 <- predict(delta.ls.5)

plot(mids.5.0.rm, lr.delta.5, pch = '.', xlim = c(5.25, 8.0))
lines(mids.5.0.rm, delta.smth.5$y, col = 'blue')
abline(h=0, lty = 2)

#lets look at the interesting areas (150 bp - 800 bp)
plot(mids.5.0.rm[mids.5.0.rm > log2(130)], lr.delta.5[mids.5.0.rm > log2(130)], 
     pch = '.',
     ylab = 'residuals log2(#reads)', 
     xlab = 'fragment length')
lines(mids.5.0.rm[mids.5.0.rm > log2(130)], 
      delta.smth.5$y[mids.5.0.rm > log2(130)], col = 'blue')
lines(mids.atac.0.rm[mids.5.0.rm > log2(130)], 
      delta.smth.atac$y[mids.5.0.rm > log2(130)], col = 'red')
abline(h=0, v = c(log2(151)), lty = 2)


###table overview with all counts----
map.stats[,'mapped'] #library size, whole read counts
atac.norm.wo.length # peak counts
gp.norm # gap counts

blckl.cts <- qCount(proj, query = blckl) #counts black listed regions
blckl.norm <- t(t(blckl.cts[,-1]) * min(map.stats[,'mapped'])/map.stats[,'mapped'])

chr.cts <- qCount(proj, query = chr.gr) #counts non auto or gonomsomes
nchr.cts <- map.stats[,'mapped'] - colSums(chr.cts[,-1])
nchr.norm <- t(t(nchr.cts) * min(map.stats[,'mapped'])/map.stats[,'mapped'])

table.all <- cbind(map.stats[,'mapped'], colSums(atac.norm.wo.length), colSums(gp.norm), colSums(blckl.norm), colSums(t(nchr.norm)))
colnames(table.all) <- c("library size", "#reads in peaks", "#reads in gaps", "#blacklisted reads", "#reads in non chromosomes")
rownames(table.all) <- c("ATAC", "LAND 5 min", "LAND 10 min", "LAND 15 min", "LAND 20 min")


###read distribution over the genome - genome binning----
#create binned chromosome 1 reference without blacklisted regions
starts = seq(1, 195154279, by = 1000)
ends = starts + 999
chr1.bin <- GRanges(seqnames = Rle('chr1'),
                    ranges = IRanges(start = starts, end = ends),
                    strand = Rle('*'))

chr1.bin.gl <- GRangesList(chr1.bin)
chr.bin.bl <- lapply(chr1.bin.gl, 
                     function(x){subsetByOverlaps(x, blckl, invert = TRUE)})

#count the reads of all samples in the chr1 regions
#count reads of proj aligned in the regions of chr1.bin
bin.cts <- qCount(proj, query = chr1.bin.bl[[1]])

#normalize it
bin.norm <- t(t(bin.cts[,-1])*min(map.stats[,'mapped'])/map.stats[,'mapped'])

#histogram
plot(log2(bin.norm[,1]+100), 
     pch = '.', xlab = "1k bp", 
     ylab = 'log2(norm #reads +1)', 
     main = 'read distribution ATAC')
plot(log2(bin.norm[,5]+100), 
     pch = '.', 
     xlab = "1k bp", 
     ylab = 'log2(norm #reads +1)', 
     main = 'read distribution LAND 20min')

#delta LAND and ATAC
bin.del <- log2((bin.norm[,c(2:5)]+100)/(bin.norm[,1] + 100))
plot(bin.del[,2], 
     pch = '.', 
     xlab = "1k bp", 
     ylab = 'delta log2 LAND10 ATAC', 
     main ='delta')

rownames(bin.del) <- paste0(start(chr1.bin.bl[[1]]), '-', end(chr1.bin.bl[[1]]))

bin.del.srt.5 <- bin.del[order(bin.del[,1], decreasing = TRUE), 1]
head(bin.del)
  
#plot with vertical lines
plot(bin.del[,1], 
     pch = '.', 
     xlab = "1k bp", 
     ylab = 'delta log2 ATAC and LAND10')
for (i in c(1:4)){abline(v = which(bin.del == bin.del.srt.5[i]), 
                         lty = 2, col = 'red')}
for (j in c((length(bin.del.srt.5)-3):length(bin.del.srt.5))){
    abline(v = which(bin.del == bin.del.srt.5[j]), 
           lty = 2, col = 'blue')}


###expected read distribution with LAND treatment----
wgs.bam <- '/omics/groups/OE0538/internal/shared_data/Sylviane_ATAC/AS-1132234_B6_CAST_EiJ_sortedByCoord_markDup_haplotagged.bam'
wgs.df <- data.frame('FileName' = wgs.bam, 'SampleName' = "whole_genome")
write.table(wgs.df, 
            file = '/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/wgs.txt', 
            sep = '\t', 
            quote = FALSE, 
            row.names = FALSE )

wgs <- qAlign(sampleFile = '/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/wgs.txt', 
              genome='BSgenome.Mmusculus.UCSC.mm39',
              paired='fr')

chr1.wg.ct <- qCount(wgs, query = chr1.bin.bl[[1]])
wgs.stats <- alignmentStats(wgs)

#create a matrix with all counts
chr1.wgs.all <- cbind(bin.cts[,-1], chr1.wg.ct[,2])
rownames(chr1.wgs.all) <- 
  paste0(start(chr1.bin.bl[[1]]), '-', end(chr1.bin.bl[[1]]))
colnames(chr1.wgs.all) <- 
  c("ATAC", "LAND5min", "LAND10min", "LAND15min", "LAND20min", "WGS")
wgs.all.stats <- colSums(chr1.wgs.all)


#norm to library size
chr1.wgs.norm <- t(t(chr1.wgs.all)*min(wgs.all.stats)/wgs.all.stats)


#plot wgs alignment
smoothScatter(log2(chr1.wgs.norm[,6] +1), 
              pch = '.', 
              xlab = "1k bp", 
              ylab = 'log2(norm #reads +1)', 
              main = 'ATAC whole genome sequencing')

#plot comparison to ATAC
chr1.wgs.atac <- log2((chr1.wgs.norm[,1] + 1)/(chr1.wgs.norm[,6] + 1))
smoothScatter(chr1.wgs.atac, 
              pch = '.', 
              xlab = "1k bp", 
              ylab = 'delta', 
              main = 'delta(ATAC-wgs)')
    
head(chr1.wgs.norm)
#plot comparison to LAND10
chr1.wgs.LAND10 <- log2((chr1.wgs.norm[,3] + 100)/(chr1.wgs.norm[,6] + 100))
smoothScatter(chr1.wgs.LAND10, 
              pch = '.', 
              xlab = "1k bp", 
              ylab = 'delta', 
              main = 'delta(LAND10-wgs)')


plot(log2(chr1.wgs.norm[,1] + 1), log2(chr1.wgs.norm[,6] + 1), 
     pch = '.', main = 'ATAC WGS')
plot(log2(chr1.wgs.norm[,3] + 1), log2(chr1.wgs.norm[,6] + 1), 
     pch = '.', main = 'LAND10 WGS')

pairs(log2(chr1.wgs.norm + 1), pch = ".", xlim = c(0, 10), ylim = c(0,10))

#compute CV = standard deviation/mean (sigma/mü)
sd.chr1 <- apply(chr1.wgs.norm, 2, function(x){sd(log2(x+1))})
mean.chr1 <- apply(chr1.wgs.norm, 2, function(x){mean(log2(x+1))})
cv.chr1 <- sd.chr1/mean.chr1*100
cv.chr1


##ECDF----
#chromosome 1
head(chr1.wgs.norm)
plot(ecdf(chr1.wgs.norm[,1]), main = 'ATAC', xlab = 'coverage')
plot(ecdf(chr1.wgs.norm[,3]), main = 'LAND10', xlab = 'coverage')
plot(ecdf(chr1.wgs.norm[,6]), main = 'WGS', xlab = 'coverage')

plot(ecdf(log2(chr1.wgs.norm[,1]+1)), main = 'log2 ATAC', xlab = 'coverage')
plot(ecdf(log2(chr1.wgs.norm[,3]+1)), main = 'log2 LAND10', xlab = 'coverage')
plot(ecdf(log2(chr1.wgs.norm[,6]+1)), main = 'log2 WGS', xlab = 'coverage')

#whole genome
plot(ecdf(chr.wgs.comb[,3]), main = 'ATAC', xlab = 'coverage')
plot(ecdf(chr.wgs.comb[,5]), main = 'LAND10', xlab = 'coverage')
plot(ecdf(chr.wgs.comb[,8]), main = 'WGS', xlab = 'coverage')

plot(ecdf(log2(chr.wgs.comb[,3]+1)), main = 'log2 ATAC', xlab = 'coverage')
plot(ecdf(log2(chr.wgs.comb[,5]+1)), main = 'log2 LAND10', xlab = 'coverage')
plot(ecdf(log2(chr.wgs.comb[,8]+1)), main = 'log2 WGS', xlab = 'coverage')

