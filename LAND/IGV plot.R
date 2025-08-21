setwd('/omics/groups/OE0538/internal/users')
save.image("./sylviane_B270/projects/LAND/04_analysis/plot_environment.RData")
load("./sylviane_B270/projects/LAND/04_analysis/plot_environment.RData")

library(Gviz)
library(biomaRt)
library(AnnotationHub)
library(GenomicFeatures)

##postions of the observed gene----
#191.1e6 and 191.2e6 with 'chr1' max FC
#GAPDH: chr6:125,137,811-125,144,450
#Vinculin: chr14:20,972,745-20,985,848
#nothing in ATAC but sth in LAND: chr1 11,445,001-11,446,000
#Rn7sk: chr9:78,081,584-78,083,915

from <- 11.440e6
to <- 11.447e6
chrsm <- 'chr1'
gen <- 'mm39'
col <- c("#F4A582","#B3CDE3", "#8C96C6", "#8856A7", "#810F7C")
snames = c('ATAC', 'LAND\n5min', 'LAND\n10min', 'LAND\n15min', 'LAND\n20min')

##Reference genome----
gtrack <- GenomeAxisTrack()
#itrack <- IdeogramTrack(genome = gen, chromosome = chrsm) #to be continued

plotTracks(gtrack, from = from, to = to)

#gene annotation
#query(AnnotationHub(),c("TxDb", "Mus musculus", "mm39"))
ga.mm39 <- AnnotationHub()[['AH84139']]
geneTrack <- GeneRegionTrack(ga.mm39, genome = 'mm39', chromosome = chrsm,
                       start = from, end = to,
                       name = "Mouse Genes (mm39)", 
                       col = 'black', fill = 'grey65')


##peak data in gene range with counts----
gene.complete.gr <- GRanges(seqnames = Rle(chrsm),
                   ranges = IRanges(start = from, end = to),
                   strand = Rle('*'))
gene.gr <- tile(gene.complete.gr, width = 1e2)
gene.cts <- qCount(proj, query = gene.gr[[1]])
gene.norm <- t(gene.cts[,-1])*min(map.stats[,'mapped'])/map.stats[,'mapped']


strack <- list()
for (i in c(1:5)){
  strack[[i]] <- DataTrack(data = gene.norm[i, ], 
                           range = gene.gr, type = 'h', 
                           col = col[i], name = snames[i])
}

##plot with counts----
#this is where the plot is coded - everything is else is up
plotTracks(list(gtrack, strack[[1]], strack[[2]], strack[[3]], 
                strack[[4]], strack[[5]], geneTrack), 
           from = from, to = to, ylim = c(0.1,45), title.width = 0.5)

##data with BAM files----
bam.dir <- c('./sylviane_B270/projects/LAND/02_align/AS61/AS61_align.sorted.bam',
             './sylviane_B270/projects/LAND/02_align/AS63/AS63_align.sorted.bam',
             './sylviane_B270/projects/LAND/02_align/AS65/AS65_align.sorted.bam',
             './sylviane_B270/projects/LAND/02_align/AS67/AS67_align.sorted.bam',
             './sylviane_B270/projects/LAND/02_align/AS69/AS69_align.sorted.bam')

bamtrack <- list()

for (i in c(1:5)){
  bamtrack[[i]] 
    <- DataTrack(range = bam.dir[i], genome = 'mm39', type = 'h',
                 name = snames[i], chromosome = chrsm, window = 0.5,
                 col = col[i], 
                 transformation = function(x) {
                   x*min(map.stats[,'mapped'])/map.stats[,'mapped'][i]})
  }

plotTracks(list(gtrack, geneTrack, bamtrack[[1]], bamtrack[[2]], bamtrack[[3]], 
                bamtrack[[4]], bamtrack[[5]]),
           from = from, to = to, ylim = c(0.1, 18), title.width = 0.6)



    
