##This is the script to align the data from the DEFND-seq Experiment 
#all samples are from mice of the 36 wks cohort 
#39 -> CasB6-DEN
#42 -> B6Cas-DEN
#45 -> B6Cas-neg
#48 -> CasB6-neg

#connect to the cluster and check what modules are avialable
ssh s205e@bsub01.lsf.dkfz.de
module avail

####We will first analyze the RNA data as we know it####
#cellranger count needs the smaples to be named after a certain pattern. 
#since we can't modify the name, we're copying it to our directory
cd /omics/groups/OE0538/internal/shared_data/Sylviane_ATAC/DEFND_liver/seqID_43693/raw/data/250624_VH00211_445_AACWCYNHV
#39
cp ./AS-1617751-LR-80745/fastq/AS-1617751-LR-80745_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617751-LR-80745_S1_L001_R1_001.fastq.gz
cp ./AS-1617751-LR-80745/fastq/AS-1617751-LR-80745_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617751-LR-80745_S1_L001_R2_001.fastq.gz

#42
cp ./AS-1617753-LR-80745/fastq/AS-1617753-LR-80745_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617753-LR-80745_S1_L001_R1_001.fastq.gz
cp ./AS-1617753-LR-80745/fastq/AS-1617753-LR-80745_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617753-LR-80745_S1_L001_R2_001.fastq.gz

#45
cp ./AS-1617755-LR-80745/fastq/AS-1617755-LR-80745_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617755-LR-80745_S1_L001_R1_001.fastq.gz
cp ./AS-1617755-LR-80745/fastq/AS-1617755-LR-80745_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617755-LR-80745_S1_L001_R2_001.fastq.gz

#48
cp ./AS-1617757-LR-80745/fastq/AS-1617757-LR-80745_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617757-LR-80745_S1_L001_R1_001.fastq.gz
cp ./AS-1617757-LR-80745/fastq/AS-1617757-LR-80745_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/RNA/AS-1617757-LR-80745_S1_L001_R2_001.fastq.gz

#let's load the module CellRanger to analyze single-cell RNA-seq data
module load CellRanger/8.0.1

#lets do counting 
#The skripts are coded in separate files. This allows us to submit them to the 
#cluster more easily
#the refenrence transcriptome was previously downloaded from the 10X website
#39
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/01_count
chmod +x ./count_39.sh
bsub -n 20 -R "rusage[mem=150G]" -q long-debian -J "counting 39" ./count_39.sh

#42
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/01_count
chmod +x ./count_42.sh
bsub -n 10 -R "rusage[mem=80G]" -q long-debian -J "counting 42" ./count_42.sh

#45
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/01_count
chmod +x ./count_45.sh
bsub -n 20 -R "rusage[mem=90G]" -q long-debian -J "counting 45" ./count_45.sh

#48
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/01_count
chmod +x ./count_48.sh
bsub -n 20 -R "rusage[mem=90G]" -q long-debian -J "counting 48" ./count_48.sh



####Here the DNA sequences are analyzed####
#to count anything we need to get a reference genome.
#the reference genome is available at the 10X website
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq/01_count/DNA
bsub -n 10 -R "rusage[mem=30G]" -q long-debian wget  "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCm39-2024-A.tar.gz"
tar -xvzf refdata-cellranger-arc-GRCm39-2024-A.tar.gz

#cellranger-atac count needs the samples to be named after a certain pattern. since we can't modify the name, we're copying it to our directory
#62
cd /omics/groups/OE0538/internal/shared_data/Sylviane_ATAC/DEFND_liver/seqID_43694/raw/data/250624_VH00693_237_AACWFK3HV

cp ./AS-1617760-LR-80746/fastq/AS-1617760-LR-80746_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617760-LR-80746/AS-1617760-LR-80746_S1_L001_R1_001.fastq.gz
cp ./AS-1617760-LR-80746/fastq/AS-1617760-LR-80746_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617760-LR-80746/AS-1617760-LR-80746_S1_L001_R2_001.fastq.gz
cp ./AS-1617760-LR-80746/fastq/AS-1617760-LR-80746_I2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617760-LR-80746/AS-1617760-LR-80746_S1_L001_I2_001.fastq.gz

#because I don't want to get crazy, the following three lines will just be re-used for all files
cp ./AS-1617790-LR-80746/fastq/AS-1617790-LR-80746_R1.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617790-LR-80746/AS-1617790-LR-80746_S1_L001_R1_001.fastq.gz
cp ./AS-1617790-LR-80746/fastq/AS-1617790-LR-80746_R2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617790-LR-80746/AS-1617790-LR-80746_S1_L001_R2_001.fastq.gz
cp ./AS-1617790-LR-80746/fastq/AS-1617790-LR-80746_I2.fastq.gz ../../../../../../../users/sylviane_B270/projects/DEFND-seq/00_raw_QC/AS-1617790-LR-80746/AS-1617790-LR-80746_S1_L001_I2_001.fastq.gz


#let's create BAM files. We use CellRanger-ATAC/2.1.0
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq
chmod +x ./01_count/DNA_count_39.sh
bsub -n 20 -R "rusage[mem=50G]" -q long-debian -J "counting" ./01_count/DNA_count_39.sh

chmod +x ./01_count/DNA_count_39_2.sh
bsub -n 20 -R "rusage[mem=50G]" -q long-debian -J "counting" ./01_count/DNA_count_39_2.sh

chmod +x ./01_count/scripts/DNA_count_AS-1617766-LR-80746.sh
bsub -n 20 -R "rusage[mem=50G]" -q long-debian -J "counting" ./01_count/scripts/DNA_count_AS-1617766-LR-80746.sh

#This is going to be reused several times to submit the jobs separately:
chmod +x ./01_count/scripts/DNA_count_AS-1617788-LR-80746.sh
bsub -n 20 -R "rusage[mem=70G]" -q long-debian -J "counting" ./01_count/scripts/DNA_count_AS-1617788-LR-80746.sh

