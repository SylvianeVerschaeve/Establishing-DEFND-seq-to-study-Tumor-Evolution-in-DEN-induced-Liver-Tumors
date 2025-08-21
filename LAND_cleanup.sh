### 004-6 ATAC-seq LAND Buffer
##This is the bash skript for processing our LAND ATAC-seq data from 004-6. We will do QC, trimming, alignment, peak calling and finally prepare the data for IGV analysis.

#connect to the server 
ssh s205e@bsub01.lsf.dkfz.de
## As a first step: Lets copy our data to our folder 00_raw. It's important to use the -r function if we want to copy all folders. The * allows us to copy multiple folders simultneously
cp -r /omics/groups/OE0538/internal/shared_data/Sylviane_ATAC/data/250430_VH00211_432_AAGHYCKM5/*79912 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/00_raw

##FastQC
#The raw data already has some FastQC files with them. To make it more concise we want to use MultiQC. Unfortunately the python version is not the right one in the environment we've used so far. So we need to create a new environment with the right python version:
micromamba create -p /omics/groups/OE0538/internal/users/sylviane_B270/environments/multiqc_env python=3.11
#Since it is new, we also have to initialize it bfeore we activate it:
micromamba init -s bash
source ~/.bashrc
micromamba activate /omics/groups/OE0538/internal/users/sylviane_B270/environments/multiqc_env
# now we install the multiqc package and test it
micromamba install multiqc
multiqc --help
#This is how we use it (I love the output):
multiqc /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/00_raw/*79912/fastq/*fastqc.zip -o /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/00_raw/QC_report.html 
#Here we see that adaptors clearly need to get trimmed. So let's do that!
#Afterwards, we deactivate the active environment
micromamba deactivate

##Trimming with fastp
#First, let's go back to our environment and activate our module. 
micromamba activate /omics/groups/OE0538/internal/users/sylviane_B270/environments/mamba1
module load fastp/0.23.4-GCC-14.1.0
#and use it to trim our samples. For this we first create new folder strucutres and then try to code a loop structure that will make things easier.
mkdir /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531361-LR-79912/fastq

for readOne in /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/00_raw/*79912/fastq/*_R1.fastq.gz
  do
  readTwo="${readOne/R1/R2}"
  outOne="${readOne/00_raw/01_trim}"
  outTwo="${readTwo/00_raw/01_trim}"
  fastp -i $readOne -o $outOne -I $readTwo -O $outTwo
  done

#For quality control we will create FastQC files 
module load Temurin/21.0.3_9
module load FastQC/0.12.1 
fastqc /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/*79912/fastq/*.fastq.gz

#And again we create a MultiQC report. Remember to activate the environment we created earlier!
micromamba activate /omics/groups/OE0538/internal/users/sylviane_B270/environments/multiqc_env
multiqc /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/*79912/fastq/*fastqc.zip -o /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/QC_trim.html
micromamba deactivate

##Alignment with bowtie2
#First I want to test Paul's code
module load Bowtie2/2.5.4-GCC-14.1.0
module load SAMtools/1.20-GCC-14.1.0
bowtie2 -p 10 --end-to-end -x /omics/groups/OE0538/internal/users/sylviane_B270/genomes/mm39/mm39_index -1 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531361-LR-79912/fastq/AS-1531361-LR-79912_R1.fastq.gz -2 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531361-LR-79912/fastq/AS-1531361-LR-79912_R2.fastq.gz | samtools view -b - > /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/AS61/AS61_align.bam

#and now see whether I can submit code to the server:
#before we can use the cluster, we need to first connect with it. For alternatives see https://wiki.odcf.dkfz.de/pub/cluster/quickstart/basics. Confirm and type in your usual s205e password:
ssh s205e@bsub02.lsf.dkfz.de
#be careful. The packages available in the cluster are different from the ones outside. 
module load bowtie2/2.3.5.1
module load samtools/1.9
bsub -n 10 -R "rusage[mem=30G]" -q long -J "align_ATAC1_LAND" "bowtie2 -p 10 --end-to-end -x /omics/groups/OE0538/internal/users/sylviane_B270/genomes/mm39/mm39_index -1 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531361-LR-79912/fastq/AS-1531361-LR-79912_R1.fastq.gz -2 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531361-LR-79912/fastq/AS-1531361-LR-79912_R2.fastq.gz | samtools view -b - > /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/AS61/AS61_align.bam"

#If this works out we can try to do it in an associated skript:
chmod +x /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/align_single.sh
bsub -n 10 -R "rusage[mem=30G]" -q long -J "align_ATAC2_LAND" /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/align_single.sh

#If this is successful we can use a loop structure to do all the remaining alignments at once:
chmod +x /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/align_all.sh
bsub -n 10 -R "rusage[mem=30G]" -q long -J "align_ATAC3_LAND" /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/align_all.sh

#using the following command we can check what jobs are currently running:
bjobs -w

##summary of output:
#AS61: 98.30% overall alignment rate
#AS63: 98.68% overall alignment rate
#AS65: 98.64% overall alignment rate
#AS67: 98.59% overall alignment rate
#AS69: 98.57% overall alignment rate

###Peak calling with macs2
module load macs2/2.1.2.1
#first for the AS 61 file, to test it. To not have such long directories, we move to the basic project one:
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND
bsub -n 10 -R "rusage[mem=10G]" -q long -J "peak_AS61" "macs2 callpeak -t ./02_align/AS61/AS61_align.bam -f BAMPE -g 'mm' --outdir ./03_peak/AS61 -n AS61_macs2"

#now we can repeat it for all the other samples using a skript 
chmod +x /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/03_peak/LAND_peak.sh
bsub -n 10 -R "rusage[mem=10G]" -q long -J "peak_All" /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/03_peak/LAND_peak.sh

#AS69 was lagging behind, but here it is:
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND
bsub -n 10 -R "rusage[mem=10G]" -q long -J "peak_AS69" "macs2 callpeak -t ./02_align/AS69/AS69_align.bam -f BAMPE -g 'mm' --outdir ./03_peak/AS69 -n AS69_macs2"

##preparation for IGV data manipulation
#for this we need to sort the bam file and create .bai files. Again first for the AS61 file
#module load samtools/1.9
#cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND
bsub -n 10 -R "rusage[mem=10G]" -q long -J "sort_A61" "samtools sort ./02_align/AS61/AS61_align.bam -o ./02_align/AS61/AS61_align.sorted.bam"
bsub -n 10 -R "rusage[mem=10G]" -q long -J "index_AS61" "samtools index ./02_align/AS61/AS61_align.sorted.bam"
#for AS69
bsub -n 10 -R "rusage[mem=10G]" -q long -J "sort_A69" "samtools sort ./02_align/AS69/AS69_align.bam -o ./02_align/AS69/AS69_align.sorted.bam"

#and now for all the others
chmod +x /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/03_peak/LAND_IGV_proc.sh
bsub -n 10 -R "rusage[mem=10G]" -q long -J "sort_index_all" /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/03_peak/LAND_IGV_proc.sh
#the sorting did not work so were doing it the simple way:
bsub -n 10 -R "rusage[mem=10G]" -q long -J "index_AS63" "samtools index ./02_align/AS63/AS63_align.sorted.bam"
bsub -n 10 -R "rusage[mem=10G]" -q long -J "index_AS65" "samtools index ./02_align/AS65/AS65_align.sorted.bam"
bsub -n 10 -R "rusage[mem=10G]" -q long -J "index_AS67" "samtools index ./02_align/AS67/AS67_align.sorted.bam"
bsub -n 10 -R "rusage[mem=10G]" -q long -J "index_AS69" "samtools index ./02_align/AS69/AS69_align.sorted.bam"
###And I think the rest of the actual analysis will be done in R...


#we need the cluster real quick
chmod + /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/04_analysis/LAND_readcount_matrix.sh
bsub -n 10 -R "rusage[mem=30G]" -q long -J "matrix_R"/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/04_analysis/LAND_readcount_matrix.sh
