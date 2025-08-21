#!/bin/bash
module load bowtie2/2.3.5.1
module load samtools/1.9
bowtie2 -p 10 --end-to-end -x /omics/groups/OE0538/internal/users/sylviane_B270/genomes/mm39/mm39_index -1 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531369-LR-79912/fastq/AS-1531369-LR-79912_R1.fastq.gz -2 /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-1531369-LR-79912/fastq/AS-1531369-LR-79912_R2.fastq.gz | samtools view -b - > /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/AS69/AS69_align.bam
