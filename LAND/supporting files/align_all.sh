module load bowtie2/2.3.5.1
module load samtools/1.9
ref="/omics/groups/OE0538/internal/users/sylviane_B270/genomes/mm39/mm39_index"
for nr in 65 67 69; do
read1="/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-15313${nr}-LR-79912/fastq/*R1.fastq.gz"
read2="/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/01_trim/AS-15313${nr}-LR-79912/fastq/*R2.fastq.gz"
out="/omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND/02_align/AS${nr}/AS${nr}_bowtie2.bam"
bowtie2 -p 10 --end-to-end -x $ref -1 $read1 -2 $read2 | samtools view -b - > $out
done
