module load samtools/1.9
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND
for nr in 63 65 67; do
samtools sort ./02_align/AS${nr}/AS${nr}_bowtie2.bam -o ./02_align/AS${nr}/AS${nr}_align.sorted.bam
samtools index ./02_align/AS${nr}_align.sorted.bam &
done
