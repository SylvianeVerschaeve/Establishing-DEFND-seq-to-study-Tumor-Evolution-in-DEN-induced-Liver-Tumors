module load macs2/2.1.2.1
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/LAND
for nr in 63 65 67; do
macs2 callpeak -t ./02_align/AS${nr}/AS${nr}_bowtie2.bam -f BAMPE --keep-dup all --nomodel --shift -100 --extsize 200 -g 'mm' --outdir ./03_peak/AS${nr} -n AS${nr}_macs2 &
done
