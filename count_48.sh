#!/bin/bash
module load CellRanger/8.0.1
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq
cellranger count --id AS-626600-LR-56594_ct --transcriptome ./01_count/refdata-gex-GRCm39-2024-A --create-bam true --fastqs ./00_raw_QC --sample AS-626600-LR-56594 --chemistry=ARC-v1

