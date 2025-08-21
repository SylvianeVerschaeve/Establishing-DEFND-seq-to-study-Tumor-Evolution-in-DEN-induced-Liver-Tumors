#!/bin/bash
module load CellRanger-ATAC/2.1.0
cd /omics/groups/OE0538/internal/users/sylviane_B270/projects/DEFND-seq
cellranger-atac count --id AS-1617762-LR-80746 --reference ./01_count/DNA/refdata-cellranger-arc-GRCm39-2024-A --fastqs ./00_raw_QC/AS-1617762-LR-80746 --sample AS-1617762-LR-80746 --chemistry=ARC-v1

