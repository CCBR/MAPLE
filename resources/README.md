##Resources

This folder, `resources/`, is meant to contain all resources necessary for running the workflow. 
- The `cluster.yaml` file can make specific resource requests to biowulf via slurm.
- The `index.yaml` file points to the paths of all Bowtie2 index files (already built).
- The `tools.yaml` file contains the version numbers of all relevent tools used within the pipeline.

# bed files generated from gene lists
## proteinCoding
## created by WG
/data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/hg19_protein-coding_genes.bed

## 1000_ALUrich
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/1000.Genes.ALU-rich.txt
/data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/1000.Genes.ALU-rich.bed

## 2500_ALUdepleted
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/2500.Genes.ALU-neg.txt
/data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/2500.Genes.ALU-depleted.bed

## GENES_2000_Length
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/_GENES_2000-Legth__Nov 29.xlsx
#create_bed /data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/GENES_2000_Length.txt /data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/GENES_2000_Length.bed

## GENES_2000_AT
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/_GENES_2000-AT__Nov 29.xlsx
#create_bed /data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/GENES_2000_AT.txt /data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/GENES_2000_AT.bed

## GENES_2000_GC
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/_GENES_2000-GC__Nov 29.xlsx
#create_bed /data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/GENES_2000_GC.txt /data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/GENES_2000_GC.bed

## GENES_2000_ALU
/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/_GENES_2000-ALU__Nov 29_sorted.xlsx
#create_bed /data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/GENES_2000_ALU.txt /data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/GENES_2000_ALU.bed