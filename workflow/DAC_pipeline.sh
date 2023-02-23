
module load cutadapt
module load pear
module load bowtie
module load samtools
module load bedtools



path=$1
R1=$2
R2=$3
output=$4


## PART I ###########

## Trim adapters
cutadapt -j 32 -b file:TruSeq_and_nextera_adapters.fa -B file:TruSeq_and_nextera_adapters.fa --trim-n -m 50 -o ${output}_R1.trimmed.fastq.gz -p ${output}_R2.trimmed.fastq.gz ${path}/${R1} ${path}/${R2}

## Assemble read pairs ###########
pear -p 0.0001 -f ${output}_R1.trimmed.fastq.gz -r ${output}_R2.trimmed.fastq.gz -o ${output} -j 32

## Compress and delete unused files ###########
gzip -f ${output}.assembled.fastq
rm *discarded* *unassembled*

## Map and discard unnmaped ###########
bowtie2 -p 32 -x hg19 -U ${output}.assembled.fastq.gz -S ${output}.assembled.hg19.bam
samtools view -b -F 260 ${output}.assembled.hg19.bam > ${output}.mapped.hg19.bam
bedtools bamtobed -i ${output}.mapped.hg19.bam > ${output}.mapped.hg19.bed

## Make histogram of fragment lengths ###########
Rscript hist.r ${output}.mapped.hg19.bed ${output}.mapped.hg19.length_hist.csv


## OPTIONAL ###########
# bedtools -a ${output}.mapped.hg19.bed -b intervals_of_interests.bed > ${output}.selected.hg19.bed
# awk '{ if ($3-$2 >= 140 && $3-$2 <= 160) print $0}' ${output}.selected.hg19.bed > ${output}.selected.hg19.140-160.bed 

## PART II ###########

limit=$5
max_dist=$6

## Find fragment centers (DYADs) and make histogram (Occurrences) ###########
python3 WeigthedDYADposition.py ${output}.selected.hg19.140-160.bed ${output}.DYADs
sort -k1,1 -k2n,2 ${output}.DYADs > ${output}.DYADs.sorted
python Uniq_Position.py ${output}.DYADs.sorted ${output}.DYADs.hist

## Compute auto-correlation ###########
python DAC.py ${output}.DYADs.hist ${limit} ${max_dist} ${output}.DAC.csv

average_length=$(awk '{ SUM += ($3-$2); n++} END {print(int(SUM/n))}' intervals_of_interests.bed)
python DAC_denominator.py ${output}.DYADs.his ${limit} ${max_dist} ${average_length} ${output}.DAC.csv

## COUNT NSM FRAGMENTS IN GENES ##


awk '{ if ($3-$2 >= 120 && $3-$2 < 130) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.120-130.hg19.bed
bedtools intersect -wo -a ${output}.mapped.120-130.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.120-130.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.120-130.InGenes.hg19.bed ${output}.120-130.InGenes ${output}.120-130.InGenes.hg19.counts.csv

awk '{ if ($3-$2 >= 130 && $3-$2 < 140) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.130-140.hg19.bed 
bedtools intersect -wo -a ${output}.mapped.130-140.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.130-140.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.130-140.InGenes.hg19.bed ${output}.130-140.InGenes ${output}.130-140.InGenes.hg19.counts.csv

awk '{ if ($3-$2 >= 140 && $3-$2 < 150) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.140-150.hg19.bed 
bedtools intersect -wo -a ${output}.mapped.140-150.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.140-150.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.140-150.InGenes.hg19.bed ${output}.140-150.InGenes ${output}.140-150.InGenes.hg19.counts.csv

awk '{ if ($3-$2 >= 150 && $3-$2 < 160) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.150-160.hg19.bed 
bedtools intersect -wo -a ${output}.mapped.150-160.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.150-160.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.150-160.InGenes.hg19.bed ${output}.150-160.InGenes ${output}.150-160.InGenes.hg19.counts.csv

awk '{ if ($3-$2 >= 160 && $3-$2 < 170) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.160-170.hg19.bed 
bedtools intersect -wo -a ${output}.mapped.160-170.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.160-170.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.160-170.InGenes.hg19.bed ${output}.80-140.InGenes ${output}.160-170.InGenes.hg19.counts.csv

awk '{ if ($3-$2 >= 170 && $3-$2 < 180) print $0}' ${output}.mapped.hg19.bed > ${output}.mapped.170-180.hg19.bed 
bedtools intersect -wo -a ${output}.mapped.170-180.hg19.bed  -b hg19_protein-coding_genes.bed > ${output}.170-180.InGenes.hg19.bed
python3 COUNT_NSMS_IN_GENES.py ${output}.170-180.InGenes.hg19.bed ${output}.80-140.InGenes ${output}.170-180.InGenes.hg19.counts.csv

python merge_counts.py 






###################################
sh run --runmode=init --workdir=/data/sevillas2/ccbr1214/test


# contrast manifest can have additional contrasts
## DEFAULT is two
contrast1   contrast2
Sample1 Sample2

## three contrasts
contrast1   contrast2   contrast3
Sample1 Sample2 Sample3
