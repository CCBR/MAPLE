rule fragment_analysis:
    '''
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
    '''
    input:
        sample_bed=rules.alignment.output.bed,
        protein_bed=config["master_bed_file"],
    envmodules:
        TOOLS["python37"]["version"],
        TOOLS["bedtools"]["version"]
    threads: getthreads("fragment_analysis")
    params:
        rname="fragment_analysis",
        sample_id='{sample_id}',
        dac_script=join(WORKDIR,"scripts","COUNT_NSMS_IN_GENES.py"),
        localtmp=join(RESULTSDIR,'tmp','fragment_analysis'),
    output:
        csv=join(RESULTSDIR,'03_aligned','04_counts','{sample_id}.{species}.{min_max_list}.InGenes.counts.csv'),
    shell:
        """
        # set tmp
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        min_frag=`echo {wildcards.min_max_list} | cut -f1 -d"_"`
        max_frag=`echo {wildcards.min_max_list} | cut -f2 -d"_"`
        
        awk '{{ if ($3-$2 >= $min_frag && $3-$2 < $max_frag) print $0}}' {input.sample_bed} > $tmp_dir/mapped.bed
        bedtools intersect -wo -a $tmp_dir/mapped.bed  -b {input.protein_bed} > $tmp_dir/InGenes.bed
        python3 {params.dac_script} $tmp_dir/InGenes.bed {params.sample_id}.$min_frag-$max_frag {output.csv}
        """
