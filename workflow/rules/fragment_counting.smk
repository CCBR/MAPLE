rule create_fragemented_bed:
    """
    creates a bed file from mapped bed to include only fragement sizes specified
    """
    input:
        sample_bed=rules.alignment.output.bed,
        protein_bed=config["master_bed_file"],
    envmodules:
        TOOLS["bedtools"]["version"]
    threads: getthreads("fragment_analysis")
    params:
        rname="fragment_analysis",
        sample_id='{sample_id}',
        count_nsms_script=join(WORKDIR,"scripts","COUNT_NSMS_IN_GENES.py"),
        localtmp=join(RESULTSDIR,'tmp','fragment_analysis'),
    output:
        frag_bed=join(RESULTSDIR,'03_aligned','02_bed',"{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.fragmented_bed")
    shell:
        """
        # set tmp
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        min_frag={wildcards.min_length}
        max_frag={wildcards.max_length}
        
        awk -v "min_frag=$min_frag" -v "max_frag=$max_frag" '{{ if ($3-$2 >= min_frag && $3-$2 < max_frag) print $0}}' {input.sample_bed} > $tmp_dir/mapped.bed
        bedtools intersect -wo -a $tmp_dir/mapped.bed  -b {input.protein_bed} > {output.frag_bed}
        """

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
        bed=rules.create_fragemented_bed.output.frag_bed
    envmodules:
        TOOLS["python37"]["version"],
        TOOLS["bedtools"]["version"]
    threads: getthreads("fragment_analysis")
    params:
        rname="fragment_analysis",
        sample_id='{sample_id}',
        count_nsms_script=join(WORKDIR,"scripts","COUNT_NSMS_IN_GENES.py"),
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
        
        python3 {params.count_nsms_script} {input.bed} {params.sample_id}.$min_frag-$max_frag {output.csv}
        """

rule merges_frag_csv:
    input:
        counts_list=expand(join(RESULTSDIR,'03_aligned','04_counts','{sample_id}.{sp}.{min_max_list}.InGenes.counts.csv'),sample_id=SAMPLES, sp=SPECIES,min_max_list=MIN_MAX_LIST)
    envmodules:
        TOOLS["python37"]["version"],
    threads: getthreads("merges_frag_csv")
    params:
        rname="merges_frag_csv",
        protein_bed=config["master_bed_file"],
        counts_dir=join(RESULTSDIR,'03_aligned','04_counts'),
        merged_script=join(WORKDIR,"scripts","merge_counts.py")
    output:
        csv=join(RESULTSDIR,'03_aligned','04_counts','{sample_id}.{sp}.merged.counts.csv'),
    shell:
        """
        python {params.merged_script} {params.counts_dir} {params.protein_bed} {output.csv}
        """