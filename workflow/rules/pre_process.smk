rule trim_adaptors:
    '''
    Read samples.tsv to determine fastq files, trim adaptors

    cutadapt -j 32 -b file:TruSeq_and_nextera_adaptors.fa -B file:TruSeq_and_nextera_adaptors.fa --trim-n -m 50 -o ${output}_R1.trimmed.fastq.gz -p ${output}_R2.trimmed.fastq.gz ${path}/${R1} ${path}/${R2}
    '''
    input:
        f1=join(WORKDIR,"fastqs",'{sample_id}.R1.fastq.gz'),
        f2=join(WORKDIR,"fastqs",'{sample_id}.R2.fastq.gz')
    envmodules:
        TOOLS["cutadapt"]["version"]
    threads: getthreads("trim_adaptors")
    params:
        rname="trim_adaptors",
        adaptors=config["adaptors"]
    output:
        f1=temp(join(RESULTSDIR,'01_trim','{sample_id}.R1.trimmed.fastq.gz')),
        f2=temp(join(RESULTSDIR,'01_trim','{sample_id}.R2.trimmed.fastq.gz'))
    shell:
        """
        cutadapt -j 32 -b file:{params.adaptors} -B file:{params.adaptors} --trim-n -m 50 -o {output.f1} -p {output.f2} {input.f1} {input.f2}
        """

rule assembly:
    '''
    Assemble read pairs 
    #pear -p 0.0001 -f ${output}_R1.trimmed.fastq.gz -r ${output}_R2.trimmed.fastq.gz -o ${output} -j 32
    '''
    input:
        f1=rules.trim_adaptors.output.f1,
        f2=rules.trim_adaptors.output.f2
    envmodules:
        TOOLS["pear"]["version"]
    threads: getthreads("assembly")
    params:
        rname="assembly",
        localtmp=join(RESULTSDIR,'tmp','assembly'),
        sp='{sample_id}'
    output:
        merged_fq=temp(join(RESULTSDIR,'02_assembled','{sample_id}.{species}.assembled.fastq.gz'))
    shell:
        """
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        # run assembly
        pear -p 0.0001 -f {input.f1} -r {input.f2} -o $tmp_dir/{params.sp} -j 32
        
        # gzip and move
        gzip -f $tmp_dir/{params.sp}.assembled.fastq
        mv $tmp_dir/{params.sp}.assembled.fastq.gz {output.merged_fq}
        """

rule alignment:
    '''
    Map and discard unnmaped
    #bowtie2 -p 32 -x hg19 -U ${output}.assembled.fastq.gz -S ${output}.assembled.hg19.bam
    #samtools view -b -F 260 ${output}.assembled.hg19.bam > ${output}.mapped.hg19.bam
    #bedtools bamtobed -i ${output}.mapped.hg19.bam > ${output}.mapped.hg19.bed
    '''
    input:
        assembled=rules.assembly.output.merged_fq
    envmodules:
        TOOLS["bowtie2"]["version"],
        TOOLS["samtools"]["version"],
        TOOLS["bedtools"]["version"],
    threads: getthreads("alignment")
    params:
        rname="alignment",
        index_dir=INDEX_LOCATION,
    output:
        bam=join(RESULTSDIR,'03_aligned','01_bam','{sample_id}.{species}.assembled.bam'),
        mapped_bam=join(RESULTSDIR,'03_aligned','01_bam','{sample_id}.{species}.mapped.bam'),
        bed=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.{species}.mapped.bed')
    shell:
        """
        bowtie2 -p 32 -x {params.index_dir} -U {input.assembled} -S {output.bam}
        samtools view -b -F 260 {output.bam} > {output.mapped_bam}
        bedtools bamtobed -i {output.mapped_bam} > {output.bed}
        """

rule all_hist_frags:
    '''
    Make histogram of fragment lengths
     #Rscript hist.r ${output}.mapped.hg19.bed ${output}.mapped.hg19.length_hist.csv
    '''
    input:
        bed=rules.alignment.output.bed,
    envmodules:
        TOOLS["R"]["version"]
    threads: getthreads("all_hist_frags")
    params:
        rname="all_hist_frags",
        rscript=join(WORKDIR,"scripts","hist.r"),
        localtmp=join(RESULTSDIR,'tmp','histo_frags'),
    output:
        hist=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.{species}.length_hist.all.csv'),
        png=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.{species}.length_hist.all.png')
    shell:
        """
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        # subset bed file to include only columns needed
        awk -v OFS='\t' '{{print $1,$2,$3}}' {input.bed} > $tmp_dir/tmp.bed

        # calculate length
        awk '{{$4 = ($3-$2); print}}' $tmp_dir/tmp.bed > $tmp_dir/tmp1.bed

        Rscript {params.rscript} $tmp_dir/tmp1.bed {output.hist} {output.png}
        """