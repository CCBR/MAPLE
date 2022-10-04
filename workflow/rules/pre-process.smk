rule create_bed_file:
    '''
    If a bed file has not already been created for target gene_list, then create it
    If it has been created then copy it to the working dir
    '''
    input:
        master=config["master_bed_file"],
    threads: getthreads("create_bed_file")
    params:
        rname="create_bed_file",
        create_flag=config["run_create_bed"],
        select_flag=config["run_select_bed"],
        gene_list=config["genes_of_interest"],
        pi_created_selected_bed=config["pi_created_selected_bed"]
    output:
        selected_bed=join(RESULTSDIR,'00_selected_bed','selected_bed.bed'),
    shell:
        """
        if [[ $create_flag == "N" ]] & [[ $select_flag == "Y" ]]; then
            grep -Fwf {params.gene_list} {input.master} | awk -v OFS='\t' '{{print \$1,\$2,\$3}}'> {output.selected_bed}
        else
            cp {params.pi_created_selected_bed} {output.selected_bed}
        fi
        """

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
        f1=join(RESULTSDIR,'01_trim','{sample_id}.R1.trimmed.fastq.gz'),
        f2=join(RESULTSDIR,'01_trim','{sample_id}.R2.trimmed.fastq.gz')
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
        merged_fq=join(RESULTSDIR,'02_assembled','{sample_id}.assembled.fastq.gz')
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
        species=config["species"],
        index_dir=INDEXDIR,
    output:
        bam=join(RESULTSDIR,'03_aligned','01_bam','{sample_id}.assembled.bam'),
        mapped_bam=join(RESULTSDIR,'03_aligned','01_bam','{sample_id}.mapped.bam'),
        bed=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.mapped.bed')
    shell:
        """
        bowtie2 -p 32 -x {params.index_dir}/{params.species} -U {input.assembled} -S {output.bam}
        samtools view -b -F 260 {output.bam} > {output.mapped_bam}
        bedtools bamtobed -i {output.mapped_bam} > {output.bed}
        """

rule hist_frags:
    '''
    Make histogram of fragment lengths
     #Rscript hist.r ${output}.mapped.hg19.bed ${output}.mapped.hg19.length_hist.csv
    '''
    input:
        bed=rules.alignment.output.bed,
    envmodules:
        TOOLS["R"]["version"]
    threads: getthreads("hist_frags")
    params:
        rname="hist_frags_all",
        rscript=join(WORKDIR,"scripts","hist.r")
    output:
        hist=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.length_hist_all.csv'),
        png=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.length_hist_all.png')
    shell:
        """
        Rscript {params.rscript} {input.bed} {output.hist} {output.png}
        """

rule select_bed:
    '''
    create select bed, if required
    check access for intervals bed 

    # bedtools -a ${output}.mapped.hg19.bed -b intervals_of_interests.bed > ${output}.selected.hg19.bed
    # awk '{ if ($3-$2 >= 140 && $3-$2 <= 160) print $0}' ${output}.selected.hg19.bed > ${output}.selected.hg19.140-160.bed

    '''
    input:
        bed=rules.alignment.output.bed,
        interval_bed=rules.create_bed_file.output.selected_bed
    envmodules:
        TOOLS["R"]["version"]
    threads: getthreads("hist_frags")
    params:
        rname="hist_frags",
        rscript=join(WORKDIR,"scripts","hist.r"),
        localtmp=join(RESULTSDIR,'tmp','selectbed'),
        f_min=config["fragment_length_min"],
        f_max=config["fragment_length_max"],
    output:
        selected_bed=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.mapped.selected.bed')
    shell:
        """
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        bedtools -a {input.bed} -b {input.interval_bed} > $tmp_dir/tmp.bed
        awk '{{ if (\$3-\$2 >= {params.f_min} && \$3-\$2 <= {params.f_max}) print \$0}}' $tmp_dir/tmp.bed > {output.selected_bed}
        """

rule selected_hist_frags:
    '''
    Make histogram of selected fragment lengths
     #Rscript hist.r ${output}.mapped.hg19.bed ${output}.mapped.hg19.length_hist.csv
    '''
    input:
        bed=rules.select_bed.output.selected_bed,
    envmodules:
        TOOLS["R"]["version"]
    threads: getthreads("hist_frags")
    params:
        rname="hist_frags_select",
        rscript=join(WORKDIR,"scripts","hist.r")
    output:
        hist=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.length_hist_selected.csv'),
        png=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.length_hist_selected.png')
    shell:
        """
        Rscript {params.rscript} {input.bed} {output.hist} {output.png}
        """