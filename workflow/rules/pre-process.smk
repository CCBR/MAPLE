rule trim_adapters:
    '''
    Read samples.tsv to determine fastq files, trim adapters
    '''
    input:
        f1=join(fq_dir,'{sample_id}_R1.fastq.gz'),
        f2=join(fq_dir,'{sample_id}_R2.fastq.gz')
    envmodules:
        TOOLS["cutadapt"]["version"]
    container: config["masterdocker"]    
    threads: getthreads("trim_adapters")
    params:
        rname="trim_adapters",
        adapters=config[adapters]
    output:
        f1=join(RESULTSDIR,'trim','{sample_id}_R1.trimmed.fastq.gz'),
        f2=join(RESULTSDIR,'trim','{sample_id}_R1.trimmed.fastq.gz')
    shell:
    """
        #cutadapt -j 32 -b file:TruSeq_and_nextera_adapters.fa -B file:TruSeq_and_nextera_adapters.fa --trim-n -m 50 -o ${output}_R1.trimmed.fastq.gz -p ${output}_R2.trimmed.fastq.gz ${path}/${R1} ${path}/${R2}
        cutadapt -j 32 -b file:{params.adapters} -B file:{params.adapters} --trim-n -m 50 -o {output.f1} -p {output.f2} {input.f1} {input.f2}
    """

rule assembly:
    '''
    Assemble read pairs 
    '''
    input:
        f1=rules.trim_adapters.output.f1,
        f2=rules.trim_adapters.output.f2
    envmodules:
        TOOLS["pear"]["version"]
    container: config["masterdocker"]    
    threads: getthreads("assembly")
    params:
        rname="assembly"
    output:
        merged_fq=join(RESULTSDIR,'assembled','{sample_id}.assembled.fastq.gz')
    shell:
    """
        tmp_dir="/lscratch/${{SLURMJOB}}"
        export $tmp_dir

        #pear -p 0.0001 -f ${output}_R1.trimmed.fastq.gz -r ${output}_R2.trimmed.fastq.gz -o ${output} -j 32
        pear -p 0.0001 -f {input.f1} -r {input.f2} -o $tmp_dir -j 32
        
        # gzip and move
        gzip -f $tmp_dir/*.assembled.fastq
        mv $tmp_dir/*.assembled.fastq.gz {output.merged_fq}
    """


rule alignment:
    '''
    Map and discard unnmaped
    '''
    input:
        assembled=rules.assembly.output.merged_fq
    envmodules:
        TOOLS["bowtie2"]["version"]
        TOOLS["samtools"]["version"]
        TOOLS["bedtools"]["version"]
    container: config["masterdocker"]    
    threads: getthreads("alignment")
    params:
        rname="alignment",
        species=config['species']
    output:
        bam=join(RESULTSDIR,'aligned','{sample_id}.assembled.bam')
        mapped_bam=join(RESULTSDIR,'aligned','{sample_id}.mapped.bam')
        bed=join(RESULTSDIR,'aligned','{sample_id}.mapped.bed')
    shell:
    """
        #bowtie2 -p 32 -x hg19 -U ${output}.assembled.fastq.gz -S ${output}.assembled.hg19.bam
        #samtools view -b -F 260 ${output}.assembled.hg19.bam > ${output}.mapped.hg19.bam
        #bedtools bamtobed -i ${output}.mapped.hg19.bam > ${output}.mapped.hg19.bed

        bowtie2 -p 32 -x {params.species} -U {input.assembled} -S {output.bam}
        samtools view -b -F 260 {output.bam} > {output.mapped_bam}
        bedtools bamtobed -i {output.mapped_bam} > {output.bed}
    """

rule select_bed:
    '''
    create select bed - FLAG FOR IF
    check access for intervals bed 
    '''
    input:
        bed=rules.alignment.output.bed
        interval_bed=config["intervals_of_interest"]
    envmodules:
        TOOLS["R"]["version"]
    container: config["masterdocker"]    
    threads: getthreads("hist_frags")
    params:
        rname="hist_frags"
        rscript=join(SCRIPTS,"hist.r")
    output:
        bed=join(RESULTSDIR,'aligned','{sample_id}.mapped.140_to_160.bed')
    shell:
    """
        tmp_dir="/lscratch/${{SLURMJOB}}"
        export $tmp_dir

        # bedtools -a ${output}.mapped.hg19.bed -b intervals_of_interests.bed > ${output}.selected.hg19.bed
        # awk '{ if ($3-$2 >= 140 && $3-$2 <= 160) print $0}' ${output}.selected.hg19.bed > ${output}.selected.hg19.140-160.bed

        bedtools -a {input.bed} -b {input.interval_bed} > $tmp_dir/tmp.bed
        awk '{ if ($3-$2 >= 140 && $3-$2 <= 160) print $0}' $tmp_dir/tmp.bed > {output.selected_bed}
    """

# if a selected bed file is created, then use this file for histograms, otherwise use the 'everything' bed file
def get_bed_file():
    if select_bed is chosen:
        bed_file=join(RESUTLSDIR,'aligned',"{wildcards.sample_id}.mapped.140_to_160.bed")
    else
        bed_file=join(RESUTLSDIR,'aligned',"{wildcards.sample_id}.mapped.bed")
    fi
    return(bed_file)


rule hist_frags:
    '''
    Make histogram of fragment lengths
    '''
    input:
        bed=get_bed_file()
    envmodules:
        TOOLS["R"]["version"]
    container: config["masterdocker"]    
    threads: getthreads("hist_frags")
    params:
        rname="hist_frags"
        rscript=join(SCRIPTS,"hist.r")
    output:
        hist=join(RESULTSDIR,'{sample_id}.length_hist.csv')
    shell:
    """
        #Rscript hist.r ${output}.mapped.hg19.bed ${output}.mapped.hg19.length_hist.csv
        Rscript {params.rscript} {input.bed} {output.hist}
    """
