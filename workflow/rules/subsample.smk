def get_selected_bed(wildcards):
    # read in bed_list
    bed_list=join(WORKDIR,"resources",config["bed_list_name"])
    bed_df = pd.read_csv(bed_list)
    
    # subset for selected_shorthand
    df_sub = bed_df[(bed_df['selected_shorthand']==wildcards.selected_shorthand)]

    # return full path of bedfile
    return(df_sub.iloc[0]['selected_bed'])

rule select_bed:
    '''
    Create a sub-selected bed file based on user input

    # bedtools -a ${output}.mapped.hg19.bed -b intervals_of_interests.bed > ${output}.selected.hg19.bed
    # awk '{ if ($3-$2 >= 140 && $3-$2 <= 160) print $0}' ${output}.selected.hg19.bed > ${output}.selected.hg19.140-160.bed        
    '''
    input:
        bed=rules.alignment.output.bed,
        interval_bed=get_selected_bed
    envmodules:
        TOOLS["R"]["version"],
        TOOLS["bedtools"]["version"]
    threads: getthreads("hist_frags")
    params:
        rname="select_bed",
        rscript=join(WORKDIR,"scripts","hist.r"),
        localtmp=join(RESULTSDIR,'tmp','selectbed'),
        f_min=config["fragment_length_min"],
        f_max=config["fragment_length_max"]
    output:
        selected_bed=temp(join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.{species}.{min_length}-{max_length}.{selected_shorthand}.bed'))
    shell:
        """
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        bedtools intersect -a {input.bed} -b {input.interval_bed} > $tmp_dir/tmp.bed
        awk '{{ if ($3-$2 >= {params.f_min} && $3-$2 <= {params.f_max}) print $0}}' $tmp_dir/tmp.bed > {output.selected_bed}
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
        rscript=join(WORKDIR,"scripts","hist.r"),
        localtmp=join(RESULTSDIR,'tmp','select_frags'),
    output:
        hist=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.{species}.length_hist.{min_length}-{max_length}.{selected_shorthand}.csv'),
        png=join(RESULTSDIR,'03_aligned','03_histograms','{sample_id}.{species}.length_hist.{min_length}-{max_length}.{selected_shorthand}.png')
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