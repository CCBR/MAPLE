# set output location
output_contrast_location=config["output_contrast_location"]

rule calculate_DYADs:
    '''
    Find fragment centers (DYADs) and make histogram (Occurrences)

    python3 WeigthedDYADposition.py ${output}.selected.hg19.140-160.bed ${output}.DYADs
    sort -k1,1 -k2n,2 ${output}.DYADs > ${output}.DYADs.sorted
    python Uniq_Position.py ${output}.DYADs.sorted ${output}.DYADs.hist
    '''
    input:
        bed=rules.select_bed.output.selected_bed,
    envmodules:
        TOOLS["python37"]["version"]
    threads: getthreads("calculate_DYADs")
    params:
        rname="calculate_DYADs",
        position_script=join(WORKDIR,"scripts","WeigthedDYADposition.py"),
        hist_script=join(WORKDIR,"scripts","uniq_position.py"),
    output:
        dyads=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.{species}.{min_length}-{max_length}.{selected_shorthand}.DYADs'),
        sorted=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.{species}.{min_length}-{max_length}.{selected_shorthand}.sorted.DYADs'),
        hist=join(RESULTSDIR,'04_dyads','02_histograms','{sample_id}.{species}.{min_length}-{max_length}.{selected_shorthand}.DYADs.hist'),
    shell:
        """
        # calculate DYADS
        python3 {params.position_script} {input.bed} {output.dyads}
        
        # sort
        sort -k1,1 -k2n,2 {output.dyads} > {output.sorted}
        
        #create histogram
        python {params.hist_script} {output.sorted} {output.hist}
        """

rule dyad_analysis:
    '''
    ## Compute auto-correlation ###########
    python ALU_DAC.py ${output}.DYADs.hist ${limit} ${max_dist} ${output}.DAC.csv
    '''
    input:
        hist=rules.calculate_DYADs.output.hist,
        source_bed=rules.create_bed_file.output.selected_bed
    envmodules:
        TOOLS["python37"]["version"]
    threads: getthreads("dyad_analysis")
    params:
        rname="dyad_analysis",
        dac_script=join(WORKDIR,"scripts","DAC.py"),
        dac_corrected_script=join(WORKDIR,"scripts","DAC_denominator.py"),
        max_d=config["max_distance"]
    output:
        csv=join(RESULTSDIR,'04_dyads','03_CSV','{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.{selected_shorthand}.DAC.csv'),
        corrected_csv=join(RESULTSDIR,'04_dyads','03_CSV','{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.{selected_shorthand}.DAC.corrected.csv'),
    shell:
        """
        # set limit
        limit=`echo {output.csv} | awk -F 'lim' '{{ print $2 }}' | cut -f1 -d"."`

        # compute auto-correlation
        python {params.dac_script} {input.hist} $limit {params.max_d} {output.csv}
        
        average_length=$(awk '{{ SUM += ($3-$2); n++}} END {{print(int(SUM/n))}}' {input.source_bed})
        echo "The average length is: $average_length"
        
        # correct based on length of genome
        python {params.dac_corrected_script} {input.hist} $limit {params.max_d} $average_length {output.corrected_csv}
        """

rule merge_DACs:
    '''
    merge all DAC files into one file
    '''
    input:
        dacs=CONTRAST_FILES
    envmodules:
        TOOLS["R"]["version"]
    threads: getthreads("merge_DACs")
    params:
        rname="merge_DACs",
        localtmp=join(RESULTSDIR,'tmp','merged'),
    output:
        merged=join(output_contrast_location,'final_' + CONTRASTS_CLEAN_LIST + '.{min_length}-{max_length}.lim{limit}.{selected_shorthand}.DAC.csv')
    shell:
        """
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        file_list=({input.dacs})
        counter=0
        header=""

        # for each output DAC file, join 
        for f in ${{file_list[@]}}; do
            if [[ $counter -eq 0 ]]; then
                echo "working on $f" 
                cat $f > $tmp_dir/join.tmp
            elif [[ $counter -gt 0 ]]; then
                echo "now for $f"
                join -a1 -a2 -e 1  -t $',' -o auto $tmp_dir/join.tmp ${{file_list[$counter]}} > $tmp_dir/join.tmp.1
                mv $tmp_dir/join.tmp.1 $tmp_dir/join.tmp
            fi

            # remove file from file_list and continue; increase counter
            file_list=("${{file_list[@]/${{file_list[$counter]}}}}")
            counter=$((counter+1))

            # clean file name to only include sample name
            clean_f=`echo ${{f##*/}} | cut -f1 -d"."`
            header="$header,$clean_f"
        done

        # echo the header and join file to final file
        echo "Dist$header" > $tmp_dir/header.tmp
        cat $tmp_dir/header.tmp $tmp_dir/join.tmp > $tmp_dir/join.tmp.1
        cat $tmp_dir/join.tmp.1 | grep -v "Dist,DAC" > {output.merged}
        """