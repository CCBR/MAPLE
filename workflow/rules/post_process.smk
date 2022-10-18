# set selection shorthand
selection_shorthand=config["selection_shorthand"]

# set min / max lengths of fragment
min_length=config["fragment_length_min"]
max_length=config["fragment_length_max"]


 # set ouput location
output_contrast_location=config["output_contrast_location"]

def get_dyad_input(wildcards):
    if (config["run_select_bed"]=="Y"):
        dyad_input=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.{species}.' + min_length + "-" + max_length + "." + selection_shorthand +".bed")
    else:
        dyad_input=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.{species}.mapped.bed')
    return(dyad_input)

print("OK")
print(join(output_contrast_location,'final_' + CONTRASTS_CLEAN_LIST + "."))

rule dyad_analysis:
    '''
    Find fragment centers (DYADs) and make histogram (Occurrences)

    python3 WeigthedDYADposition.py ${output}.selected.hg19.140-160.bed ${output}.DYADs
    sort -k1,1 -k2n,2 ${output}.DYADs > ${output}.DYADs.sorted
    python Uniq_Position.py ${output}.DYADs.sorted ${output}.DYADs.hist

    ## Compute auto-correlation ###########
    python ALU_DAC.py ${output}.DYADs.hist ${limit} ${max_dist} ${output}.DAC.csv
    '''
    input:
        bed=get_dyad_input
    envmodules:
        TOOLS["python37"]["version"]
    threads: getthreads("find_dyads")
    params:
        rname="dyad_analysis",
        position_script=join(WORKDIR,"scripts","WeigthedDYADposition.py"),
        hist_script=join(WORKDIR,"scripts","uniq_position.py"),
        csv_script=join(WORKDIR,"scripts","DAC.py"),
        limit=config["limit"],
        max_d=config["max_distance"]
    output:
        dyads=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.{species}.{min_length}-{max_length}.{selection_shorthand}.DYADs'),
        sorted=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.{species}.{min_length}-{max_length}.{selection_shorthand}.sorted.DYADs'),
        hist=join(RESULTSDIR,'04_dyads','02_histograms','{sample_id}.{species}.{min_length}-{max_length}.{selection_shorthand}.DYADs.hist'),
        csv=join(RESULTSDIR,'04_dyads','03_CSV','{sample_id}.{species}.{min_length}-{max_length}.{selection_shorthand}.DAC.csv'),
    shell:
        """
        python3 {params.position_script} {input.bed} {output.dyads}
        sort -k1,1 -k2n,2 {output.dyads} > {output.sorted}
        python {params.hist_script} {output.sorted} {output.hist}
        python {params.csv_script} {output.hist} {params.limit} {params.max_d} {output.csv}
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
        merged=join(output_contrast_location,'final_' + CONTRASTS_CLEAN_LIST + '.{min_length}-{max_length}.{selection_shorthand}.DAC.csv')
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
            echo "working $f"
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