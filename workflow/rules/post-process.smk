def get_dyad_input(wildcards):
    if (config["run_select_bed"]=="Y"):
        dyad_input=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.mapped.selected.bed')
    else:
        dyad_input=join(RESULTSDIR,'03_aligned','02_bed','{sample_id}.mapped.bed')
    return(dyad_input)

rule dyad_analysis:
    '''
    Find fragment centers (DYADs) and make histogram (Occurrences)

    python3 WeigthedDYADposition.py ${output}.selected.hg19.140-160.bed ${output}.DYADs
    sort -k1,1 -k2n,2 ${output}.DYADs > ${output}.DYADs.sorted
    python Uniq_Position.py ${output}.DYADs.sorted ${output}.DYADs.hist
    '''
    input:
        bed=get_dyad_input
    envmodules:
        TOOLS["python37"]["version"]
    threads: getthreads("find_dyads")
    params:
        rname="find_dyads",
        position_script=join(WORKDIR,"scripts","WeigthedDYADposition.py"),
        hist_script=join(WORKDIR,"scripts","Uniq_Position.py"),
        csv_script=join(WORKDIR,"scripts","ALU_DAC.py"),
        limit=config["limit"],
        max_d=config["max_distance"]
    output:
        dyads=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.DYADs'),
        sorted=join(RESULTSDIR,'04_dyads','01_DYADs','{sample_id}.sorted.DYADs'),
        hist=join(RESULTSDIR,'04_dyads','02_histograms','{sample_id}.DYADs.hist'),
        csv=join(RESULTSDIR,'04_dyads','03_CSV','{sample_id}.DAC.csv'),
    shell:
        """
        python3 {params.position_script} {input.bed} {output.dyads}
        sort -k1,1 -k2n,2 {output.dyads} > {output.sorted}
        python {params.hist_script} {output.sorted} {output.hist}
        python {params.csv_script} {output.hist} {params.limit} {params.max_d} {output.csv}
        """