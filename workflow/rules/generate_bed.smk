selected_shorthand=config["selected_shorthand"]

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
        gene_list=config["genes_of_interest"],
        selected_bed=config["selected_bed"],
        selected_shorthand=config["selected_shorthand"]
    output:
        selected_bed=join(RESULTSDIR,'00_selected_bed',selected_shorthand + '.bed'),
    shell:
        """
        # if the selected bed file does not exist, then create it
        if [[ ! -f {params.selected_bed} ]]; then
            grep -Fwf {params.gene_list} {input.master} | awk -v OFS='\t' '{{print $1,$2,$3}}'> {output.selected_bed}
        else
            cp {params.selected_bed} {output.selected_bed}
        fi
        """
