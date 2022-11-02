def get_selected_bed(wildcards):
    # read in bed_list
    bed_list=join(WORKDIR,"resources","bed_lists.csv")
    bed_df = pd.read_csv(bed_list)
    
    # subset for selected_shorthand
    df_sub = bed_df[(bed_df['selected_shorthand']==wildcards.selected_shorthand)]

    # return full path of bedfile
    return(df_sub.iloc[0]['selected_bed'])

rule create_bed_file:
    '''    
    Copy the bed file to the working dir
    '''
    input:
        master=config["master_bed_file"],
        bed_file=get_selected_bed
    threads: 
        getthreads("create_bed_file")
    params:
        rname="create_bed_file",
    output:
        selected_bed=temp(join(RESULTSDIR,'00_selected_bed','{selected_shorthand}.bed'))
    shell:
        """
        cp {input.bed_file} {output.selected_bed}
        """