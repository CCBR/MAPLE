# set species
species=config["species"]

# set reference source
reference_source=config["reference_source"]

rule build_index_files:
    '''
    If the index has not been created, then create it
    '''
    threads: 
        getthreads("create_index_file")
    envmodules:
        TOOLS["bowtie2"]["version"],
    params:
        rname="create_index_file",
        index_dir=INDEX_LOCATION,
        index_shorthand=INDEX[species][reference_source]["shorthand"],
        source_link=INDEX[species][reference_source]["source"],
        sp = config["species"],
        ref_type = config["reference_source"]
    output:
        built_index=GENOME_INDEX_FILE
    shell:
        """
        # download the zip file to the dir
        wget -P {params.index_dir} {params.source_link}
        
        # unzip the dir
        cd {params.index_dir}
        gunzip {params.index_dir}/${params.index_shorthand}.zip
        
        # run the sh file
        cd {params.index_dir}/{params.index_shorthand}
        sh make_${params.sp}.sh

        # remove fastas
        rm *fa
        """