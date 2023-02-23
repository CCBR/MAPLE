localrules: create_master_gene_lists, create_master_table

def increase_time(wildcards, attempt):
    if attempt == 1: 
        clean_time="01-00:00:00"
    else:
        clean_time="0"+attempt+"-00:00:00"
    return clean_time

rule create_master_gene_lists:
    """
    Processing individual genes for master list is time consuming. This step will break the gene list into useable chunks for use 
    the rule create_indiv_master_table
    """
    input:
        master_bed=config["bed_list_name"],
    threads: getthreads("create_master_gene_lists")
    params:
        rname="create_master_gene_lists",
        gene_list_n=NUMBEROFGENESLISTS,
        prefix=join(RESULTSDIR,'04_dyads','04_master_table',"{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.gene_list_")
    output:
        master_gene_list=join(RESULTSDIR,'04_dyads','04_master_table',"{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.master_gene_list.csv"),
    shell:
        """
        #pull gene names from input gene list
        awk \'{{print $4}}\' {input.master_bed} | sort | uniq > {output.master_gene_list}

        # determine number of genes
        total_num_genes=`cat {output.master_gene_list} | wc -l`

        # define number of gene lists
        n={params.gene_list_n}

        # determine number of genes in each list; round
        div=$((total_num_genes/n))
        number_of_genes=`echo $div | awk \'{{print int($1+1.5)}}\'`
        echo "number of genes included in each split: $number_of_genes"

        # split the file
        split -l $number_of_genes --numeric-suffixes=1 --additional-suffix .txt {output.master_gene_list} {params.prefix}

        # fix leading zero
        for i in {{1..9}}; do
            old_filename="{params.prefix}0$i.txt"
            new_filename="{params.prefix}$i.txt"
            if [[ -f $old_filename ]]; then
                mv $old_filename $new_filename
            fi
        done
        """

rule create_indiv_master_table:
    '''
    Create list of genes based on input bed file
    For each gene calculate the DYAD and generate corrected, master_table output
    NROW=max_distance and NCOL=number of genes in bed file
    '''
    input:
        master_bed=config["bed_list_name"],
        sample_bed=rules.alignment.output.bed,
        gene_lists=rules.create_master_gene_lists.output.master_gene_list,
    envmodules:
        TOOLS["bedtools"]["version"],
        TOOLS["python37"]["version"]
    threads: getthreads("create_indiv_master_table")
    resources:
        time=increase_time
    params:
        rname="create_indiv_master_table",
        localtmp=join(RESULTSDIR,'tmp','merged'),
        position_script=join(WORKDIR,"scripts","WeigthedDYADposition.py"),
        hist_script=join(WORKDIR,"scripts","uniq_position.py"),
        dac_script=join(WORKDIR,"scripts","DAC.py"),
        max_d=config["max_distance"],
        dac_corrected_script=join(WORKDIR,"scripts","DAC_denominator.py"),
        split_file=join(RESULTSDIR,'04_dyads','04_master_table',"{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.gene_list_{n}.txt")
    output:
        n_master_table=temp(join(RESULTSDIR,'04_dyads','04_master_table','{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.split_table_{n}.DAC.corrected.csv')),
    shell:
        """
        # create tmp dir to hold data
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        # for each gene calculate DYAD
        while IFS=\',\' read -a gene 
            do
                # make dir for gene tmp files
                echo "$gene"
                mkdir $tmp_dir/$gene
                
                # subset master bed file for gene only info
                gene_bed="$tmp_dir/$gene/$gene.bed"
                cat {input.master_bed} | grep -w $gene > $gene_bed
                
                # create bed files of only gene info
                echo "--bed"
                mapped_gene_bed="$tmp_dir/$gene/$gene.mapped.bed"
                bedtools intersect -a {input.sample_bed} -b $gene_bed > $mapped_gene_bed

                ################################################################################
                # calculate DYADS
                ################################################################################
                echo "--dyads"
                dyads="$tmp_dir/$gene/DYADs"
                python3 {params.position_script} $mapped_gene_bed $dyads

                # sort
                echo "--sorting"
                sorted="$tmp_dir/$gene/sorted.DYADs"
                sort -k1,1 -k2n,2 $dyads > $sorted

                #create histogram
                echo "--histo"
                hist="$tmp_dir/$gene/DYADs.hist"
                python {params.hist_script} $sorted $hist

                # set limit
                limit=`echo {output.n_master_table} | awk -F 'lim' '{{ print $2 }}' | cut -f1 -d"."`

                # compute auto-correlation
                echo "--computing"
                csv="$tmp_dir/$gene/DAC.csv"
                python {params.dac_script} $hist $limit {params.max_d} $csv
                average_length=$(awk '{{ SUM += ($3-$2); n++}} END {{print(int(SUM/n))}}' $gene_bed)
                echo "The average length is: $average_length"

                echo "--correcting"
                corrected_csv="$tmp_dir/$gene/DAC.corrected.csv"
                python {params.dac_corrected_script} $hist $limit {params.max_d} $average_length $corrected_csv

                # if the mater file doesn't exist, create it with the row names, otherwise, only add gene column
                echo "--cleaning"
                merged_table="$tmp_dir/merged_master_table.csv"
                tmp_table="$tmp_dir/tmp_master_table.csv"
                tmp_col="$tmp_dir/$gene/col.csv"
                tmp_col2="$tmp_dir/$gene/col2.csv"
        
                # round to 5 decimals
                # remove the header
                # pull the calculation col
                # multiple this col by 10
                # round to 4 decimals
                awk -F"," \'{{$2=sprintf("%.5f",$2)}}1\' $corrected_csv | grep -v "Dist" | awk \'{{print $2}}\' | awk \'{{$1=$1*10; print $0}}\' | awk -F"," \'{{$1=sprintf("%.4f",$1)}}1\' > $tmp_col

                # save files for genes
                if [[ $gene == "C3" ]] || [[ $gene == "TTC3" ]] || [[ $gene == "TTC34" ]]; then
                    cp $tmp_col /data/sevillas2/tmp/master_table/$gene.csv
                fi

                # if the master file doesn't exist, create it with the row names, otherwise, only add gene column
                if [[ ! -f $merged_table ]]; then
                    echo "$gene" > $merged_table
                    cat $tmp_col >> $merged_table
                else
                    # set the tmp table to be merged
                    cp $merged_table $tmp_table

                    # add header then col file
                    echo "$gene" > $tmp_col2
                    cat $tmp_col >> $tmp_col2

                    # add calc column to output
                    paste -d" " $tmp_table $tmp_col2 > $merged_table
                fi

                #cleanup
                rm -rf $tmp_dir/$gene
        done <  {params.split_file}

        # move final file
        cp $tmp_dir/merged_master_table.csv {output.n_master_table}
        """

rule create_master_table:
    """
    Processing indiv_master_tables into one master table
    """
    input:
        indiv_tables=expand(join(RESULTSDIR,'04_dyads','04_master_table','{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.split_table_{n}.DAC.corrected.csv'),sample_id=SAMPLES, species=SPECIES,min_length=MIN_LENGTH, max_length=MAX_LENGTH, limit=LIMITSIZE,n=RANGEOFGENELISTS)
    threads: getthreads("create_master_table")
    params:
        rname="create_master_table",
        localtmp=join(RESULTSDIR,'tmp','merged'),
        outDIR=join(RESULTSDIR,'04_dyads','04_master_table'),
        gene_lists=join(RESULTSDIR,'04_dyads','04_master_table',"{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.gene_list_"),
        sort_script=join(WORKDIR,"scripts","sort.py"),
    output:
        master_table=join(RESULTSDIR,'04_dyads','04_master_table','{sample_id}.{species}.{min_length}-{max_length}.lim{limit}.master_table.DAC.corrected.csv'),
    shell:
        """
        # create tmp dir to hold data
        tmp_dir="/lscratch/${{SLURM_JOB_ID}}"
        if [[ ! -d $tmp_dir ]]; then
            tmp_dir={params.localtmp}
            if [[ -d $tmp_dir ]]; then rm -r $tmp_dir; fi 
            mkdir -p $tmp_dir
        fi

        # split merged files into array
        master_table_list=({input.indiv_tables})

        # for each of the merged files, merge into one final file
        tmp_master=$tmp_dir/tmp_master.csv
        tmp_col=$tmp_dir/tmp_col.csv
        merged_master=$tmp_dir/merged_master.csv
        for f in ${{master_table_list[@]}}; do
            if [[ ! -f $merged_master ]]; then
                    cp $f $merged_master
                else
                    # set the tmp table to be merged
                    cp $merged_master $tmp_master 
                    
                    # paste remaining columns together
                    paste -d"," $tmp_master $f > $merged_master
            fi
        done

        # reformat
        sed -i -e 's/\s\+/,/g' $merged_master 
        sed -i -e 's/\t/,/g' $merged_master

        # create final merged file
        ## create header
        echo "Dist" > $tmp_dir/tmp_header
        head -n1 $merged_master > $tmp_dir/tmp_header2
        paste -d"," $tmp_dir/tmp_header $tmp_dir/tmp_header2 > $tmp_dir/final_header
        
        ## create distances only file
        cat $merged_master | grep -v "A" > $tmp_dir/final_dists
        
        ## create rownum col
        awk \'{{print NR}}\' $tmp_dir/final_dists > $tmp_dir/final_rownums

        # create final file
        ## paste rownums to distances
        paste -d"," $tmp_dir/final_rownums $tmp_dir/final_dists > $tmp_dir/final_rowsnums_dists
        ## add row names
        cat $tmp_dir/final_header > $tmp_dir/cleaned_output
        ## add distances
        cat $tmp_dir/final_rowsnums_dists >> $tmp_dir/cleaned_output

        # transpose the file, cleanup
        awk -F"," \'
        {{ 
            for (i=1; i<=NF; i++)  {{
                a[NR,i] = $i
            }}
        }}
        NF>p {{ p = NF }}
        END {{    
            for(j=1; j<=p; j++) {{
                str=a[1,j]
                for(i=2; i<=NR; i++){{
                    str=str" "a[i,j];
                }}
                print str
            }}
        }}\' $tmp_dir/cleaned_output > $tmp_dir/transposed_output.csv
        sed -i -e 's/\s\+/,/g' $tmp_dir/transposed_output.csv
        sed -i -e 's/\t/,/g' $tmp_dir/transposed_output.csv

        # sort the file via python
        ## NOTE: sort command in bash does not sort the same way that the reference master_table was created
        ## using this script to sort them the same way for ease of use by PI
        python {params.sort_script} $tmp_dir/transposed_output.csv {output.master_table}

        # cleanup gene lists
        if [[ -f ${params.gene_lists}* ]]; then
            rm ${params.gene_lists}*
        fi
        """