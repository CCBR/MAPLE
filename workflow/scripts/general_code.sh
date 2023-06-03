flag_stage=$1
# create_bed
# pipe_prep
# pipe_run
# mv_samples
# create_links

# create gene list
if [[ $flag_stage == "create_genelist" ]]; then
    # file=$(cat 1000.Genes.ALU-rich.bed)
    for line in $file; do
        start=`echo $line | awk -F" " '{print $2}'`
        echo $start
    done

    awk -F"	" '{print $1,$2,$3}' 1000.Genes.ALU-rich.bed

    #reverse engineer gene list
    awk 'FNR==NR{arr[$1,$2];next} (($1,$2) in arr)' 1000.Genes.ALU-rich.bed /data/CCBR_Pipeliner/MAPLE/bed_files/hg19_protein-coding_genes.bed > 1000.Genes.ALU-rich.txt

    #create bed file
    grep -Fwf 1000.Genes.ALU-rich.txt /data/CCBR_Pipeliner/MAPLE/bed_files/hg19_protein-coding_genes.bed | awk -v OFS='\t' '{print $1,$2,$3}'> final.bed; head final.bed
fi

# create bed file
if [[ $flag_stage == "create_bed" ]]; then
    gene_list=/data/CCBR_Pipeliner/Pipelines/MAPLE/gene_list/$2
    final_bed=/data/CCBR_Pipeliner/Pipelines/MAPLE/bed_files/$3

    # clean file
    #sed -e "s/\r//g" $gene_list > $gene_list.tmp

    #create bed file
    awk -v OFS='\t' '{print $1,$2,$3}' $gene_list > $final_bed

    #preview and cleanup
    head $final_bed
    wc -l $gene_list
    wc -l $final_bed
   
    #rm $gene_list.tmp 
fi

#########################################################################################################################
# prep fastq
#########################################################################################################################
sub_id="kid_batch5through8"
# ### rawdata / links
if [[ $flag_stage == "create_links" ]]; then
    project_id="MAPLE";\
    raw_dir="/data/CCBR/rawdata/$project_id/$sub_id"; \
    
    proj_dir="/data/Zhurkin-20/rawdata/$sub_id"
    if [[ ! -d $proj_dir ]]; then mkdir -p $proj_dir; fi
    
    # create links
    for f in $raw_dir/*/*/*fastq.gz; do \
        ln -s $f $proj_dir/$(basename "$f");\
    done;\
    ls $proj_dir/*

    for f in $proj_dir/*.gz; do \
        new_name=`echo $f | sed -s "s/_S[0-9][0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_S[0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_001//g"`
        new_name=`echo $new_name | sed -s "s/_R/.R/g"`
        mv $f $new_name
    done;\
    ls $proj_dir/*

    # #create rename file
    #col1 is old and col2 is new
    ls $proj_dir/*fastq.gz | cut -f6 -d "/" > $proj_dir/file_rename.csv
fi

# ### reorganize data
if [[ $flag_stage == "reorg_fq" ]]; then
    project_id="ccbr1214";\
    raw_dir="/data/CCBR/rawdata/$project_id/$sub_id"; \
    
    proj_dir="/data/Zhurkin-20/rawdata/$sub_id"
    if [[ ! -d $proj_dir ]]; then mkdir -p $proj_dir; fi

    #mv files
    for f in $raw_dir/*/*/*; do
        mv $f $proj_dir/$(basename "$f")
    done
    ls $proj_dir/*

    for f in $proj_dir/*.gz; do \
        new_name=`echo $f | sed -s "s/_S[0-9][0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_S[0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_001//g"`
        new_name=`echo $new_name | sed -s "s/_R/.R/g"`
        mv $f $new_name
    done;\
    ls $proj_dir/*

    # #create rename file
    #col1 is old and col2 is new
    ls $proj_dir/*fastq.gz | cut -f6 -d "/" > $proj_dir/file_rename.csv
fi

if [[ $flag_stage == "rename" ]]; then
    proj_dir="/data/Zhurkin-20/rawdata/$sub_id"; \

    #rename from file col1 is old and col2 is new
    echo "" >> $proj_dir/file_rename.csv
    sed -e "s/\r//g" $proj_dir/file_rename.csv > $proj_dir/file_rename_clean.csv;\
    while IFS=',' read -a files 
    do
        mv "$proj_dir/${files[0]}" "$proj_dir/${files[1]}"
    done < $proj_dir/file_rename_clean.csv; \
    
    rm $proj_dir/file_rename.csv; \
    ls $proj_dir/*
fi

#########################################################################################################################
# pipeline init, dryrun, run
#########################################################################################################################
# pipeline prep
#sample_list=("1_Kid_norm_5u" "2_Kid_norm_5_200u" "3_Kid_tumor_5u" "4_Kid_tumor_5_200")
# sample_list=("1_Kid2_norm_5u" "2_Kid2_norm_5_200u" "3_Kid2_tumor_5u" "4_Kid2_tumor_5_200")
# sample_list=("5_Kid3_norm_5u" "6_Kid3_norm_5_200u" "7_Kid3_tumor_5u" "8_Kid3_tumor_5_200")
# sample_list("10_Kid4_tum1_5_200_20" "11_Kid4_tum2_5u_20" "12_Kid4_tum2_5_200_20" "13_Kid4_norm_1u_10" 
# sample_list=("1_K5_nor_5u_5" "2_K5_nor_5_200_5" "3_K5_tum_5u_5" "4_K5_tum_5_200_5" "5_K6_nor_5u_5")
# sample_list=( "6_K6_nor_5u" "7_K6_nor_5_200" "8_K6_tum_5u_5" "9_K6_tum_5u" "10_K6_tum_5_200")
# sample_list=("11_K7_nor_5u" "12_K7_nor_5_200" "13_K7_tum_5u" "14_K7_tum_5_200")
# sample_list=("15_K8_nor_5u" "16_K8_nor_5_200"  "17_K8_tum_5u" "18_K8_tum_5_200")
# sample_list=("1_Kid4_norm_5u_10" "2_Kid4_norm_5_200_10" "3_Kid4_tum1_5u_10" "4_Kid4_tum1_5_200_10" "1_K5_nor_5u_5" "2_K5_nor_5_200_5" "3_K5_tum_5u_5" "4_K5_tum_5_200_5" "15_K8_nor_5u" "16_K8_nor_5_200" "17_K8_tum_5u" "18_K8_tum_5_200")
sample_list=("1_Kid4_norm_5u_10")

#contrast_list=("1000_ALUrich" "1000_ALUrich" "2500_ALUdepleted" "GENES_2000_ALU" "GENES_2000_AT" "GENES_2000_GC" "GENES_2000_Length" "proteinCoding")
#contrast_list=("CDK11A" "CDK11B")
#contrast_list=("TP73" "P58")
contrast_list=("TTC34")

range="140-160"
min="140"
max="160"

lim="1000000"

bed_list_old="bed_lists_230208"
bed_list_new="bed_lists_230212"

############################################################################
# PIPELINE TASKS
############################################################################
if [[ $flag_stage == "prep_config" ]]; then

    sample_list="${sample_list[0]} ${sample_list[4]}"
    for f in ${sample_list[@]}; do
        echo "**$f**"

        # set dir
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        config_file="$analysis_dir/config.yaml"
        
        # update bp lengths
        sed -i "s/fragment_length_min: \\\"120\\\"/fragment_length_min: \\\"140\\\"/g" $config_file
        sed -i "s/fragment_length_max: \\\"140\\\"/fragment_length_max: \\\"160\\\"/g" $config_file

        # second pass
        sed -i "s/pipeline_phase: \\\"first_pass\\\"/pipeline_phase: \\\"second_pass\\\"/g" $config_file

        # set master_table info
        sed -i "s/master_table: \\\"Y\\\"/master_table: \\\"N\\\"/g" $config_file
        sed -i "s/data\\/CCBR_Pipeliner\\/Pipelines\\/ccbr1214\\/bed_files\\/hg19_protein-coding_genes.bed/data\\/Zhurkin-20\\/analysis\\/$f\\/resources\\/bed_lists_alu.csv/g" $config_file

        # bed list
        cp /home/sevillas2/git/MAPLE/resources/bed_lists_alu.csv $analysis_dir/resources/

        # verify
        cat $config_file | grep "fragment_length_min"
        cat $config_file | grep "fragment_length_max"
        cat $config_file | grep "master_table:" | grep -v "#"
        cat $config_file | grep "bed_list_name:"
        cat $config_file | grep "pipeline_phase"

        echo -e "\n\n"
    done
fi

# pipeline_init
if [[ $flag_stage == "pipe_init" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # initialize
        analysis_dir="/data/Zhurkin-20/analysis/$f"        
        if [[ ! -d $analysis_dir ]]; then
           /home/sevillas2/git/MAPLE/run_maple --runmode=init --workdir=$analysis_dir
        fi
        
        # set type
        type=`echo $f | cut -f2 -d"_"`
        if [[ $type == "Bst" ]]; then
            type="bst"
        elif [[ $type == "NB26" ]]; then
            type="nb26"
        elif [[ $type == "Kid" ]]; then
            type="kid"
        elif [[ $type == "Kid2" || $type == "Kid3" ]]; then
            type="kid_batch2"
        elif [[ $type == "Kid4" ]]; then
            type="kid_batch4"
        elif [[ $type == "K5" || $type == "K6" || $type == "K7" || $type == "K8" ]]; then
            type="kid_batch5through8"
        else
            type="rwpe"
        fi

        # set fqs
        fq_dir="/data/Zhurkin-20/rawdata/$type"
        fq1=$fq_dir/$f.R1.fastq.gz
        fq2=$fq_dir/$f.R2.fastq.gz

        echo -e "sampleName\ttype\tpath_to_R1_fastq\tpath_to_R2_fastq" > $analysis_dir/manifests/samples.tsv
        echo -e "$f\tkidney\t$fq1\t$fq2" >> $analysis_dir/manifests/samples.tsv
        cat $analysis_dir/manifests/samples.tsv
    done
fi

# copy the scripts and resources file to the workding dir
if [[ $flag_stage == "pipe_second_init" ]]; then
    sample_list="${sample_list[0]} ${sample_list[4]}"
    for f in ${sample_list[@]}; do
        echo "--$f"

        analysis_dir="/data/Zhurkin-20/analysis/$f"        
        cp /home/sevillas2/git/MAPLE/resources/* $analysis_dir/resources/
        cp /home/sevillas2/git/MAPLE/workflow/scripts/* $analysis_dir/scripts/
    done
fi

# pipeine dryrun
if [[ $flag_stage == "pipe_dry" ]]; then
    sample_list="${sample_list[0]} ${sample_list[4]}"
    for f in ${sample_list[@]}; do
        echo "--$f"

        # dry run
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        /home/sevillas2/git/MAPLE/run_maple --runmode=dryrun --workdir=$analysis_dir/
    done
fi

# run the pipeline
if [[ $flag_stage == "pipe_run" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        /home/sevillas2/git/MAPLE/run_maple --runmode=run --workdir=$analysis_dir/
    done
fi

# run local
if [[ $flag_stage == "pipe_local" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        /home/sevillas2/git/MAPLE/run_maple --runmode=local --workdir=$analysis_dir/
    done
fi

# unlock
if [[ $flag_stage == "pipe_unlock" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        /home/sevillas2/git/MAPLE/run_maple --runmode=unlock --workdir=$analysis_dir/
    done
fi

############################################################################
# PREP TASKS
############################################################################
# move from zhurkin to rawdata
if [[ $flag_stage == "mv_samples" ]]; then
    for f in ${project_list[@]}; do
        proj=`echo $f | cut -f3 -d "_"`;\
        if [[ $proj == "ARPE" ]]; then
            id="arpe";\
        else
            id="mel";\
        fi

        mkdir -p /data/CCBR/rawdata/$id/$f;\
        mv /data/Zhurkin-20/rawdata/$f/* /data/CCBR/rawdata/$id/$f;\
    done
fi
############################################################################ second pass
# copies bed list to resource/dir for each sample; updates config to use this file
if [[ $flag_stage == "prep_bedlist" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"
        master_list="/data/Zhurkin-20/analysis/${sample_list[0]}/resources/$bed_list_new.csv"
        resource_dir="/data/Zhurkin-20/analysis/$f/resources"
        config_file="/data/Zhurkin-20/analysis/$f/config.yaml"
        cp $master_list $resource_dir

        sed -i "s/$bed_list_old.csv/$bed_list_new.csv/" $config_file
        sed -i "s/pipeline_phase: \"first_pass\"/pipeline_phase: \"second_pass\"/" $config_file
        sed -i "s/pipeline_phase: \"third_pass\"/pipeline_phase: \"second_pass\"/" $config_file
        sed -i "s/fragment_length_min: \"1[0-9]0\"/fragment_length_min: \"$min\"/" $config_file
        sed -i "s/fragment_length_max: \"1[0-9]0\"/fragment_length_max: \"$max\"/" $config_file
        sed -i "s/master_table: \"Y\"/master_table: \"N\"/" $config_file
    done
fi
############################################################################ third pass
# prepare contrast files in first samples /manifest/dir
if [[ $flag_stage == "prep_contrasts" ]]; then
    for c in ${contrast_list[@]}; do
        contrast_manifest="/data/Zhurkin-20/analysis/${sample_list[0]}/manifests/contrast_$c.tsv"
        if [[ -f $contrast_manifest ]]; then rm $contrast_manifest; fi
        echo "DAC_files" > $contrast_manifest
        
        for f in ${sample_list[@]}; do
            echo "/data/Zhurkin-20/analysis/$f/results/04_dyads/03_CSV/$f.hg19.$range.lim$lim.$c.DAC.corrected.csv" >> $contrast_manifest
        done

        cat $contrast_manifest
    done
fi

# run contrasts
if [[ $flag_stage == "run_contrasts" ]]; then

    for c in ${contrast_list[@]}; do
        contrast_manifest="/data/Zhurkin-20/analysis/${sample_list[0]}/manifests/contrast_$c.tsv"
        if [[ -f $contrast_manifest ]]; then
            sed -i "s/second_pass/third_pass/g" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
            sed -i "s/CONTRAST_FILL_IN/contrast_$c/g" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
            sed -i "s/SHORTHAND_FILL_IN/$c/" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
            
            ../.././run_maple --runmode=runlocal --workdir=/vf/users/Zhurkin-20/analysis/${sample_list[0]}

            sed -i "s/contrast_$c/CONTRAST_FILL_IN/g" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
            sed -i "s/$c/SHORTHAND_FILL_IN/" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
        fi
    done
fi

############################################################################
# OTHER TASKS
############################################################################
if [[ $flag_stage == "cleanup" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # dir
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        rm $analysis_dir/results/04_dyads/04_master_table/*gene_list*
    done
fi

if [[ $flag_stage == "check" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # check
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        
        ls $analysis_dir/results/04_dyads/03_CSV/
        head $analysis_dir/results/04_dyads/03_CSV/*corrected*
    done
fi

if [[ $flag_stage == "check_master" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # check
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        
        # master table checks
        head -n3 $analysis_dir/results/04_dyads/04*/*120-140*master_table* | awk -F"," '{print $1,$2,$3}'
        cat $analysis_dir/results/04_dyads/04*/*120-140*master_table* | wc -l
        awk -F"," '{print NF}' $analysis_dir/results/04_dyads/04*/*120-140*master_table*  | sort -nu | tail -n 1
    done
fi

if [[ $flag_stage == "update_scripts" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # copy scripts
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        cp * $analysis_dir/scripts
    done
fi

if [[ $flag_stage == "tabtocomma" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # check
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        
        sed -i "s/\t/,/g" $analysis_dir/results/04_dyads/04_master_table/$f.hg19.140-160.lim1000000.master_table.DAC.corrected.csv
    done
fi

if [[ $flag_stage == "output" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # check
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        
        for f in $analysis_dir/results/04_dyads/03_CSV/*alu*; do
            echo "----$f"
        done
    done
fi