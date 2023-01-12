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
    awk 'FNR==NR{arr[$1,$2];next} (($1,$2) in arr)' 1000.Genes.ALU-rich.bed /data/CCBR_Pipeliner/ccbr1214/bed_files/hg19_protein-coding_genes.bed > 1000.Genes.ALU-rich.txt

    #create bed file
    grep -Fwf 1000.Genes.ALU-rich.txt /data/CCBR_Pipeliner/ccbr1214/bed_files/hg19_protein-coding_genes.bed | awk -v OFS='\t' '{print $1,$2,$3}'> final.bed; head final.bed
fi

# create bed file
if [[ $flag_stage == "create_bed" ]]; then
    gene_list=/data/CCBR_Pipeliner/Pipelines/ccbr1214/gene_list/$2
    final_bed=/data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/$3

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
# ### rawdata / links
if [[ $flag_stage == "create_links" ]]; then
    project_id="ccbr1214";\
    sub_id="kid_batch2"
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

if [[ $flag_stage == "rename" ]]; then
    project_id="ccbr1214";\
    sub_id="kid_batch2"
    raw_dir="/data/CCBR/rawdata/$project_id/$sub_id"; \
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
#sample_list=("1_Kid2_norm_5u" "2_Kid2_norm_5_200u" "3_Kid2_tumor_5u" "4_Kid2_tumor_5_200" "5_Kid3_norm_5u" "6_Kid3_norm_5_200u" "7_Kid3_tumor_5u" "8_Kid3_tumor_5_200")
#sample_list=("5_Kid3_norm_5u" "6_Kid3_norm_5_200u" "7_Kid3_tumor_5u" "8_Kid3_tumor_5_200")
#sample_list=("1_Kid2_norm_5u" "2_Kid2_norm_5_200u" "3_Kid2_tumor_5u" "4_Kid2_tumor_5_200" "5_Kid3_norm_5u" "6_Kid3_norm_5_200u" "7_Kid3_tumor_5u" "8_Kid3_tumor_5_200" "1_Kid_norm_5u" "2_Kid_norm_5_200u" "3_Kid_tumor_5u" "4_Kid_tumor_5_200")
sample_list=("5_Kid3_norm_5u" "6_Kid3_norm_5_200u" "7_Kid3_tumor_5u" "8_Kid3_tumor_5_200")

#contrast_list=("1000_ALUrich" "1000_ALUrich" "2500_ALUdepleted" "GENES_2000_ALU" "GENES_2000_AT" "GENES_2000_GC" "GENES_2000_Length" "proteinCoding")
#contrast_list=("CDK11A" "CDK11B")
#contrast_list=("TP73" "P58")
contrast_list=("TTC34")

range="160-180"
min="160"
max="180"

lim="1000000"

bed_list_old="bed_lists_221219"
bed_list_new="bed_lists_221223"

if [[ $flag_stage == "pipe_init" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"

        # initialize
        analysis_dir="/data/Zhurkin-20/analysis/$f"        
        if [[ ! -d $analysis_dir ]]; then
           /home/sevillas2/git/ccbr1214/run --runmode=init --workdir=$analysis_dir
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

# pipeine dryrun
if [[ $flag_stage == "pipe_dry" ]]; then
    sample_list="${sample_list[0]} ${sample_list[4]}"
    for f in ${sample_list[@]}; do
        echo "--$f"

        # dry run
        analysis_dir="/data/Zhurkin-20/analysis/$f"        
        /home/sevillas2/git/ccbr1214/run --runmode=dryrun --workdir=$analysis_dir/
    done
fi
# run the pipeline
if [[ $flag_stage == "pipe_run" ]]; then
    for f in ${sample_list[@]}; do
        echo "--$f"
        analysis_dir="/data/Zhurkin-20/analysis/$f"
        /home/sevillas2/git/ccbr1214/run --runmode=run --workdir=$analysis_dir/
    done
fi

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
            
            ../.././run --runmode=runlocal --workdir=/vf/users/Zhurkin-20/analysis/${sample_list[0]}

            sed -i "s/contrast_$c/CONTRAST_FILL_IN/g" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
            sed -i "s/$c/SHORTHAND_FILL_IN/" /vf/users/Zhurkin-20/analysis/${sample_list[0]}/config.yaml
        fi
    done
fi