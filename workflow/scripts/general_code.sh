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

    reverse engineer gene list
    awk 'FNR==NR{arr[$1,$2];next} (($1,$2) in arr)' 1000.Genes.ALU-rich.bed /data/CCBR_Pipeliner/ccbr1214/bed_files/hg19_protein-coding_genes.bed > 1000.Genes.ALU-rich.txt

    create bed file
    grep -Fwf 1000.Genes.ALU-rich.txt /data/CCBR_Pipeliner/ccbr1214/bed_files/hg19_protein-coding_genes.bed | awk -v OFS='\t' '{print $1,$2,$3}'> final.bed; head final.bed
fi

# create bed file
if [[ $flag_stage == "create_bed" ]]; then
    gene_list=$2
    final_bed=$3

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

#prep samples
sample_list=("1_Kid_norm_5u" "2_Kid_norm_5_200u" "3_Kid_tumor_5u" "4_Kid_tumor_5_200")
if [[ $flag_stage == "pipe_prep" ]]; then
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

# ### rawdata / links
if [[ $flag_stage == "create_links" ]]; then
    project_id="ccbr1214";\
    sub_id="kid"
    raw_dir="/data/CCBR/rawdata/$project_id/$sub_id"; \
    proj_dir="/data/Zhurkin-20/rawdata/$sub_id"; \
    if [[ ! -d $proj_dir ]]; then mkdir -p $proj_dir; fi; \
    for f in $raw_dir/*/*fastq.gz; do \
        ln -s $f $proj_dir/$(basename "$f");\
    done;\
    ls $proj_dir/*

    project_id="ccbr1214";\
    for f in $proj_dir/*; do \
        new_name=`echo $f | sed -s "s/_S[0-9][0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_S[0-9]//g"`
        new_name=`echo $new_name | sed -s "s/_001//g"`
        new_name=`echo $new_name | sed -s "s/_R/.R/g"`
        mv $f $new_name
    done;\
    ls $proj_dir/*
fi