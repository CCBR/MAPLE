## A_SATS_REMOVED
# code changes tracked: https://github.com/CCBR/samantha_log/issues/184

flag=$1

#120-140
# sample_list=("2_Kid4_norm_5_200_10" "4_Kid4_tum1_5_200_10" "2_K5_nor_5_200_5" "4_K5_tum_5_200_5" "16_K8_nor_5_200" "18_K8_tum_5_200")

#140-160
# sample_list=("1_Kid4_norm_5u_10" "3_Kid4_tum1_5u_10" "1_K5_nor_5u_5" "3_K5_tum_5u_5" "15_K8_nor_5u" "17_K8_tum_5u")

sample_list=("1_Kid4_norm_5u_10")
min_frag="140"
max_frag="160"
limit=100000
max_frag_dist=1100

sats_bed="/data/CCBR_Pipeliner/Pipelines/ccbr1214/bed_files/a_sats_only.bed"
pipeline_dir="/home/sevillas2/git/MAPLE/workflow/scripts"

# base_dir="/data/Zhurkin-20/analysis/"
base_dir="/data/sevillas2/victor"


if [[ $flag == "check_samples" ]]; then
    for f in ${sample_list[@]}; do
        echo "**$f**"
        
        analysis_dir="$base_dir/$f/results"
        bed_dir="$analysis_dir/03_aligned/02_bed"
        mapped_bed="$bed_dir/$f.hg19.mapped.bed"

        ls -la $mapped_bed
    done
fi

if [[ $flag == "create_sh" ]]; then
    for f in ${sample_list[@]}; do
        echo "**$f**"

        # set dir
        analysis_dir="$base_dir/$f/results"
        tmp_dir="$analysis_dir/a_sat_tmp"
        bed_dir="$analysis_dir/03_aligned/02_bed"
        dyad_dir="$analysis_dir/04_dyads/01_DYADs"
        hist_dir="$analysis_dir/04_dyads/02_histograms"
        csv_dir="$analysis_dir/04_dyads/03_CSV"
        script_dir="$base_dir/$f/scripts"
        if [[ ! -d $tmp_dir ]]; then mkdir -p $tmp_dir; fi

        # set files
        sh_file="$script_dir/$f.asatsremoved.$min_frag-$max_frag.sh"

        mapped_bed="$bed_dir/$f.hg19.mapped.bed"
        sats_removed_bed="$bed_dir/$f.hg19.sats_removed.bed"
        sats_removed_sub_bed="$bed_dir/$f.hg19.${min_frag}-${max_frag}.lim${limit}.sats_removed.bed"

        dyads="$dyad_dir/$f.hg19.$min_frag-$max_frag.sats_removed.DYADs"
        dyads_sorted="$dyad_dir/$f.hg19.$min_frag-$max_frag.sats_removed.DYADs.sorted"
        dyads_hist="$hist_dir/$f.hg19.$min_frag-$max_frag.sats_removed.DYADs.hist"
        
        csv="$csv_dir/$f.hg19.$min_frag-$max_frag.lim1000000.sats_removed.DAC.csv"

        echo "#!/bin/sh
        
        module load bedtools python
        
        # remove sats
        echo "--remove sats"
        bedtools intersect -v -a ${mapped_bed} -b ${sats_bed} > ${sats_removed_bed}
        awk -v min_frag="$min_frag" -v max_frag="$max_frag" '{ if (\$3-\$2 >= min_frag && \$3-\$2 <= max_frag) print \$0}' ${sats_removed_bed} > ${sats_removed_sub_bed}
        head ${sats_removed_sub_bed}

        # Find fragment centers (DYADs) and make histogram (Occurrences) ###########
        echo "--weighted pos"
        python3 $pipeline_dir/WeigthedDYADposition.py ${sats_removed_bed} ${dyads}
        sort -k1,1 -k2n,2 ${dyads} > ${dyads_sorted}
        head ${dyads_sorted}
        echo "--uniq pos"
        python $pipeline_dir/Uniq_Position.py ${dyads_sorted} ${dyads_hist}
        head ${dyads_hist}

        # Compute auto-correlation ###########
        echo "--DAC"
        python $pipeline_dir/DAC.py ${dyads_hist} ${limit} ${max_frag_dist} ${csv}

        head ${csv}" > $sh_file

    done
fi

if [[ $flag == "submit_sh_local" ]]; then
    for f in ${sample_list[@]}; do
        echo "**$f**"

        # set dir
        script_dir="$base_dir/$f/scripts"

        ## sh file
        sh_file="$script_dir/$f.asatsremoved.$min_frag-$max_frag.sh"
        
        sh $sh_file
    done
fi

if [[ $flag == "submit_sh_cluster" ]]; then
    for f in ${sample_list[@]}; do
        echo "**$f**"

        # set dir
        script_dir="$base_dir/$f/scripts"

        ## sh file
        sh_file="$script_dir/$f.asatsremoved.$min_frag-$max_frag.sh"
        
        # submit
        sbatch --cpus-per-task=32 --verbose \
        --output=$base_dir/$f/logs/%j.out \
        --mem=200g --gres=lscratch:450 --time 10:00:00 \
        --error=$base_dir/$f/logs/%j.err \
        $sh_file
    done
fi