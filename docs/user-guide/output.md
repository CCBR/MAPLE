#4. Expected Outputs
The following directories are created under the output_directory, dependent on the Pass of the pipeline

## First Pass (first_pass)

- 01_trim: this directory includes trimmed FASTQ files
- 02_assembled: this directory includes assembled FASTQ files
- 03_aligned: this directory includes aligned BAM files and BED files
    - 01_bam: BAM files after alignment
    - 02_bed: converted bed files
    - 03_histograms: histograms of bed files

## Second Pass (second_pass)
- 04_dyads: this directory contains DYAD calculated files
    - 01_DYADs: this includes direct DYAD calculations
    - 02_histograms: this includes histogram occurances
    - 03_CSV: this includes the occurance data in CSV format

## Third Pass (third_pass)
- /path/to/output/contrast: this includes the contrast file for each sample provided in the contrasts.tsv manifest

## All Passes
- log: this includes log files
    - [date of run]: the slurm output files of the pipeline sorted by pipeline start time; copies of config and manifest files used in this specific pipeline run; error reporting script