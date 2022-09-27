# 1 output log file prefix
# 2 samples manifest
# 3 contrast manifest
# rm /data/sevillas2/ccbr1214/test/logs/manifest_*; python workflow/scripts/check_manifest.py /data/sevillas2/ccbr1214/test/logs/manifest_log_ config/samples.tsv config/contrasts.tsv; cat /data/sevillas2/ccbr1214/test/logs/manifest_log_*

import pandas as pd
import re
import sys
from datetime import datetime


def check_header(input_df,expected_list,check_file,error_log):
    check =  all(item in list(input_df.columns.values) for item in expected_list)

    if check is False:
        error_log.append("The file {} does not include all required elements:{}\n".format(check_file,set(expected_list).difference(set(list(input_df.columns.values)))))

    return(error_log)

def check_samples(input_df,error_log):
    for select_cols in input_df.columns.values:
        tmp_list = input_df[select_cols].tolist()

        #check if all cols are unique; except type
        if select_cols != "type":
            if len(tmp_list) != len(set(tmp_list)):
                error_log.append("File names must be unique: {}\n".format(set([x for x in tmp_list if tmp_list.count(x) > 1])))
        
        #file_name
        if select_cols == "path_to_R1_fastq" or select_cols == "path_to_R2_fastq":
            #check all filenames end in .fastq.gz
            for item in tmp_list:
                if not (item.endswith('.fastq.gz')):
                    error_log.append("All items in file_name column must end in fastq.gz - please update filename: {}\n".format(item))

        #sampleName
        if select_cols == "sampleName":
            #check values are alpha/numeric or _
            regex = re.compile(r'[A-Z]_[a-z][0-9]')
            for item in tmp_list:
                if(regex.search(item) != None):
                    error_log.append("{} values can only contain alpha/numeric values or _ - please review: {}\n".format(select_cols,item))
    return(error_log)

def check_contrasts(input_df,error_log):
    #check if there are any NA values, otherwise check alpha/num
    check_na = input_df.isnull().values.any()

    if check_na == False:
      for select_cols in input_df.columns.values:
          tmp_list = input_df[select_cols].tolist()
          #check values are alpha/numeric or _
          regex = re.compile(r'[A-Z]_[a-z][0-9]')
          for item in tmp_list:
            if(regex.search(item) != None):
              error_log.append("{} values can only contain alpha/numeric values or _ - please review: {}".format(select_cols,item))
    else:
      error_log.append("Contrast manifest cannot contain blanks or NA's on any line. Please review")
    return(error_log)

def check_sample_to_contrast(input_df_s,input_df_c,error_log,de_log):
    #create list of samples
    ls = list(input_df_s["sampleName"].unique())

    # set error counter
    check=0

    # for each contrast, check it's sample name / fastq files were given
    for select_cols in input_df_c.columns.values:
        if select_cols != "sampleName":
            #create list of requested contrasts
            lc = list(input_df_c[select_cols].unique())

            #make sure all contrasts are samples found in workflow
            l_ctos = list(set(lc) - set(ls))
            
            #if any samples are missing, print error
            if len(l_ctos) !=0:
                error_log.append("The following sample(s) was/were not found in the sample_manifest tsv but were found in the contrasts manifest: {}. Please review manifests.".format(l_ctos))
                check=check+1

    # check row contains unique values
    for index, row in input_df_c.iterrows():
        if(len(set(row)) != len(row)):
            error_log.append("Duplicate contrasts were found in the contrast manifest. Each row must contain unique sample names. Please review and update the manifest {}.".format(row))
            check=check+1

    if check == 0:
        de_log.append("The following sample to sample comparisons will be run:\n{}.".format(input_df_c))

    return (error_log, de_log)


#create log
error_log = []
de_log = []

#Check samples file
check_file = sys.argv[2]
s_df = pd.read_csv(check_file,sep="\t")
s_req = ['sampleName','type','path_to_R1_fastq', 'path_to_R2_fastq']
error_log = check_header(s_df,s_req,check_file,error_log)
error_log = check_samples(s_df,error_log)

#Check contrast file
# will check for the minimum of two contrasts, but more may be added
check_file = sys.argv[3]
c_df = pd.read_csv(check_file,sep="\t")
c_req = ['contrast1','contrast2']
error_log = check_header(c_df,c_req,check_file,error_log)
error_log = check_contrasts(c_df,error_log)
      
#Check concordance between sample and contrast files
error_log,de_log = check_sample_to_contrast(s_df,c_df,error_log,de_log)

# if there is no error, create empty file
# if there is an error, create error log
if len(error_log)==0:
    #create empty no_error file
    new_path = str(sys.argv[1]) + "no_errors.txt"
    open(new_path,"w+").write('\n'.join(error_log))
else:
    date_str = '%s%s_%s' % (sys.argv[1],'contains_errors',datetime.now().strftime('%Y%m%d'))
    new_path =  date_str + '.txt'
    open(new_path, 'w+').write('\n'.join(error_log))

#print out DE log
if len(de_log)!=0:
    new_path = str(sys.argv[1]) + "manifest_log_comparisons.txt"
    open(new_path,"w+").write('\n'.join(de_log))


