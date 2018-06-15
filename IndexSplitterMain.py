#!/usr/bin/env python

# This script would parse the input directory to idenify eligible flowcell ids, and then construct the 
# eligible input files and from there would submit index_splitter for each of the flowcell ids

import argparse
import shutil
import re
import os
import os.path
from subprocess import call

parser = argparse.ArgumentParser(description = "Process inputs for index splitter", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir', '-i', dest = 'indir', type = str, required = True, help = "Input directory")
parser.add_argument('--outdir', '-o', dest = 'outdir', type = str, required = True, help = "Output directory")
parser.add_argument('--dict1_file', dest = 'dict1file', type = str, required = True, help = "index 1 dict file")
parser.add_argument('--dict2_file', dest = 'dict2file', type = str, required = True, help = "index 2 dict file")
parser.add_argument('--suffix_s1', dest = "suffix_s1", required = False, default = ".1.fastq")
parser.add_argument('--suffix_s2', dest = "suffix_s2", required = False, default = ".2.fastq")

parser.add_argument("-m", "--memory", type=int, default=8000,
                    help="Amount of memory to allocate")
parser.add_argument('--no_qsub', dest = 'use_qsub', action = 'store_false', default = True, help = 'Does not submit qsub jobs.' )
parser.add_argument('--total_run', dest = 'total_run', type = int, default = -1, help = 'Run for top read_count number of reads (if -1, run for all reads)' )
parser.add_argument('--run_time_ind', dest = 'run_time_ind', type = int, default = 24, help = 'Run time for individual jobs (default: 24 hours)' )

args = parser.parse_args()

indir = args.indir
outdir = args.outdir
dict1file = args.dict1file
dict2file = args.dict2file
memory = args.memory
use_qsub = args.use_qsub
total_run = args.total_run
suffix_s1 = args.suffix_s1
suffix_s2 = args.suffix_s2
run_time_ind = args.run_time_ind

ldelim = '/'

print("Input dir: " + indir)
print("Output dir: " + outdir)
print("index 1 dict file" + dict1file)
print("index 2 dict file" + dict2file)
print("total allocated memory " + str(memory))

allowed_mem = int(memory /2.5)

# Get all the files that ends with fastq.gz

prefix_set = set()

for filename in os.listdir(indir):
    lsuffix_s1 = suffix_s1
    lsuffix_s2 = suffix_s2
    if filename.endswith("gz"):
        lsuffix_s1 += ".gz"
        lsuffix_s2 += ".gz"
    lprefix = ''
    if filename.endswith(lsuffix_s1):
        lprefix = filename.replace(lsuffix_s1, '')
    elif filename.endswith(lsuffix_s2):
        lprefix = filename.replace(lsuffix_s2, '')
    if lprefix and lprefix not in prefix_set:
        prefix_set.add(lprefix)

print(".........")
print(prefix_set)
print(".........")


# For each element in the prefix_set submit a job to the cluster

Script_dir = os.path.dirname(os.path.realpath(__file__))
index_split_cpp = Script_dir + ldelim + "dual_bc_splitter"

if not os.path.exists(outdir):
    os.makedirs(outdir)
joblist_path = outdir + ldelim + "index_bcsplit_joblist.txt"
jfile = open(joblist_path, "w")

def refine_gz_file(gz_file):
    if not os.path.isfile(gz_file):
        only_file = gz_file.replace(".gz", "")
        if not os.path.isfile(only_file):
            raise IOError('File not found: ' + gz_file)
        else:
            return only_file
    else: 
        return gz_file


for prefix in prefix_set:
    parts = prefix.split("_")
    lane = parts[0]

    index1_file_str = indir + ldelim + prefix +   ".barcode_1.fastq.gz"
    index1_file = refine_gz_file(index1_file_str)
    index2_file_str = indir + ldelim + prefix +   ".barcode_2.fastq.gz"
    index2_file = refine_gz_file(index2_file_str)
    read1_file_str = indir + ldelim + prefix +   ".1.fastq.gz"
    read1_file = refine_gz_file(read1_file_str)
    read2_file_str = indir + ldelim + prefix +  ".2.fastq.gz"
    read2_file = refine_gz_file(read2_file_str)

    out_log = outdir + ldelim + prefix + "." + lane + "_out.txt"
    err_log = outdir + ldelim + prefix + "." + lane + "_err.txt"
    lprefix = prefix + "." + lane
    print("prefix: " + lprefix)
    #print("index: " + index_file)
    #print("read1_file: " + read1_file)
    #print("read2_file: " + read2_file)
    job_str = index_split_cpp + " --bc1 " + dict1file + " --bc2 " + dict2file + \
        " --bc1_read " + index1_file + " --bc2_read " + index2_file + \
        " --read1 " + read1_file + " --read2 " + read2_file + " -p " + lprefix + \
        " -o " + outdir + " --allowed-mb " + str(allowed_mem) + " --total_run " + str(total_run) + " 1> " + out_log + " 2> " + err_log + "\n"
    print("job_str: " + job_str)
    jfile.write(job_str)
jfile.close()

print(joblist_path)

UGER_cbp = "/broad/IDP-Dx_work/nirmalya/tools/ugetools/UGE_SUBMISSIONS/UGER_cmd_batch_processor.py"
UGER_cbp_dir = outdir + ldelim + "UGER_cbp"

if use_qsub and os.path.exists(UGER_cbp_dir):
    shutil.rmtree(UGER_cbp_dir)

joblist_cmd = UGER_cbp + " --cmds_file " + joblist_path + \
                                " --batch_size 1" + \
                                " --queue long" + \
                                " --memory 16" + \
                                " --run_time " + str(run_time_ind) + \
                                " --tracking_dir " + UGER_cbp_dir + \
                                " --project_name broad --bash_header /broad/IDP-Dx_work/nirmalya/bash_header"

if use_qsub:
    call(joblist_cmd.split())    

