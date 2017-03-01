import subprocess
import string
import random
import os
import time
from collections import defaultdict
from . import io

array_dict = {'lsf': '''#!/bin/bash
#BSUB -n 1
#BSUB -R "rusage[mem={c[mem]}]"
#BSUB -J {c[job]}[1-{c[job_no]}]
#BSUB -o wrf.%I_tmp.out
#BSUB -e wrf.%I_tmp.out
#BSUB -W {c[time]}

cd {c[home_dir]}

name=$(sed -n "$LSB_JOBINDEX"p {c[namelist_file]})

{c[command]}
''', 'slurm': '''#!/bin/bash
#SBATCH --mem-per-cpu={c[mem]}
#SBATCH -J {c[job]}
#SBATCH --array=1-{c[job_no]}
#SBATCH -t {c[time]}
#SBATCH -o array_job_out_%j_tmp.out
#SBATCH -e array_job_err_%j_tmp.out
#SBATCH -n 1
#SBATCH -p serial

cd {c[home_dir]}

name=$(sed -n "$SLURM_ARRAY_TASK_ID"p {c[namelist_file]})

{c[command]}
'''}


batch_dict = {'lsf': '''

''', 'slurm': '''#!/bin/bash -l
#SBATCH --mem-per-cpu={c[mem]}
#SBATCH -J {c[job]}
#SBATCH -t {c[time]}
#SBATCH -o {c[job]}_tmp_output.txt
#SBATCH -e {c[job]}_tmp_errors.txt
#SBATCH -n 1
#

cd {c[home_dir]}

{c[command]}
'''}



def generate_id(size=8):
    """ Generate random sequences of characters for temporary file names.
    """
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.choice(chars) for _ in range(size))


def make_split_dict(chunk_iter, no_splits):
    chunk_list = list(chunk_iter)
    chunk_len = len(chunk_list)
    chunk_size = int(chunk_len/no_splits)
    split_dict = defaultdict(list)
    for ix, chunk in enumerate(chunk_list):
        if ix % chunk_size == 0:
            chunk_id = generate_id()
        split_dict[chunk_id].append(chunk)
    return split_dict


def split_seqs(seq_file, no_splits):
    seqs = io.read_fasta(seq_file)
    split_dict = make_split_dict(seqs, no_splits)
    for key, val in split_dict.items():
        seq_name = key + "_tmp.fasta"
        io.write_fasta(val, seq_name)
    return list(split_dict.keys())


def run_batch_job(batch_command, scheduler='slurm', memory=2048, run_time='02:00:00', cleanup=True):
    user = subprocess.check_output('whoami', universal_newlines=True).strip()
    home_dir = os.getcwd()

    if isinstance(batch_command, str):
        job_name = generate_id()
        batch_info = {'mem': memory, 'job': job_name, 'time': run_time,
                    'home_dir': home_dir, 'command': batch_command}
        batch = batch_dict[scheduler].format(c=batch_info)
        batch_file_name = generate_id() + "_tmp.sh"
        with open(batch_file_name, "w") as f:
            for line in batch:
                f.write(line)
        if scheduler == 'slurm':
            subprocess.call(['sbatch', batch_file_name])

    elif isinstance(batch_command, list):
        for command in batch_command:
            job_name = generate_id()
            batch_info = {'mem': memory, 'job': job_name, 'time': run_time,
                          'home_dir': home_dir, 'command': command}
            batch = batch_dict[scheduler].format(c=batch_info)
            batch_file_name = generate_id() + "_tmp.sh"
            with open(batch_file_name, "w") as f:
                for line in batch:
                    f.write(line)
            if scheduler == 'slurm':
                subprocess.call(['sbatch', batch_file_name])

    time.sleep(10)
    while True:
        jobs = subprocess.check_output(['squeue', '-u', user], universal_newlines=True).split("\n")
        if len(jobs) == 2:
            break
        print("{} jobs left.".format(len(jobs) - 2))
        time.sleep(5)

    if cleanup:
        print("Cleaning up.")
        [os.remove(tmp_file) for tmp_file in os.listdir() if "tmp" in tmp_file]
    print("Done!")


def run_array_job(seqs, batch_command, post_command=None, no_splits=1000, scheduler='slurm', memory=2048, run_time='02:00:00', cleanup=True):
    job_name = generate_id()
    user = subprocess.check_output('whoami', universal_newlines=True).strip()
    namelist = job_name + "_tmp.namelist"
    seq_ids = split_seqs(seqs, no_splits)
    job_no = len(seq_ids)
    home_dir = os.getcwd()
    batch_info = {'mem': memory, 'job': job_name, 'job_no': job_no, 'time': run_time,
                  'home_dir': home_dir, 'namelist_file': namelist, 'command': batch_command}
    array = array_dict[scheduler].format(c=batch_info)
    array_file_name = generate_id() + "_tmp.sh"
    with open(array_file_name, "w") as f:
        for line in array:
            f.write(line)
    with open(namelist, "w") as f:
        for item in seq_ids:
            f.write(item + "_tmp.fasta\n")
    if scheduler == 'slurm':
        subprocess.call(['sbatch', array_file_name])
        print("A total of {} jobs.".format(job_no))
        time.sleep(10)
        while True:
            jobs = subprocess.check_output(['squeue', '-u', user],
                                           universal_newlines=True).split("\n")
            if len(jobs) == 2:
                break
            print("{} jobs left.".format(len(jobs) - 2))
            time.sleep(5)
    if post_command:
        print("Executing the post-batch command.")
        subprocess.call(post_command, shell=True)
    if cleanup:
        print("Cleaning up.")
        [os.remove(tmp_file) for tmp_file in os.listdir() if "tmp" in tmp_file]
    print("Done!")
