#!/bin/python

import os
import pandas as pd
import subprocess

def submitSlurm(cmd, dependencies=None):
    if dependencies != None:
        # execute
        p = subprocess.Popen(['sbatch', '--dependency=afterany:' + ':'.join(dependencies), cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    else:
        p = subprocess.Popen(['sbatch', cmd], stdout=subprocess.PIPE)
        out, err = p.communicate()
    # get slurm job id to set up job submission dependency
    job_id = str(out).split(' ')[-1].replace("\\n'", '')
    return(job_id)


def writeSlurm(slurm_file, partition, cmd, jobName, stime='6:00:00', nthreads=None, mem=None):
    '''
    Submit batch job to Yale Milgram cluster (SLURM)

    required arguments:
        - slurm_file        base filepath string for writing slurm command and output
                                (e.g. /gpfs/milgram/project/holmes/Open_Data/DATA_UKBIOBANK/REPOSITORY/MRI/100024/slurm/100024_hello_world_cmd)
        - partition         short/long/scavenge
                                (e.g. short)
        - nthreads          up to 28 cpus on a node
                                (e.g. 4)
        - cmd           full string of the command to be written/run
                            (e.g. module load Apps/FREESURFER/5.3.0\\n print('Hello World'))
    '''
    slurm_name = slurm_file + '_slurm.txt'
    slurm_out  = slurm_file + '_slurmOut.txt'
    slurm_file = open(slurm_name, "w")
    slurm_file.write('#!/bin/bash\n')
    slurm_file.write('#SBATCH --partition=' +  partition + '\n')
    slurm_file.write('#SBATCH --output=' + slurm_out + '\n')
    slurm_file.write('#SBATCH --nodes=1\n')
    if mem != None:
        slurm_file.write('#SBATCH --mem=' + str(mem) + '\n')
    if nthreads != None:
        slurm_file.write('#SBATCH --ntasks=1 --cpus-per-task=' + str(nthreads) + '\n')
    slurm_file.write('#SBATCH --job-name=' + jobName + '\n')
    slurm_file.write('#SBATCH --time=' + stime + '\n')
    slurm_file.write(str(cmd))
    slurm_file.close()
    subprocess.call(['chmod', '0770', slurm_name])
    return(slurm_name)