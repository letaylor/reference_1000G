#!/usr/bin/env python3

# submits slurm jobs

import os
import sys
import math
#from itertools import izip # for py2.7 and old method
from snakemake.utils import read_job_properties

# below is a snakemake way to do this
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)


# get the parameters for the command ###########################################
cluster = job_properties.get('cluster', {})

# specify the job name
job_name = cluster.get('name', 'snkmk') # default = snkmk

# what environemnt variables from current envi do you wish to export
# common options would be ALL or NONE
export = cluster.get('export', 'ALL') # default is export all variables

# pass customized flags (e.g., -exclusive, --gres=gpu)
# * gres = select a generic resource. e.g., --gres=gpu:2 means request use of 2
#   GPUs on compute nodes
custom_flags = cluster.get('custom_flags', None)

# partition flags define a group of nodes (e.g., stanard nodes, gpu nodes)
# think of these as different job queues for a group of nodes
partition = cluster.get('partition', None)

# specify the number of nodes
# note that by adding the -exclusive you can make your job the only job running
# on a particular nodes
#nnodes = cluster.get('nodes', 1) # default = 1
nnodes = None # node count required for job ... can set min-max

# specify the number of number tasks
# * "task" in this context can be thought of as a "process" or an instance of a
#   running program.
# * a multi-process program (MPI) is comprised of multiple tasks. Worded
#   another way if your program supports communication across computers or you
#   plan on running independent tasks in parallel, request multiple tasks.
# * by contrast to MPI, a multi-threaded program is comprised of a single task,
#   which can in turn use multiple CPUs (this is set via cpus flag).
# * note that, even within the same job, multiple tasks do not necessarily run
#   on a single node.
ntasks = cluster.get('tasks', 1) # tasks = 1

# specify the number of number of cpus per task
# * if cpu without hypterthreading, this is equivalent to number of threads.
# * to request cpus for your multi-threaded program, use the --cpus-per-task
#   flag. Individual tasks cannot be split across multiple compute nodes, so
#   requesting a number of CPUs with --cpus-per-task flag will always result in
#   all your CPUs allocated on the same compute node.
ncpus = cluster.get('cpus', 1) # tasks = 1

# specify the total memory requirements in Gb
# --mem = memory limit per compute node (MB) for the  job.
#         Do not use with mem-per-cpu flag.
# --mem-per-cpu = memory required per allocated cpu (MB)
memory = cluster.get('memory', 10) # default to 10Gb
# mem is specified as total memory in Gb, but we
# need to give it to slurm as memory per cpu in Mb
memory = int(math.floor(float(memory * 1e9) / (ncpus * 1e6)))

# get the std out / err file
# this could be within several directories
# it is recommended that the the job number (%j) and the first node (%N) added
# example myjob-%j-%N.o myjob-%j-%N.e
output = cluster.get('output', '{}-%j-%N.o'.format(job_name))
if output:
    output = os.path.realpath(output)
    tdir = os.path.dirname(output)
    if not os.path.exists(tdir): os.makedirs(tdir)
error = cluster.get('error', '{}-%j-%N.e'.format(job_name))
if error:
    error = os.path.realpath(error)
    tdir = os.path.dirname(error)
    if not os.path.exists(tdir): os.makedirs(tdir)
################################################################################


# build the command ############################################################
cmd = 'sbatch'
# add standard commands
cmd = '{} --job-name="{}" --export={} --ntasks={:d}'.format(
    cmd,
    job_name,
    export,
    ntasks)
# add cpu information
cmd = '{} --cpus-per-task={:d} --mem-per-cpu={:d}'.format(
    cmd,
    ncpus,
    memory)
# add optional flags
if custom_flags:
    cmd = '{} {}'.format(cmd, custom_flags)
if partition:
    cmd = '{} --partition={}'.format(cmd, partition)
if nnodes:
    cmd = '{} --nodes={:d}'.format(cmd, nnodes)
if output:
    cmd = '{} --output={}'.format(cmd, output)
if error:
    cmd = '{} --error={}'.format(cmd, error)
cmd = '{} {}'.format(cmd, jobscript) # could add --wrap=<command string>
################################################################################


# run the command ##############################################################
os.system(cmd)
################################################################################
