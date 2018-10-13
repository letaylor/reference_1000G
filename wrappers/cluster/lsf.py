#!/usr/bin/env python3

# submits lsf jobs

import os
import sys
import math
#from itertools import izip # for py2.7 and old method
from snakemake.utils import read_job_properties

# below is a snakemake way to do this 
jobscript = sys.argv[-1]
job_properties = read_job_properties(jobscript)


cluster = job_properties.get('cluster', {})

# pass customized flags (e.g., -P bigmem)
custom_flags = cluster.get('custom_flags', None)

group = cluster.get('group', 'snakemake') # job group
queue = cluster.get('queue', None) # job queue
job_name = cluster.get('name', 'snakemake') # job name

threads = job_properties.get('threads', 1)
memory = cluster.get('memory', 10) # default to 10Gb
if memory:
    # mem is specified as total memory in Gb, but we
    # need to give it to SGE as memory per thread in Mb
    memory = int(math.floor(float(memory * 1e9) / (threads * 1e6)))
resources = 'select[mem>{}] rusage[mem={}] span[hosts=1]'.format(
    memory, memory
)

output = cluster.get('output', None)
if output:
    output = os.path.realpath(output)
    tdir = os.path.dirname(output)
    if not os.path.exists(tdir): os.makedirs(tdir)
error = cluster.get('error', None)
if error:
    error = os.path.realpath(error)
    tdir = os.path.dirname(error)
    if not os.path.exists(tdir): os.makedirs(tdir)

# build the command
cmd = 'bsub'
if custom_flags:
    cmd = '{} {}'.format(cmd, custom_flags)
if queue: 
    cmd = '{} -q {}'.format(cmd, queue)
cmd = '{} -g {} -M {} -R "{}" -J "{}"'.format(
    cmd,
    group,
    memory,
    resources,
    job_name)
if threads > 1:
    # -n n_procs: The expected resource consumption of the task.
    # -R "select[nprocs>2]": Select hosts with > 2 REAL PROCESSORS (chips).
    cmd = '{} -E "export OMP_NUM_THREADS={}"'.format(cmd, threads)
    cmd = '{} -n {}'.format(cmd, threads)
    # an alternative is to do '-R "affinity[core({})]"'.format(threads)
if output: 
    cmd = '{} -oo {}'.format(cmd, output)
if error: 
    cmd = '{} -eo {}'.format(cmd, error)
cmd = '{} {}'.format(cmd, jobscript)

# run the command
os.system(cmd)


# old way listed below
# get input args. drop 1 value (path to this script) and last value (path
# to the bash script)
# args = sys.argv[1:-1]
#
#
# # make dictionary of input
# def grouped(iterable, n):
#     return zip(*[iter(iterable)]*n)
#
# kwargs = {}
# for x, y in grouped(args, 2):
#    kwargs[x] = y
#
#
# # check for log output. if so, make sure the log directory is made
# for key in ['-o', '-oo', '-e', '-eo']:
#     if key in kwargs:
#         tdir = os.path.realpath(os.path.dirname(kwargs[key]))
#         if not os.path.exists(tdir): os.makedirs(tdir)

