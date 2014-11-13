import subprocess as sub
import sys
from brendanlib.grifflib import getcurrentjobs

currentjobs,jobids,jobstatus = getcurrentjobs()

argin = sys.argv[1]

for jobid,jobname,status in zip(jobids,currentjobs,jobstatus):
    if argin in jobname and status == 'PD':
        cmd_scancel = "scancel " + str(jobid)
	sub.call([cmd_scancel],shell=True)
        print cmd_scancel


