#!/usr/bin/env python
"""
To run:
e.g. python script_example.py --project "Mu2eII" --jobname "test" --input "../../PS/TS3/PStape.txt" --fclname "ts.fcl" --memory=4GB 
"""
import os
import sys
import glob
import subprocess
import shutil
import argparse
from optparse import OptionParser

base_dir = os.path.join("/mu2e/app/users/sophie/Mu2e-IITargets/Offline")
tarbull = "/pnfs/mu2e/resilient/users/sophie/gridexport/tmp.OETAZLyzLw/Code.tar.bz"

def clean():
  """
  remove old scripts
  """
  remove = ["rm", "-rf", "000"]
  subproc = subprocess.run(remove)
  return subproc
  
def generate_fcls_staged(project, jobname, version, inputs, fclname, memory, lifetime):
    """
    generates fcl for grid jobs
    """
    generate = ["generate_fcl", "--embed", str(fclname), "--description="+str(jobname), "--dsconf="+str(version),  "--merge-factor=10",  "--inputs="+str(inputs),] 
    subproc = subprocess.run(generate)
    return subproc
    
def loop(jobname):
  """
  Loop over gen_fcl output and store as tarbull
  """
  loop =  "ls -1 -d ??? | cut -c 1-2 | sort | uniq | while read NN; do echo $NN; tar -cjf fcllist"+str(jobname)+"_${NN}.bz2 ${NN}?; done"
  subproc = subprocess.run(loop,shell=True)
  return subproc
  
def grid_upload(jobname):
  """
  upload tarbull of job fcls to the pnfs server
  """
  upload = ["mv","fcllist"+str(jobname)+"_00.bz2","/pnfs/mu2e/scratch/users/sophie/"]
  subproc = subprocess.run(upload)
  return subproc

def send_grid_job(jobname, memory,lifetime, version):
  """
  upload tarbull of job fcls to the pnfs server
  """
  grid = "mu2eprodsys --transfer-all-files --dsconf="+str(version)+" --memory="+str(memory)+" --code="+str(tarbull)+" --fcllist=/pnfs/mu2e/scratch/users/sophie/fcllist"+str(jobname)+"_00.bz2 --expected-lifetime="+str(lifetime)
  subproc = subprocess.run(grid,shell=True)
  return subproc
   
def main(options, args):
  clean()
  generate_fcls_staged(options.project, options.jobname, options.version, options.inputs, options.fclname, options.memory, options.lifetime)
  loop(options.jobname)
  grid_upload(options.jobname)
  send_grid_job(options.jobname, options.memory, options.lifetime, options.version)
  
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
    parser.add_option('-p','--project', dest='project', default = 'Mu2eII',help='tag for project tarbull', metavar='Mu2eIIdir')
    parser.add_option('-n','--jobname', dest='jobname', default = 'name',help='code name for physics run', metavar='Mu2eIIdir')
    parser.add_option('-v','--version', dest='version', default = 'v0',help='version', metavar='Mu2eIIdir')
    parser.add_option('-i','--inputs', dest='inputs', default = 'file.txt',help='for staged input', metavar='Mu2eIIdir')
    parser.add_option('-f','--fclname', dest='fclname', default = 'name.fcl',help='fcl file name', metavar='Mu2eIIdir')
    parser.add_option('-m','--memory', dest='memory', default = '4GB',help='allocated memomry for grid job', metavar='Mu2eIIdir')
    parser.add_option('-l','--lifetime', dest='lifetime', default = '24h',help='lifetime for jobs on grid', metavar='Mu2eIIdir')
    (options, args) = parser.parse_args()
    main(options,args);
    print("Finished")
