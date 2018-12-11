from __future__ import print_function
from qwalk_objects.trialfunc import export_qwalk_trialfunc
import os
####################################################
class DMCWriter:
  def __init__(self,sys,trialfunc,timestep=0.01,nblock=100,tmoves=True,savetrace=True,averages=''):
    ''' Object for producing input into a VMC QWalk run. 
    Args:
      trialfunc (str): trial wavefunction section or object that can export_qwalk_trialfunc().
      sys (str): system section or object that can export_qwalk_sys().
      nblock (int): Number of blocks to attempt.
      averages (str): averages section in qwalk.
    '''
    self.sys=sys
    self.trialfunc=trialfunc
    self.nblock=nblock
    self.averages=averages
    self.timestep=timestep
    self.tmoves=tmoves
    self.savetrace=savetrace
    self.completed=False

  #-----------------------------------------------
  def qwalk_input(self,infile):
    ''' Generate QWalk input file and write to infile.'''
    assert self.trialfunc is not None, "Must specify trialfunc before asking for qwalk_input."
    assert self.sys is not None, "Must specify system before asking for qwalk_input."

    if type(self.sys) is not str:
      system=self.sys.export_qwalk_sys()
    else:
      system=self.sys
    if type(self.trialfunc) is not str:
      trialfunc=export_qwalk_trialfunc(self.trialfunc)
    else:
      trialfunc=self.trialfunc

    outlines=[
        "method { DMC",
        "  nblock %i"%(self.nblock),
        "  timestep %g"%self.timestep
      ]
    if self.tmoves:
      outlines+=['  tmoves']
    if self.savetrace:
      tracename = "  %s.trace"%infile
      outlines+=["  save_trace %s"%tracename]
      
    outlines+=['  '+line for line in self.averages.split('\n')]+[
        "}"
      ]
    outlines.append(system)
    outlines.append(trialfunc)

    with open(infile,'w') as f:
      f.write('\n'.join(outlines))

    self.completed=True
     
####################################################
import subprocess as sub
import json
class DMCReader:
  ''' Reads results from a DMC calculation. 

  Attributes:
    output (dict): results of calculation. 
    completed (bool): whether the run has converged to a final answer.
  '''
  def __init__(self,errtol=0.01,minblocks=15):
    self.output={}
    self.completed=False

    self.errtol=errtol
    self.minblocks=minblocks
    self.gosling="gosling"

  def read_outputfile(self,outfile):
    ''' Read output file results.

    Args:
      outfile (str): output to read.
    '''
    return json.loads(sub.check_output([self.gosling,"-json",outfile.replace('.o','.log')]).decode())

  def check_complete(self):
    ''' Check if a DMC run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    completed=True
    if len(self.output)==0:
      return False # No results yet.
    if self.output['properties']['total_energy']['error'][0] > self.errtol:
      print("DMC incomplete: (%f) does not meet tolerance (%f)"%\
          (self.output['properties']['total_energy']['error'][0],self.errtol))
      completed=False
    if self.output['total blocks']-self.output['warmup blocks'] < self.minblocks:
      print("DMC incomplete: Run completed %d blocks, but requires %d."%\
          (self.output['total blocks']-self.output['warmup blocks'],self.minblocks))
      completed=False
    return completed
          
  #------------------------------------------------
  def collect(self,outfile):
    ''' Collect results for an output file and resolve if the run needs to be resumed. 

    Args: 
      outfiles (list): list of output file names to open and read.
    Returns:
      str: status of run = {'ok','restart'}
    '''
    # Gather output from files.
    self.completed=True
    status='unknown'
    if os.path.exists(outfile):
      self.output=self.read_outputfile(outfile)
      self.output['file']=outfile

    # Check files.
    self.completed=self.check_complete()
    if not self.completed:
      status='restart'
    else:
      status='ok'
    return status
      
  #------------------------------------------------
  def write_summary(self):
    ''' Print out all the items in output. '''
    print("#### Diffusion Monte Carlo")
    for f,out in self.output.items():
      print(f,out)
