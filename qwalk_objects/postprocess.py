from __future__ import print_function
import os
import qwalk_objects as obj
import subprocess as sub
import json

####################################################
class PostprocessWriter:
  def __init__(self,sys,trialfunc,tracefn,energy=False,nskip=0,averages=''):
    ''' Object for producing input into a Postprocess QWalk run. 
    Args:
      trialfunc (str): trial wavefunction section or object that can export_qwalk_trialfunc().
      sys (str): system section or object that can export_qwalk_sys().
      tracefn (str): file name of trace to read in.
      energy (bool): Whether to compute the energy.
      nskip (int): Number of blocks to skip due to warmup.
      averages (str): averages section in qwalk.
    '''
    self.sys = sys
    self.trialfunc = trialfunc
    self.nskip = nskip
    self.tracefn = tracefn
    self.averages = averages
    self.energy = energy

    self.completed = False

  #-----------------------------------------------
  def qwalk_input(self,infile):
    ''' Generate QWalk input file and write to infile.'''
    assert self.trialfunc is not None, "Must specify trialfunc before asking for qwalk_input."
    assert self.sys is not None, "Must specify system before asking for qwalk_input."

    if type(self.sys) is not str:
      sys = self.sys.export_qwalk_sys()
    else:
      sys = self.sys
    if type(self.trialfunc) is not str:
      trialfunc = obj.trialfunc.export_qwalk_trialfunc(self.trialfunc)
    else:
      trialfunc = self.trialfunc

    outlines = [
        "method { postprocess",
        "  nskip %i"%self.nskip,
        "  readconfig %s"%self.tracefn,
    ]

    if not self.energy:
      outlines += ["  noenergy"]

    outlines += ['  '+line for line in self.averages.split('\n')]+[
        "}"
      ]
    outlines.append(sys)
    outlines.append(trialfunc)

    with open(infile,'w') as f:
      f.write('\n'.join(outlines))

    self.completed = True

####################################################
class PostprocessReader:
  def __init__(self,errtol=0.01,minblocks=15):
    ''' Object for reading and storing variance optimizer results.
    Args are only important for collect and check_complete.

    Args:
      vartol (float): Tolerance on the variance for collec.
      vardifftol (float): Tolerance of the change between the first and last variance.
      minsteps (int): minimun number of steps to attempt >= 2.
    Attributes:
      output (dict): Results for energy, error, and other information.
      completed (bool): Whether no more runs are needed.
    '''
    self.output = {}
    self.completed = False

    self.errtol = errtol
    self.minblocks = minblocks
    self.gosling = "gosling"

  def read_outputfile(self,jsonfn):
    ''' Read output file results.

    Args:
      outfile (str): output to read.
    '''
    return json.load(open(jsonfn,'r'))

  def check_complete(self):
    ''' Check if a postprocess run is complete.
    Returns:
      bool: If self.results are within error tolerances.
    '''
    completed = True
    if len(self.output)==0:
      return False # No results yet.
    return completed
          
  #------------------------------------------------
  def collect(self,outfile):
    ''' Collect results for an output file and resolve if the run needs to be resumed. 

    Args: 
      outfiles (list): list of output file names to open and read.
    Returns:
      str: status of run = {'ok','killed'}
    '''
    # Gather output from files.
    self.completed = True
    status = 'unknown'
    jsonfn = outfile.replace('.o','.json')
    if os.path.exists(jsonfn):
      self.output = self.read_outputfile(jsonfn)
      self.output['file'] = jsonfn

    # Check files.
    self.completed = self.check_complete()
    if not self.completed:
      status = 'killed'
    else:
      status = 'ok'
    return status
      
  #------------------------------------------------
  def write_summary(self):
    ''' Print out all the items in output. '''
    print("#### Variational Monte Carlo")
    for f,out in self.output.items():
      print(f,out)
