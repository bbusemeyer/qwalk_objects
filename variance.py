from __future__ import print_function
from system import find_cutoff_divider
from orbitals import find_min_exp
import os
####################################################
class VarianceWriter:
  def __init__(self,sys,trialfunc,iterations=10,macro_iterations=3):
    ''' Object for producing input into a variance optimization QWalk run. 
    Args:
      trialfunc (str): trial wavefunction section or object that can export_qwalk_trialfunc().
      sys (str): system section or object that can export_qwalk_sys().
      iterations (int): number of VMC steps to attempt.
      macro_iterations (int): Number of optimize calls to make.
    '''
    self.sys = sys
    self.trialfunc = trialfunc
    self.macro_iterations = macro_iterations
    self.iterations = iterations
    self.completed=False
    
  #-----------------------------------------------
  def qwalk_input(self,infile):
    ''' Generate QWalk input file and write to infile.'''
    assert self.trialfunc is not None, "Must specify trialfunc before asking for qwalk_input."
    assert self.sys is not None, "Must specify system before asking for qwalk_input."

    # Output all to strings.
    if type(self.sys) != str:
      sys = self.sys.export_qwalk_sys()
    else:
      sys = self.sys
    if type(self.trialfunc) != str:
      trialfunc=self.trialfunc.export_qwalk_trialfunc()
    else:
      trialfunc=self.trialfunc

    with open(infile,'w') as f:
      for j in range(self.macro_iterations):
        f.write("method { optimize iterations %i }\n"%self.iterations)
      f.write(sys),
      f.write(trialfunc)
    self.completed=True
     
####################################################
class VarianceReader:
  def __init__(self,vartol=10,vardifftol=0.1,minsteps=2):
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
    self.output={}
    self.completed=False

    self.vartol=vartol
    self.vardifftol=vardifftol
    self.minsteps=minsteps

  #------------------------------------------------
  def read_outputfile(self,outfile):
    ''' Read output file into dict.'''
    ret={}
    ret['sigma_trace']=[]
    with open(outfile,'r') as f:
      for line in f:
        if 'dispersion' in line:
          ret['sigma_trace'].append(float(line.split()[4]))
    if len(ret['sigma_trace'])>0:
      ret['sigma']=ret['sigma_trace'][-1]
    else:
      ret['sigma']=None
    return ret

  #------------------------------------------------
  def collect(self,outfile):
    ''' Collect results for output file and resolve if the run needs to be resumed based on accuracy requirements.

    Args: 
      outfile (str): Output file name to open and read.
    Returns:
      str: status of run = {'ok','restart'}
    '''
    # Gather output from files.
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
  def check_complete(self):
    ''' Check if a variance optimize run is complete.
    Returns:
      bool: If self.output are within error tolerances.
    '''
    if len(self.output['sigma_trace']) < self.minsteps:
      print(self.__class__.__name__,": Variance optimize incomplete: number of steps (%f) less than minimum (%f)"%\
          (len(self.output['sigma_trace']),self.minsteps))
      return False
    if self.output['sigma'] > self.vartol:
      print(self.__class__.__name__,": Variance optimize incomplete: variance ({}) does not meet tolerance ({})"\
          .format(self.output['sigma'],self.vartol))
      return False
    if (self.output['sigma_trace'][-1]-self.output['sigma_trace'][-2]) > self.vardifftol:
      print(self.__class__.__name__,": Variance optimize incomplete: change in variance (%f) less than tolerance (%f)"%\
          (self.output['sigma_trace'],self.vardifftol))
      return False
    return True
