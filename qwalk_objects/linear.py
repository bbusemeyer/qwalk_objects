from __future__ import print_function
import os
####################################################
class LinearWriter:
  def __init__(self,sys,trialfunc,total_nstep=2048*4,total_fit=2048):
    ''' Object for producing input into a variance optimization QWalk run. 
    Args:
      trialfunc (str): trial wavefunction section or object that can export_qwalk_trialfunc().
      sys (str): system section or object that can export_qwalk_sys().
      total_nstep (int): VMC steps to resolve Hamiltonian.
      total_fit (int): Steps to use when computing step size.
    '''
    self.sys=sys
    self.trialfunc=trialfunc
    self.total_nstep = total_nstep
    self.total_fit = total_fit
    self.completed=False

  #-----------------------------------------------
  def qwalk_input(self,infile):
    assert self.trialfunc is not None, "Must specify trialfunc before asking for qwalk_input."
    assert self.sys is not None, "Must specify system before asking for qwalk_input."

    # Output all to strings.
    if type(self.sys) != str:
      sys = self.sys.export_qwalk_sys()
    else:
      sys=self.sys
    if type(self.trialfunc) != str:
      trialfunc=self.trialfunc.export_qwalk_trialfunc()
    else:
      trialfunc=self.trialfunc

    with open(infile,'w') as f:
      f.write("method { linear \n")
      f.write("total_nstep %i \n"%self.total_nstep)
      f.write("total_fit %i \n"%self.total_fit)
      f.write("}\n")
      f.write(sys),
      f.write(trialfunc)
    self.completed=True

     
####################################################
class LinearReader:
  def __init__(self,sigtol=2.0,minsteps=2):
    ''' Object for reading, diagnosing, and storing Linear results.
    Args are only important for collect and check_complete.

    Args:
      sigtol (float): How many standard errors away from zero to you consider zero energy change?
      minsteps (int): Minimum number of steps allowed for no restart.
    Attributes:
      output (dict): Results for energy, error, and other information.
      completed (bool): Whether no more runs are needed.
    '''
    self.output={}
    self.completed=False
    self.sigtol=sigtol
    self.minsteps=minsteps

  #------------------------------------------------
  def read_outputfile(self,outfile):
    ret={}
    ret['energy_trace']=[]
    ret['energy_trace_err']=[]
    with open(outfile) as f:
      for line in f:
        if 'current energy' in line:
          ret['energy_trace'].append(float(line.split()[4]))
          ret['energy_trace_err'].append(float(line.split()[6]))
    if len(ret['energy_trace']) > 0:
      ret['total_energy']=ret['energy_trace'][-1]
      ret['total_energy_err']=ret['energy_trace_err'][-1]
    return ret

  #------------------------------------------------
  def check_complete(self):
    ''' Check if a variance optimize run is complete.
    Returns:
      bool: If self.output are within error tolerances.
    '''
    if len(self.output)==0:
      return False
    if len(self.output['energy_trace']) < self.minsteps:
      print(self.__class__.__name__,"Linear optimize incomplete: number of steps (%f) less than minimum (%f)"%\
          (len(self.output['energy_trace']),self.minsteps))
      return False
    else:
      ediff=self.output['energy_trace'][-1]-self.output['energy_trace'][-2]
      ediff_err=(self.output['energy_trace_err'][-1]**2 + self.output['energy_trace_err'][-2]**2)**0.5
      if ediff > self.sigtol*ediff_err:
        print(self.__class__.__name__,"Linear optimize incomplete: change in energy (%.5f) less than tolerance (%.2f*%.2f=%.5f)"%\
            (ediff,self.sigtol,ediff_err,self.sigtol*ediff_err))
        return False
    return True
          
  #------------------------------------------------
  def collect(self,outfile,errtol=None,minblocks=None):
    ''' Collect results for each output file and resolve if the run needs to be resumed. 

    Args: 
      outfile (str): output file to read.
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
  def write_summary(self):
    print("#### Linear optimization")
    for f,out in self.output.items():
      nruns=len(out)
      print(f,"Number of runs",nruns)
      for run in out:
        print("energy",run['energy'])
        print("energy_err",run['energy_err'])
      
      
      

