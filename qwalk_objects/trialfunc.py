''' Trial function objects facilitate combining orbitals and other options to generating trial function input.'''
from numpy import array

#######################################################################
def export_qwalk_trialfunc(wf,**kwargs):
  ''' Wrapper for making making a wave function section into a trialfunc.
  Args:
    wf (TrialFunc): wave function for QWalk.
  Returns:
    str: trialfunc section.
  '''
  outlines=['trialfunc { ']+\
      ['  '+line for line in wf.export_qwalk_wf(**kwargs).split('\n')]+\
      ['}']
  return '\n'.join(outlines)

#######################################################################
class TrialFunc:
  ''' Abstract wave function class defined minimal interface.'''
  def export_qwalk_wf(self,**kwargs):
    ''' This must be implmented. '''
    raise NotImplementedError('Trial function not completely implemented')
  export_qwalk_trialfunc = export_qwalk_trialfunc

#################################################################################################
# Simple implementation for now. 
# In the future this should actually store the data for the wave function. 
class Jastrow(TrialFunc):
  ''' Class representing a Jastrow wave function. '''
  def __init__(self,jaststr):
    self.jaststr=jaststr

  def export_qwalk_wf(self,**kwargs):
    return self.jaststr

#################################################################################################
class Slater(TrialFunc):
  ''' Class representing a slater determinant wave function. '''
  def __init__(self,orbitals,orbfile,states,weights=(1.0,),shift_downorb=False):
    '''
    Args: 
      weights (array-like): Weights of determinants for multideterminant expansion. 
      states (array-like): states[determinant][spin channel][orbital] select orbitals for determinants.
        Indicies should reference whats written in the orbfile.
      orbfile (str): where orbitals are stored on disk (see write_qwalk_orb).
      orbitals (Orbitals): Something that can export_qwalk_orbitals(orbfile).
      shift_downorb (bool): Shift states[1] by number of up orbitals. 
        Useful for unrestricted calculations.
    '''
    self.weights=weights
    self.states=states
    self.orbfn=orbfile
    self.orbitals=orbitals
    self.shift_downorb=shift_downorb*(orbitals.eigvecs[0].shape[0])

  def export_qwalk_wf(self,optimize_det=False,rotate_orbs=None):
    ''' Write out a qwalk wave function section .
    Args:
      orbitals (Orbitals): object contraining orbitals. 
      orbfile (str): path to where orbitals are writting to disk.
      states (list): nested list of dimension 3. Indices are states[determinant][spin][orbital]. 
        Use 0-based indexing, as this will be converted to 1-based indexing in the export.
      optimize_det (bool): Add optimize_det flag to optimize determinant coefficients.
      rotate_orbs (list): Adds 'optimize_data' flag if not None.
        list of lists of orb groups to rotate between. For example, [[1,2,3],[4,5,6]]
    Returns:
      str: wave function section.
    '''
    # Convert python indexing to QWalk indexing.
    states = [[array(spinchannel) for spinchannel in det] for det in self.states]
    for det in states:
      for spinchan in det:
        assert (spinchan!=0).all(),"Are you using 0-based indexing? QWalk uses 1-based index!"
    weights=array(self.weights)

    # Check input validity.
    assert len(states[0])==2 and (len(states)==weights.shape[0]),\
        "States array should be nweights({}) by nspin(2) by nelectrons. One detweight per determinant.".format(weights.shape[0])

    for stdex in range(len(states)): states[stdex][1] += self.shift_downorb
    upstatelines = [' '.join((det[0]).astype(str)) for det in states]
    downstatelines = [' '.join((det[1]).astype(str)) for det in states]

    # Det. coefficients.
    if optimize_det: optimize_det_lines=['optimize_det']
    else:            optimize_det_lines=[]

    # Orbital rotations.
    if rotate_orbs is not None: 
      optimize_mo_lines=['optimize_mo','optimize_data { ','  det { ']
      for group in rotate_orbs:
        optimize_mo_lines += [
            '    orb_group { ',
            '      ' + ' '.join(map(str,group)),
            '    }'
          ]
      optimize_mo_lines += ['  }','}']
    else: optimize_mo_lines=[]

    outlines = [
        "slater",
        self.orbitals.export_qwalk_orbitals(self.orbfn),
        "detwt {{ {} }}".format(' '.join(weights.astype(str))),
        "states {"
      ]
    for detweight,upline,downline in zip(weights,upstatelines,downstatelines):
      outlines+=[
          "  # Spin up orbitals detweight {}.".format(detweight), 
          "  " + upline,
          "  # Spin down orbitals detweight {}.".format(detweight), 
          "  " + downline
        ]
    outlines+=['}']
    outlines += optimize_det_lines
    outlines += optimize_mo_lines
    return "\n".join(outlines)

#################################################################################################
class SlaterJastrow(TrialFunc):
  ''' Class rpresenting a slater determinant wave function. '''
  def __init__(self,slater,jastrow):
    '''
    Args: 
      slater (Slater): Slater wave function to multiply.
      jastrow (Jastrow): Jastrow wave function to multiply.
    '''
    self.slater=slater
    self.jastrow=jastrow
  def export_qwalk_wf(self,**kwargs):
    ''' Write out a qwalk wave function section .
    Returns:
      str: wave function section.
    '''
    outlines=[
        'slater-jastrow',
        '  wf1 {'
        ]+['    '+line for line in self.slater.export_qwalk_wf(**kwargs).split('\n')]+[
        '  }',
        '  wf2 {'
        ]+['    '+line for line in self.jastrow.export_qwalk_wf(**kwargs).split('\n')]+[
        '  }'
      ]
    return '\n'.join(outlines)

######################################################################
def separate_jastrow(wffile,optimizebasis=True,freezeall=False):
  ''' Seperate the jastrow section of a QWalk wave function file.
  Args:
    optimizebasis (bool): whether to leave optimizebasis keyword in the section. 
    freezeall (bool): whether to add freeze to all sections with parameters. 
  Returns:
    str: Jastrow section for qwalk.
  '''
  wff=open(wffile,'r')
  tokens=wff.read().split('\n')
  in_jastrow=False
  nopen=0
  nclose=0
  jastlines=[]
  for line in tokens:

    if 'jastrow2' in line.lower() or 'jastrow3' in line.lower():
      in_jastrow=True

    if in_jastrow:

      if not optimizebasis and 'optimizebasis' in line.lower():
        line=' '.join([word for word in line.split() if word.lower() != 'optimizebasis'])
        if line=='': continue

      if 'body' in line.lower() and freezeall:
        line=line+' FREEZE'

      nopen+=line.count("{")
      nclose+=line.count("}")

      if nopen<nclose:
        in_jastrow=False
      else:
        jastlines.append(line)

  return '\n'.join(jastlines)
